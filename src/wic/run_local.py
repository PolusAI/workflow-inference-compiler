import argparse
import glob
import json
import subprocess as sub
import sys
from pathlib import Path
import shutil
import platform
import traceback
from typing import Dict, List, Optional

try:
    import cwltool.main
except ImportError as exc:
    print('Could not import cwltool.main')
    # (pwd is imported transitively in cwltool.provenance)
    print(exc)
    if exc.msg == "No module named 'pwd'":
        print('Windows does not have a pwd module')
        print('If you want to run on windows, you need to install')
        print('Windows Subsystem for Linux')
        print('See https://pypi.org/project/cwltool/#ms-windows-users')
    else:
        raise exc

from . import utils  # , utils_graphs
from .wic_types import Yaml, RoseTree


def run_local(args: argparse.Namespace, rose_tree: RoseTree, cachedir: Optional[str], cwl_runner: str, use_subprocess: bool) -> int:
    """This function runs the compiled workflow locally.

    Args:
        args (argparse.Namespace): The command line arguments
        rose_tree (RoseTree): The compiled workflow
        cachedir (Optional[str]): The --cachedir to use (if any)
        cwl_runner (str): Either 'cwltool' or 'toil-cwl-runner'
        use_subprocess (bool): When using cwltool, determines whether to use subprocess.run(...) or use the cwltool python api.

    Returns:
        retval: The return value
    """

    docker_cmd: str = args.user_space_docker_cmd
    # Check that docker is installed, so users don't get a nasty runtime error.
    cmd = [docker_cmd, 'run', 'hello-world']
    try:
        docker_cmd_exists = True
        proc = sub.run(cmd, check=False, stdout=sub.PIPE, stderr=sub.STDOUT)
        output = proc.stdout.decode("utf-8")
    except FileNotFoundError:
        docker_cmd_exists = False
    out_d = "Hello from Docker!"
    out_p = "Hello Podman World"
    if (not docker_cmd_exists or not (proc.returncode == 0 and out_d in output or out_p in output)) and not args.ignore_docker_install:
        print(f'Warning! {docker_cmd} does not appear to be installed.')
        print(f'Most workflows require containers and will fail at runtime if {docker_cmd} is not installed.')
        print('If you want to run the workflow anyway, use --ignore_docker_install')
        sys.exit(1)

    # If docker is installed, check for too many running processes. (on linux, macos)
    if docker_cmd == 'docker' and docker_cmd_exists and sys.platform != "win32":
        cmd = 'pgrep com.docker | wc -l'  # type: ignore
        proc = sub.run(cmd, check=False, stdout=sub.PIPE, stderr=sub.STDOUT, shell=True)
        output = proc.stdout.decode("utf-8")
        num_processes = int(output.strip())
        max_processes = 1000
        if num_processes > max_processes and not args.ignore_docker_processes:
            print(f'Warning! There are {num_processes} running docker processes.')
            print(f'More than {max_processes} may potentially cause intermittent hanging issues.')
            print('It is recommended to terminate the processes using the command')
            print('`sudo pkill com.docker && sudo pkill Docker`')
            print('and then restart Docker.')
            print('If you want to run the workflow anyway, use --ignore_docker_processes')
            sys.exit(1)

    yaml_path = args.yaml
    yaml_stem = Path(args.yaml).stem

    yaml_inputs = rose_tree.data.workflow_inputs_file
    stage_input_files(yaml_inputs, Path(args.yaml).parent.absolute())

    retval = 1  # overwrite if successful
    provenance: List[str] = []
    # NOTE: By default, cwltool will attempt to download schema files.
    # $schemas:
    #   - https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
    # If you have connection issues (e.g. firewall, VPN, etc) then failure to download will
    # not actually cause any problems immediately (other than a ~30 second timeout).
    # However, cwltool does not appear to cache these files, so it will attempt to download
    # them repeatedly. These ~30 second timeouts will eventually add up to >6 hours, which
    # will cause github to terminate the CI Action...
    skip_schemas = ['--skip-schemas'] if not args.no_skip_dollar_schemas else []

    yaml_stem = yaml_stem + '_inline' if args.cwl_inline_subworkflows else yaml_stem
    if cwl_runner == 'cwltool':
        parallel = ['--parallel'] if args.parallel else []
        # NOTE: --parallel is required for real-time analysis / real-time plots,
        # but it seems to cause hanging with Docker for Mac. The hanging seems
        # to be worse when using parallel scattering.
        quiet = ['--quiet'] if args.quiet else []
        cachedir_ = ['--cachedir', cachedir] if cachedir else []
        net = ['--custom-net', args.custom_net] if args.custom_net else []
        provenance = ['--provenance', 'provenance']
        docker_cmd_ = [] if docker_cmd == 'docker' else ['--user-space-docker-cmd', docker_cmd]
        write_summary = ['--write-summary', args.write_summary] if args.write_summary else []
        path_check = ['--relax-path-checks']
        # See https://github.com/common-workflow-language/cwltool/blob/5a645dfd4b00e0a704b928cc0bae135b0591cc1a/cwltool/command_line_tool.py#L94
        # NOTE: Using --leave-outputs to disable --outdir
        # See https://github.com/dnanexus/dx-cwl/issues/20
        # --outdir has one or more bugs which will cause workflows to fail!!!
        cmd = ['cwltool'] + parallel + quiet + cachedir_ + net + provenance + \
            docker_cmd_ + write_summary + skip_schemas + path_check
        cmd += ['--leave-outputs',
                # '--js-console', # "Running with support for javascript console in expressions (DO NOT USE IN PRODUCTION)"
                f'autogenerated/{yaml_stem}.cwl', f'autogenerated/{yaml_stem}_inputs.yml']
        # TODO: Consider using the undocumented flag --fast-parser for known-good workflows,
        # which was recently added in the 3.1.20220913185150 release of cwltool.

        print('Running ' + ' '.join(cmd))
        if use_subprocess:
            # To run in parallel (i.e. pytest ... --workers 8 ...), we need to
            # use separate processes. Otherwise:
            # "signal only works in main thread or with __pypy__.thread.enable_signals()"
            proc = sub.run(cmd, check=False)
            retval = proc.returncode
            return retval  # Skip copying files to outdir/ for CI
        else:
            print('via python API')
            try:
                retval = cwltool.main.main(cmd[1:])
                assert retval == 0

                if args.write_summary:
                    print(f'Final output json blob is in {args.write_summary}')

                # See https://pypi.org/project/cwltool/#import-as-a-module
                # This also works, but doesn't easily allow using --leave-outputs, --provenence, --cachedir
                # import cwltool.factory
                # fac = cwltool.factory.Factory()
                # rootworkflow = fac.make(f'autogenerated/{yaml_stem}.cwl')
                # output_json = rootworkflow(**yaml_inputs)
                # with open('primary-output.json', mode='w', encoding='utf-8') as f:
                #     f.write(json.dumps(output_json))
            except Exception as e:
                print('Failed to execute', yaml_path)
                print(f'See error_{yaml_stem}.txt for detailed technical information.')
                # Do not display a nasty stack trace to the user; hide it in a file.
                with open(f'error_{yaml_stem}.txt', mode='w', encoding='utf-8') as f:
                    traceback.print_exception(etype=type(e), value=e, tb=None, file=f)
                if not cachedir:  # if running on CI
                    print(e)

    if cwl_runner == 'toil-cwl-runner':
        if platform.python_implementation().lower() == 'pypy':
            print('Error! Toil is not compatible with pypy!')
            print('Please use the standard cpython interpreter with Toil.')
            sys.exit(1)

        # NOTE: toil-cwl-runner always runs in parallel
        net = ['--custom-net', args.custom_net] if args.custom_net else []
        provenance = ['--provenance', 'provenance']
        docker_cmd_ = [] if docker_cmd == 'docker' else ['--user-space-docker-cmd', docker_cmd]
        cmd = ['toil-cwl-runner'] + net + provenance + docker_cmd_
        cmd += ['--outdir', 'outdir_toil',
                '--jobStore', f'file:./jobStore_{yaml_stem}',  # NOTE: This is the equivalent of --cachedir
                # TODO: Check --clean, --cleanWorkDir, --restart
                '--clean', 'always',  # This effectively disables caching, but is reproducible
                f'autogenerated/{yaml_stem}.cwl', f'autogenerated/{yaml_stem}_inputs.yml']

        print('Running ' + ' '.join(cmd))
        proc = sub.run(cmd, check=False)
        retval = proc.returncode

    if retval == 0:
        print('Success! Output files should be in outdir/')
    else:
        print('Failure! Please scroll up and find the FIRST error message.')
        print('(You may have to scroll up A LOT.)')

    # Remove the annoying cachedir* directories that somehow aren't getting automatically deleted.
    # NOTE: Do NOT allow cachedir to be absolute; otherwise
    # if users pass in "/" this will delete their entire hard drive.
    cachedir_path = str(cachedir)
    if not Path(cachedir_path).is_absolute():
        for d in glob.glob(cachedir_path + '*'):
            if not d == cachedir_path:
                shutil.rmtree(d)  # Be VERY careful when programmatically deleting directories!

    # Finally, since there is an output file copying bug in cwltool,
    # we need to copy the output files manually. See comment above.
    output_json_file_prov = Path('provenance/workflow/primary-output.json')
    # NOTE: The contents of args.write_summary (as printed to stdout) is
    # slightly different than provenance/workflow/primary-output.json!
    # They are NOT the same file!
    if output_json_file_prov.exists():
        with open(output_json_file_prov, mode='r', encoding='utf-8') as f:
            output_json = json.loads(f.read())
        files = utils.parse_provenance_output_files(output_json)

        dests = set()
        for location, namespaced_output_name, basename in files:
            yaml_stem_init, shortened = utils.shorten_namespaced_output_name(namespaced_output_name)
            parentdirs = yaml_stem_init + '/' + shortened.replace('___', '/')
            Path('outdir/' + parentdirs).mkdir(parents=True, exist_ok=True)
            source = 'provenance/workflow/' + location
            # NOTE: Even though we are using subdirectories (not just a single output directory),
            # there is still the possibility of filename collisions, i.e. when scattering.
            # For now, let's use a similar trick as cwltool of append _2, _3 etc.
            # except do it BEFORE the extension.
            # This could still cause problems with slicing, i.e. if you scatter across
            # indices 11-20 first, then 1-10 second, the output file indices will get switched.
            dest = 'outdir/' + parentdirs + '/' + basename
            if dest in dests:
                idx = 2
                while Path(dest).exists():
                    stem = Path(basename).stem
                    suffix = Path(basename).suffix
                    dest = 'outdir/' + parentdirs + '/' + stem + f'_{idx}' + suffix
                    idx += 1
            dests.add(dest)
            cmd = ['cp', source, dest]
            sub.run(cmd, check=True)

    return retval


def stage_input_files(yml_inputs: Yaml, root_yml_dir_abs: Path,
                      relative_run_path: bool = True, throw: bool = True) -> None:
    """Copies the input files in yml_inputs to the working directory.

    Args:
        yml_inputs (Yaml): The yml inputs file for the root workflow.
        root_yml_dir_abs (Path): The absolute path of the root workflow yml file.
        relative_run_path (bool): Controls whether to use subdirectories or\n
        just one directory when writing the compiled CWL files to disk
        throw (bool): Controls whether to raise/throw a FileNotFoundError.

    Raises:
        FileNotFoundError: If throw and it any of the input files do not exist.
    """
    for key, val in yml_inputs.items():
        if isinstance(val, Dict) and val.get('class', '') == 'File':
            path = root_yml_dir_abs / Path(val['path'])
            if not path.exists() and throw:
                # raise FileNotFoundError(f'Error! {path} does not exist!')
                print(f'Error! {path} does not exist!')
                print('(Did you forget to use an explicit edge?)')
                print('See https://workflow-inference-compiler.readthedocs.io/en/latest/userguide.html#explicit-edges')
                sys.exit(1)

            relpath = Path('autogenerated/') if relative_run_path else Path('.')
            pathauto = relpath / Path(val['path'])  # .name # NOTE: Use .name ?
            pathauto.parent.mkdir(parents=True, exist_ok=True)

            if path != pathauto:
                cmd = ['cp', str(path), str(pathauto)]
                proc = sub.run(cmd, check=False)
