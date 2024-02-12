import json
import subprocess as sub

# https://cwl.discourse.group/t/override-docker-entrypoint-in-command-line-tool/695/2
# "Here is what the CWL standards have to say about software container entrypoints"
# "I recommend changing your docker container to not use ENTRYPOINT."

# This script will remove / overwrite the entrypoints in ALL of the polusai/ docker images on your local machine.

# The reason is that we want to simultaneously allow
# 1. release environments, where users will simply run code in a docker image and
# 2. dev/test environments, where developers can run the latest code on the host machine and/or in the CI.

# With entrypoints, there is no easy way to switch between 1 and 2; developers will
# have to manually prepend the entrypoint string to the baseCommand, and/or possibly
# modify paths to be w.r.t. their host machine instead of w.r.t. the image. (/opt/.../main.py)

# Without entrypoints, to switch between 1 and 2, simply comment out DockerRequirement ... that's it!
# The point is that we want a uniform API so that we can programmatically switch between 1 and 2 in the CI.
# We want to run the integration tests first, and only push releases to dockerhub when the tests pass (NOT vice versa!).

docker_image_cmd = ['docker', 'image', 'list', '--no-trunc', '--format', 'json']

proc = sub.run(docker_image_cmd, check=True, stdout=sub.PIPE)
output = proc.stdout.decode("utf-8")
lines = output.splitlines()

for json_str in lines:
    # print('json_str', json_str)
    image_json = json.loads(json_str)
    # print(image_json)
    repo = image_json['Repository']
    tag = image_json['Tag']
    if repo != '<none>' and tag != '<none>' and repo.startswith('polusai/'):
        print(f'{repo}:{tag}')
        with open('Dockerfile_tmp', mode='w', encoding='utf-8') as f:
            f.write(f'FROM {repo}:{tag}')
            f.write('\n')
            f.write('ENTRYPOINT []')
        docker_build_cmd = ['sudo', 'docker', 'build', '-f', 'Dockerfile_tmp', '-t', f'{repo}:{tag}', '.']
        sub.run(docker_build_cmd, check=True)
        sub.run(['rm', 'Dockerfile_tmp'], check=True)
