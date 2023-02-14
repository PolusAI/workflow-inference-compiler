from pathlib import Path
import re

import requests
from bs4 import BeautifulSoup

from ..wic_types import Json


def gromacs_selection_groups() -> Json:
    """A dictionary of all the gromacs selection group names, with their documentation.

    Returns:
        Json: A dictionary of all the gromacs selection group names, with their documentation.
    """
    groups_dict = {
        'System': 'all atoms in the system',
        'Protein': 'all protein atoms',
        'Protein-H': 'protein atoms excluding hydrogens',
        'C-alpha': 'C-alpha atoms',
        'Backbone': 'protein backbone atoms: N; C-alpha and C',
        'MainChain': 'protein main chain atoms: N; C-alpha; C and O; including oxygens in C-terminus',
        'MainChain+Cb': 'protein main chain atoms including C-beta',
        'MainChain+H': 'protein main chain atoms including backbone amide hydrogens and hydrogens on the N-terminus',
        'SideChain': """protein side chain atoms: that is all atoms except N; C-alpha; C; O;
        backbone amide hydrogens and oxygens in C-terminus and hydrogens on the N-terminus""",
        'SideChain-H': 'protein side chain atoms excluding all hydrogens',
        'Prot-Masses': 'protein atoms excluding dummy masses',
        'non-Protein': 'all non-protein atoms',
        'Water': 'water molecules',
        'SOL': 'water molecules',
        'non-Water': 'anything not covered by the Water group',
        'Ion': 'any name matching an Ion entry in residuetypes.dat',
        'NA': 'all NA atoms',
        'CL': 'all CL atoms',
        'Water_and_ions': 'combination of the Water and Ions groups',
        'DNA': 'all DNA atoms',
        'RNA': 'all RNA atoms',
        'Protein_DNA': 'all Protein-DNA complex atoms',
        'Protein_RNA': 'all Protein-RNA complex atoms',
        'Protein_DNA_RNA': 'all Protein-DNA-RNA complex atoms',
        'DNA_RNA': 'all DNA-RNA complex atoms'
    }
    return groups_dict


def gromacs_mdp_schema() -> Json:
    """The schema of all gromacs mdp options screen-scraped from\n
    https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html

    Returns:
        Json: A formal schema of the gromacs mdp options, including documentation.
    """
    # TODO: Gromacs just updated their mdp options html page for gromacs 2023.
    # The new html tags are slightly different, so this screen scraper needs updated.
    # For now, we can store the old version of the html in this repository.
    gromacs_mdp_html_file = Path(__file__).parent / 'gromacs_mdp.html'
    if gromacs_mdp_html_file.exists():
        with open(gromacs_mdp_html_file, mode='r', encoding='utf-8') as f:
            html_content = f.read()
    else:
        mdp_url = 'https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html'
        response = requests.get(mdp_url, timeout=60) # seconds
        html_content = response.text # cache this locally
        with open(gromacs_mdp_html_file, mode='w', encoding='utf-8') as f:
            f.write(html_content)

    def gettype(key: str, desc: str) -> str:
        # Use heuristic to determine most number types
        # NOTE: This does not attempt to determine lists/arrays of numbers.

        number_in_parens_pattern = '.*\\(-?[0-9]+\\.?[0-9]*\\).*'
        m = re.match(number_in_parens_pattern, desc)

        # TODO: Triple-check that all of the mdp keys which ought to be numeric
        # are actually numeric. If not, add them to this explicit list.
        # NOTE: Some of these entries have default values in scientific notation,
        # so it is unclear if they should be encoded as numbers or strings.
        explicit = ['fourier-nx', 'fourier-ny', 'ewald-rtol', 'ewald-rtol-lj',
                    'ref-t', 'tau-t', 'compressibility', 'ref-p', 'wall-density',
                    'pull-constr-tol', 'pull-coord1-dx', 'pull-coord1-kB',
                    'awh1-dim1-diffusion', 'nstexpanded', 'mc-temperature',
                    'density-guided-simulation-force-constant',
                    'qmmm-cp2k-qmcharge', 'qmmm-cp2k-qmmultiplicity']

        if m or key in explicit:
            return 'number'
        return 'string'

    mdp: Json = {}

    # Look for the root-level tags
    # <dl class="std mdp">
    soup = BeautifulSoup(html_content, "html.parser")
    resultset = soup.find_all(name='dl', attrs={'class': 'std mdp'})
    for result in resultset:
        #print(result.prettify())

        prefix = 'mdp-'
        name = result.dt['id'][len(prefix):] # remove prefix
        desc = result.dd

        # Escape html tags in description with double quotes for json serialization.
        # NOTE: This may be related to the following issue:
        # https://github.com/redhat-developer/vscode-yaml/issues/381
        mdp[name] = {'type': gettype(name, str(desc)), 'description': f'"{desc}"'} # placeholder schema

        # Look for sub-tags / enumerated values
        # <dl class="std mdp-value">
        subsoup = BeautifulSoup(str(result), "html.parser")
        subresultset = subsoup.find_all(name='dl', attrs={'class': 'std mdp-value'})

        values_schemas = []
        for subresult in subresultset:
            #print(subresult.prettify())

            prefix = f'mdp-value-{name}-'
            valuename = subresult.dt['id'][len(prefix):] # remove prefix
            valuedesc = subresult.dd

            # It looks like all the numeric types are at the root level.
            # i.e. we probably don't need gettype() here.
            schema = {'type': gettype(valuename, str(valuedesc)), # 'title': '',
                      'const': valuename, 'description': f'"{valuedesc}"'}
            values_schemas.append(schema)

        # Escape html tags in description with double quotes for json serialization.
        # NOTE: This may be related to the following issue:
        # https://github.com/redhat-developer/vscode-yaml/issues/381
        if len(values_schemas) > 0:
            # NOTE: Use anyOf instead of enum so we can add descriptions to each value
            mdp[name] = {'anyOf': values_schemas, 'description': f'"{desc}"'}
            #mdp[name] = {'type': 'string', 'enum': [s['const'] for s in values_schemas]}

        # Overwrite some cases that are handled incorrectly above
        if name == 'nstlist':
            mdp[name] = {'type': 'number', 'description': f'"{desc}"'}
        if name == 'rot-fit-method0':
            mdp[name] = {'type': 'string', 'enum': ['rmsd', 'norm', 'potential']}
        if 'userint' in name or 'userreal' in name:
            mdp[name] = {'type': 'number', 'description': f'"{desc}"'}

    return mdp
