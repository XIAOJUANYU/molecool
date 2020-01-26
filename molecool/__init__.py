"""
molecool
a python package for analyzing and virualizing xyz file for molssi best practices workshop
"""

# Add imports here
from .functions import *
from .measure import calculate_angle, calculate_distance
from .visulize import draw_molecule, bond_histogram
from .molecule import build_bond_list, calculate_molecular_mass
from .atom_data import atomic_colors, atomic_weights

import molecool.io

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
