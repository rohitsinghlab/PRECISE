# global_vars.py: Global variables used by MaSIF -- mainly pointing to environment variables of programs used by MaSIF.
# Pablo Gainza - LPDI STI EPFL 2018-2019
# Released under an Apache License 2.0

epsilon = 1.0e-6

# These are now handled by verify_executables() in msms_pipeline.py
# and accessed directly via os.environ in other modules.
msms_bin = None
pdb2pqr_bin = None
apbs_bin = None
multivalue_bin = None


class NoSolutionError(Exception):
    pass
