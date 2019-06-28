"""
Some useful tools for writing relevant validation notebooks.
"""
from importlib import import_module


def print_dep_versions(extras=None):
    """
    Prints versions of all important "active" modules.

    This includes modules that are not explicitly imported, as they *may* be used as
    deps of used packages. It will skip any module that isn't installed at all (since
    obviously this is not being used).

    :param extras: Any extra modules that may be useful for this particular notebook.
    """
    MODULES = [
        "pyuvdata",
        "hera_stats",
        "hera_sim",
        "hera_qm",
        "hera_pspec",
        "linsolve",
        "uvtools",
        "numpy",
        "healvis",
        "healpy",
        "h5py"
    ]

    if extras is not None:
        MODULES += extras

    for module in MODULES:
        try:
            _mdl = import_module(module)
        except ModuleNotFoundError:
            pass

        if hasattr(_mdl, 'version'):
            gh = getattr(_mdl.version, 'git_hash', None)

        print("Module {:<11}....\tVersion {:<7}.......\tGit {}".format(module, _mdl.__version__, gh))
