from importlib.metadata import version as _version, PackageNotFoundError
from .hmf2smf import *

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass