from .version import __version__
from .system import System

# if somebody does "from somepackage import *", this is what they will
# be able to access:
__all__ = [
    'System'
]
