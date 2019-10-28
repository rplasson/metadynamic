from .version import __version__
from .polym import System

# if somebody does "from somepackage import *", this is what they will
# be able to access:
__all__ = [
    'shout_and_repeat',
    'my_add',
]
