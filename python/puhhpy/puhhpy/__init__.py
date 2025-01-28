# Expose some main function as `puhhpy.some_function()` instead of having to use `puhhpy.subdir.some_function()`.
#from .main.initial_fields import create_initial_fields

# Expose sub-directories as `import puhhpy; puhhpy.subdir.some_function()`
# NOTE: this only exposes what is defined in the subdirectory `__init__.py`.
from .spatial import *
from .thermo import *
from .interpolate import *