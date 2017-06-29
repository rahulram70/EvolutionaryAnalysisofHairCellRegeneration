
import pkgutil

__name__ = "Procomp Master Package"
__doc__ = """ here is the master doc """

__all__ = []
for loader, module_name, is_pkg in  pkgutil.walk_packages(__path__):
    __all__.append(module_name)
    module = loader.find_module(module_name).load_module(module_name)
    exec('%s = module' % module_name)

