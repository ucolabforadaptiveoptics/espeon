from importlib import import_module

from .simulation import PhotonicLanternOptics
from .utils import *


_BACKEND_EXPORTS = {
    "save_lantern_design": (
        ".design", "lightbeam", "espeon[lightbeam]"
    ),
    "save_lantern_design_cbeam": (
        ".design_cbeam", "cbeam", "espeon[cbeam]"
    ),
}


def __getattr__(name):
    if name not in _BACKEND_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, backend, extra = _BACKEND_EXPORTS[name]
    try:
        value = getattr(import_module(module_name, __name__), name)
    except ModuleNotFoundError as error:
        raise ImportError(
            f"{name} requires the optional {backend} backend. "
            f"Install it with `pip install '{extra}'`."
        ) from error
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_BACKEND_EXPORTS))
