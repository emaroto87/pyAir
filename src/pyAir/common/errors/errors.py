# -*- coding: utf-8 -*-
from packaging import version
import importlib.metadata
'''
BASIC ERROR CHECKS
'''
__version__ = '1.0.0'
__author__ = 'E.Maroto'
__license__ = 'MIT'


def check_type(obj, obj_class):
    if not (isinstance(obj, obj_class)):
        error_msg = 'Object {0} is not a {1}type'.format(
            str(type(obj)),
            str(type(obj_class))
        )
        raise TypeError(error_msg)


def check_option(option: str, options: list[str]):
    check_type(option, str)
    check_type(options, list)
    if not (option in options):
        raise ValueError(
            f'{option} is not one of the valid values included in {options}'
        )


def typerror_msg(obj: object, expected_obj: object) -> str:
    msg = '{0} object is not a {1}'.format(
        type(obj),
        expected_obj
    )
    return TypeError(msg)


class error_messages:
    def __init__(self):
        self.invalid_option = '[Error] Incorret option. Please try again.'


class ModuleVersionError(ImportError):
    """
    Raise when the version of the imported module does not match the requirements.
    """
    pass


def check_module_version(module, min_version: str):
    """
    Check that the installed version of the module is equal or greater to the 
    minimum allowable version.
    """

    try:
        installed_version = module.__version__
    except AttributeError:
        raise ModuleVersionError(
            f"Module {module.__name.__} does not has attribute __version__"
        )

    if version.parse(installed_version) < version.parse(min_version):
        raise ModuleVersionError(
            f"It is required that {module.__name__} version >= {min_version}"
            f"but found {installed_version}"
        )
    else:
        return True
