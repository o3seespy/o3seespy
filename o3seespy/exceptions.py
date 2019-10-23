import warnings


class ModelError(Exception):
    pass


class ModelWarning(Warning):
    pass


def deprecation(message):
    warnings.warn(message, stacklevel=3)
