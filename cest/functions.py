import functools
import inspect
import sys

import spire

from . import tasks

# Wrap all tasks as functions. The resulting function will
# - be named after the class module: cest.wassr.WASSR → wassr
# - carry the docstring of the class
# - correctly display the arguments when documented (based on the `wrapper`
#   decorator)

def wrapper(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    return wrapper

for member in tasks.__dict__.values():
    if inspect.isclass(member) and issubclass(member, spire.TaskFactory):
        name = member.__module__.rsplit(".", 1)[-1]
        wrapped = wrapper(member.action)
        wrapped.__name__ = name
        wrapped.__doc__ = member.__doc__
        setattr(sys.modules[__name__], name, wrapped)
