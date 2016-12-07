#! /usr/bin/env python


"""
Pipeline.

porting pipeline to python
"""
###############################################################################
import sys
###############################################################################

###############################################################################
# PACKAGE METADATA
import collections
version_info = collections.namedtuple("pypiret_version_info",
                                      ["major", "minor", "micro",
                                       "releaselevel"])(major=1, minor=0,
                                                        micro=0)
__project__ = "PyPiReT"
__version__ = ".".join(str(s) for s in version_info[:4] if s != "")


def name():
    """Name.

    from dendropy
    """
    return "{} {}{}".format(__project__, __version__, revision_description())


def _get_revision_object():
    from dendropy.utility import vcsinfo
    __revision__ = vcsinfo.Revision(repo_path=homedir())
    return __revision__


def revision_description():
    """
    Revision description.

    revision
    """
    __revision__ = _get_revision_object()
    if __revision__.is_available:
        revision_text = " ({})".format(__revision__)
    else:
        revision_text = ""
    return revision_text


def homedir():
    """
    Get home directory.

    get homedir
    """
    import os
    try:
        try:
            __homedir__ = __path__[0]
        except AttributeError:
            __homedir__ = os.path.dirname(os.path.abspath(__file__))
        except IndexError:
            __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except OSError:
        __homedir__ = None
    except:
        __homedir__ = None
    return __homedir__


def description(dest=None):
    """
    Contain descripiton.

    gies out different attributes of the package
    """
    import sys
    import site
    if dest is None:
        dest = sys.stdout
    fields = collections.OrderedDict()
    fields["PyPiReT version"] = name()
    fields["PyPiReT location"] = homedir()
    fields["Python version"] = sys.version.replace("\n", "")
    fields["Python executable"] = sys.executable
    try:
        fields["Python site packages"] = site.getsitepackages()
    except:
        pass
    max_fieldname_len = max(len(fieldname) for fieldname in fields)
    for fieldname, fieldvalue in fields.items():
        dest.write("{fieldname:{fieldnamewidth}}: {fieldvalue}\n".format(
            fieldname=fieldname,
            fieldnamewidth=max_fieldname_len + 2,
            fieldvalue=fieldvalue))

if __name__ == "__main__":
    description(sys.stdout)
