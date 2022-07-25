from ._version import __version__

from . import _muscle
from ._muscle import (
    Sequence,
    MultiSequence,
    Alignment,
    _AlignmentSequences,
    Aligner,
)

__doc__ = _muscle.__doc__
__all__ = [
    "Sequence",
    "MultiSequence",
    "Alignment",
    "_AlignmentSequences",
    "Aligner",
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pymuscle5.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
