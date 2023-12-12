import archs4py.data
import archs4py.download
import archs4py.meta
import archs4py.utils
import archs4py.align

import importlib
importlib.reload(archs4py.data)
importlib.reload(archs4py.download)
importlib.reload(archs4py.meta)
importlib.reload(archs4py.utils)
importlib.reload(archs4py.align)

from archs4py.utils import versions
from archs4py.utils import normalize
from archs4py.utils import ls

__version__="0.2.12"
