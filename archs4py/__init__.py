import archs4py.data
import archs4py.download
import archs4py.meta
import archs4py.utils

import importlib
importlib.reload(archs4py.data)
importlib.reload(archs4py.download)
importlib.reload(archs4py.meta)
importlib.reload(archs4py.utils)

from archs4py.utils import versions
from archs4py.utils import normalize


