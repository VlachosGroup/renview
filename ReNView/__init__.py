# -*- coding: utf-8 -*-
"""
ReNView
"""

####
#
# setuptools likes to see a name for the package,
# and it's best-practices to have the __version__
# present, too:
#
name = 'ReNView'
__version__ = '1.1'

import GraphGenerator
import example
import Legend
import Reactions

print(__version__)