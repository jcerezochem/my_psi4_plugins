#
#@BEGIN LICENSE
#
# scf_pc by Psi4 Developer, a plugin to:
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def run_scf_pc(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    scf_pc can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('scf_pc')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    #psi4.set_global_option('BASIS', 'sto-3g')
    psi4.set_local_option('SCF_PC', 'PRINT', 1)
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('scf_pc.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)


# Integration with driver routines
procedures['energy']['scf_pc'] = run_scf_pc


def exampleFN():
    # Your Python code goes here
    pass
