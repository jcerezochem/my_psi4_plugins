#
#@BEGIN LICENSE
#
# vibronic by Psi4 Developer, a plugin to:
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

def run_vibronic(state1,state2,method1='ccsd',method2='eom-ccsd'):
    """
    Performs the different steps needed to compute the spectrum
    """

    lowermethod1 = method1.lower()
    lowermethod2 = method2.lower()

    # Perform optimizations and save states 
    psi4.set_active_molecule(state1)
    print "Optimizing Ground State..."
    optimize(lowermethod1)
    E1=psi4.get_variable("CURRENT ENERGY")

    psi4.set_active_molecule(state2)
    print "Optimizing Excited State..."
    optimize(lowermethod2)
    E2=psi4.get_variable("CURRENT ENERGY")
    
    # Compute Hessian on first state (for AS model)
    psi4.set_active_molecule(state1)
    frequencies(lowermethod1)
    
    # Select states for vibronic:
    #  active molecule is initial state
    #  secondary molecule is final state
    psi4.set_active_molecule(state1)
    psi4.set_secondary_molecule(state2)
    
    # Compute adiavatic energy
    DE = (E2-E1)
    psi4.set_variable("CURRENT ENERGY",DE)

    # Call the plugin
    psi4.plugin('vibronic.so')

def spectrum():
    # Simple call to the plugin
    psi4.plugin('vibronic.so')


