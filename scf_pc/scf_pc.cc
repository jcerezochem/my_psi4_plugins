/*
 *@BEGIN LICENSE
 *
 * scf_pc by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <boost/python.hpp>

#include <cstdlib>
#include <cstdio>
#include <string>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libmints/writer.h>
#include <libmints/writer_file_prefix.h>

#include <libscf_solver/rhf.h>
#include <libscf_solver/rohf.h>
#include <libscf_solver/uhf.h>
#include <libscf_solver/cuhf.h>
#include <libscf_solver/ks.h>

#include "rhf_pc.h"
#include "PointCharges.h"

INIT_PLUGIN

using namespace boost;
using namespace std;

namespace psi{ namespace scf {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "SCF_PC"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- Type of reference function -*/
        options.add_str("REFERENCE", "RHF");
        /*- Print basis -*/
        options.add_bool("PRINT_BASIS", false);
        /*- Print MOLDEN -*/
        options.add_bool("MOLDEN_WRITE", false);
        /*- Print MOs -*/
        options.add_bool("PRINT_MOS", false);
        /*- External file with point charges -*/
        options.add_str_i("POINTCHARGES_FILE", "pointcharges.dat");
    }

    return true;
}

extern "C" 
PsiReturnType scf_pc(Options& options, PyObject* pre, PyObject* post)
{
    tstart();

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();

    string reference = options.get_str("REFERENCE");
    boost::shared_ptr<Wavefunction> scf_pc;
    double energy;
    //bool parallel = options.get_bool("PARALLEL");


    if (reference == "RHF") {
        scf_pc = boost::shared_ptr<Wavefunction>(new scf::RHF_PC(options, psio));
    }
//    else if (reference == "ROHF") {
//        scf_pc = boost::shared_ptr<Wavefunction>(new scf::ROHF(options, psio));
//    }
//    else if (reference == "UHF") {
//        scf_pc = boost::shared_ptr<Wavefunction>(new scf::UHF(options, psio));
//    }
//    else if (reference == "CUHF") {
//        scf_pc = boost::shared_ptr<Wavefunction>(new scf::CUHF(options, psio));
//    }
//    else if (reference == "RKS") {
//        scf_pc = boost::shared_ptr<Wavefunction>(new scf::RKS(options, psio));
//    }
//    else if (reference == "UKS") {
//        scf_pc = boost::shared_ptr<Wavefunction>(new scf::UKS(options, psio));
//    }
    else {
        throw InputException("Unknown reference for PointCharges " + reference, "REFERENCE", __FILE__, __LINE__);
        energy = 0.0;
    }

    // print the basis set
    if ( options.get_bool("PRINT_BASIS") ) {
        boost::shared_ptr<BasisSet> basisset = BasisSet::pyconstruct_orbital(Process::environment.molecule(),
            "BASIS", options.get_str("BASIS"));
        basisset->print_detail();
    }

    // Set this early because the callback mechanism uses it.
    Process::environment.set_wavefunction(scf_pc);

    if (pre)
        scf_pc->add_preiteration_callback(pre);
    if (post)
        scf_pc->add_postiteration_callback(post);

    energy = scf_pc->compute_energy();

//    // Print a molden file
//    if ( options.get_bool("MOLDEN_WRITE") ) {
//       boost::shared_ptr<MoldenWriter> molden(new MoldenWriter(scf_pc));
//       std::string filename = get_writer_file_prefix() + ".molden";
//       RHF_PC* hf = (RHF_PC*)scf_pc.get();
//       SharedVector occA = hf->occupation_a();
//       SharedVector occB = hf->occupation_b();
//       molden->write(filename, scf_pc->Ca(), scf_pc->Cb(), scf_pc->epsilon_a(), scf_pc->epsilon_b(),occA,occB);
//    }
//
    // Print molecular orbitals
    if ( options.get_bool("PRINT_MOS") ) {
       boost::shared_ptr<MOWriter> mo(new MOWriter(scf_pc,options));
       mo->write();
    }

    // Set some environment variables
    Process::environment.globals["SCF TOTAL ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = energy;

    // Shut down psi.

    tstop();


    return Success;
}

}} // End namespaces

