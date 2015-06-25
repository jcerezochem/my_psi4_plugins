/*
 *@BEGIN LICENSE
 *
 * vibronic by Psi4 Developer, a plugin to:
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
// My includes
#include "constants.h"
#include "rotation.h"
#include "tdspectra.h"
//#include "fcc_utils.h"

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace vibronic {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "VIBRONIC"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        // ==> GENERAL OPTIONS <== //
        /*- Temperature in K -*/
        options.add_double("TEMPERATURE", 300.0);
        // ==> APPROACH OPTIONS <==
        /* Model PES */
        options.add_str("MODEL_PES", "AS");
        /* Dipole approximaiton */
        options.add_str("DIPOLE_APPROX", "FC");
        // ==> TRANSITION OPTIONS <==
        /* Transition Type */
        options.add_str("TRANSITION", "ABS");
        /* Chiral spectroscopy */
        options.add_str("CHIRAL", "NO");
        // ==> CORRELATION FUNCITON <== //
        /* Number of points for the correlation function -*/
        options.add_int("CORR_NPOINTS", 1024);
        /* Final time of the correlation functions (fs) -*/
        options.add_int("CORR_TIME", 1000.0);
        // ==> INHOMOGENEOUS BROADENING <== //
        /* HWHM of the convoluted inhomogeneous broadening (eV) */
        options.add_double("BROAD_HWHM", 0.01);
        /* Inhomogeneous broadening type (GAUSSIAN/LORENTZIAN) */
        options.add_str("BROAD_TYPE", "GAUSSIAN");
    }

    return true;
}

extern "C" 
PsiReturnType vibronic(Options& options)
{
    int print = options.get_int("PRINT");

    /* Your code goes here */
    double temp               = options.get_double("TEMPERATURE");
    std::string opdip         = options.get_str("DIPOLE_APPROX");
    std::string optrans       = options.get_str("TRANSITION");
    std::string optrans1      = options.get_str("CHIRAL");
    std::string opmod         = options.get_str("MODEL_PES");
    int corr_npoints          = options.get_int("CORR_NPOINTS");
    double corr_time          = options.get_double("CORR_TIME");
    double inh_brd_hwhm       = options.get_double("BROAD_HWHM");
    std::string inh_brd_type  = options.get_str("BROAD_TYPE");

       
    // Get difference energy (from "CURRENT ENERGY", expected in Hartree) 
    double de = Process::environment.globals["CURRENT ENERGY"];

    outfile->Printf("\n\n"  );
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("                                VIBRONIC                          \n");
    outfile->Printf("              Compute vibrationally resolved electronic spectra   \n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("\n");
    outfile->Printf("     Options \n"                                  );
    outfile->Printf("     -------------------------------------\n"     );
    outfile->Printf("     Transition Type  : %s\n",optrans.c_str()     );
    outfile->Printf("     Chiral option    : %s\n",optrans1.c_str()    );
    outfile->Printf("     Temperature (K)  : %f\n",temp                );
    outfile->Printf("     Model PES        : %s\n",opmod.c_str()       );
    outfile->Printf("     Dipole Aprox.    : %s\n",opdip.c_str()       );
    outfile->Printf("     Evaluation of correlation function\n"        );
    outfile->Printf("        Points        : %i\n",corr_npoints        );
    outfile->Printf("        Time (fs)     : %f\n",corr_time           );
    outfile->Printf("     Ihnomogeneous broadening\n"                  );
    outfile->Printf("        Type          : %s\n",inh_brd_type.c_str());
    outfile->Printf("        HWHM (eV)     : %f\n",inh_brd_hwhm        );
    outfile->Printf("     Adiabatic En.(eV): %f\n",de*autoev           );
    outfile->Printf("     -------------------------------------\n"     );
    outfile->Printf("\n\n"  );
   

    //
    // Get state geometries and vibrations
    // state1 is environment.molecule 
    // state2 is environment.molecule2
    //
    // State1
    shared_ptr<Molecule> state1  = Process::environment.molecule();
    int Nat = state1->natom();
    shared_ptr<Vector> frequencies1 = Process::environment.wavefunction()->frequencies();
    frequencies1->set_name("Frequencies (State 1)");
    int Nvib = frequencies1->dim();
    shared_ptr<Vector> normalmodes1 = Process::environment.wavefunction()->normalmodes();
    // Transform nm vector into matrix
    Matrix T1("Normal Modes (State 1)",3*Nat,Nvib);
    int k=0;
    for (int i=0; i<Nvib; i++) {
        for (int j=0; j<3*Nat; j++) {
            T1(j,i) = normalmodes1->get(k);
            k++;
        }
    }

    // State2
    shared_ptr<Molecule> state2  = Process::environment.molecule2();
    // Only AS for the moment...
    shared_ptr<Vector> frequencies2;
    Matrix T2;
    if (opmod == "AS") {
        frequencies2 = Process::environment.wavefunction()->frequencies();
        T2.copy(T1);
    } 
    else {
        throw PSIEXCEPTION("Only AS model available");
    }
    
    // Print data
    outfile->Printf("  ==> Geometries and frequencies <==\n\n");
    outfile->Printf("    -------------\n");
    outfile->Printf("     State 1 \n"      );
    outfile->Printf("    -------------\n");
    // Geom
    state1->print();
    // Freq
    frequencies1->print();
    // Normal Modes
    T1.print();
    //
    outfile->Printf("    -------------\n");
    outfile->Printf("     State 2 \n"      );
    outfile->Printf("    -------------\n");
    // Geom
    state2->print();
    // Freq
    frequencies2->print();
    // Normal Modes
    T2.print();
    //

    // Preapare for vibronic computation:
    // Set to com
    state1->move_to_com();
    state2->move_to_com();
    // Mass-weight coordinates
    Matrix geo1(Nat,3);
    Matrix geo2(Nat,3);
    for (int i=0; i<Nat; i++) {
        geo1(i,0) = state1->x(i)*sqrt(state1->mass(i)/autoamu);
        geo1(i,1) = state1->y(i)*sqrt(state1->mass(i)/autoamu);
        geo1(i,2) = state1->z(i)*sqrt(state1->mass(i)/autoamu);
        geo2(i,0) = state2->x(i)*sqrt(state2->mass(i)/autoamu);
        geo2(i,1) = state2->y(i)*sqrt(state2->mass(i)/autoamu);
        geo2(i,2) = state2->z(i)*sqrt(state2->mass(i)/autoamu);
    }
    // Convert Freq to AU
    frequencies1->scale(1./autown);
    //frequencies2->scale(1./autown);

    // Rotate state1
    // Compute rotation matrix with Quaternion formalism
    outfile->Printf("  ==> Rotate (RMSD fit) <==\n\n");
    Matrix rot = rot_to_fit_rmsd(geo1, geo2);
    outfile->Printf("     Rotation Matrix\n");
    outfile->Printf("     ----------------\n");
    rot.print();
    // Rotate T1
    rotate3d_3Nmat(T1, rot);
  
    // Compute Duschinsky matrix and disp vector
    Matrix GM("Duschinsky matrix",Nvib,Nvib);
    Vector GK("Normal mode displacements",Nvib);
    GM.gemm(true,false,1.,&T1,&T2,1.);
    Vector vec(3*Nat);
    for (int i=0; i<Nat; i++) {
        int k=3*i;
        vec(k  ) = geo2(i,0)-geo1(i,0);
        vec(k+1) = geo2(i,1)-geo1(i,1);
        vec(k+2) = geo2(i,2)-geo1(i,2);
    }
    for (int i=0; i<Nvib; i++) {
        for (int j=0; j<3*Nat; j++) {
            GK(i) += vec(j)*T1(j,i);
        }
    }
    
    //GK.gemv(true,1.,&T1,&vec,1.);
    GM.print();
    GK.print();


    // Temporary assigment
    de = 2.0/autoev;
    double dip0[3] = {1.0};
    double dipm0[3] = {1.0};
    double** dipq  = new double*[Nvib];
    double** dipmq = new double*[Nvib];
    for (int i=0; i<Nvib; i++) {
        dipq[i]  = new double[3];
        dipmq[i] = new double[3];
    }

    tdspectra(de,*frequencies1,*frequencies2,GM,GK,
              dip0,dipm0,dipq,dipmq,
              options);
    
    return Success;
}

}} // End namespaces

