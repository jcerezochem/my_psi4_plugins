/*
 *@BEGIN LICENSE
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>
#include <physconst.h>

#include <libfunctional/superfunctional.h>
#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <liboptions/liboptions_python.h>
#include <psifiles.h>
#include <libfock/jk.h>
#include <libpsipcm/psipcm.h>


#include <libmints/basisset_parser.h>
#include <libmints/mints.h>
#include <libfock/jk.h>
#include "libtrans/integraltransform.h"
#include "libdpd/dpd.h"

#include "libscf_solver/rhf.h"
#include "rhf_pc.h"
#include "PointCharges.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi { namespace scf {

RHF_PC::RHF_PC(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
    : RHF(options, psio, chkpt)
{
    common_init();
}

RHF_PC::RHF_PC(Options& options, boost::shared_ptr<PSIO> psio)
    : RHF(options, psio)
{
    common_init();
}

RHF_PC::~RHF_PC()
{
}

void RHF_PC::common_init()
{
    rhf_pc_ = static_cast<SharedPointCharges>(new PointCharges(options_, psio_, nirrep_, basisset_));
}

void RHF_PC::form_H()
{
    T_ = SharedMatrix(factory_->create_matrix(PSIF_SO_T));
    V_ = SharedMatrix(factory_->create_matrix(PSIF_SO_V));

    // Assumes these have already been created and stored
    T_->load(psio_, PSIF_OEI);
    V_->load(psio_, PSIF_OEI);

    if (debug_ > 2)
        T_->print("outfile");

    if (debug_ > 2)
        V_->print("outfile");

    if (perturb_h_) {
      if(perturb_ == embpot || perturb_ == sphere || perturb_ == dx) { // embedding potential read from file
        if(nirrep_ > 1)
          throw PSIEXCEPTION("RHF_embed: embedding, dx, and spherical potentials require 'symmetry c1'.");
        int nso = 0;
        for(int h=0; h < nirrep_; h++) nso += nsopi_[h];
        int nao = basisset_->nao();

        // Set up AO->SO transformation matrix (u)
        MintsHelper helper(options_, 0);
        SharedMatrix aotoso = helper.petite_list(true)->aotoso();
        int *col_offset = new int[nirrep_];
        col_offset[0] = 0;
        for(int h=1; h < nirrep_; h++)
          col_offset[h] = col_offset[h-1] + aotoso->coldim(h-1);

        double **u = block_matrix(nao, nso);
        for(int h=0; h < nirrep_; h++)
          for(int j=0; j < aotoso->coldim(h); j++)
            for(int i=0; i < nao; i++)
              u[i][j+col_offset[h]] = aotoso->get(h, i, j);
        delete[] col_offset;

        double *phi_ao, *phi_so, **V_eff;
        phi_ao = init_array(nao);
        phi_so = init_array(nso);
        V_eff = block_matrix(nso, nso);

        if(perturb_ == embpot) {

          FILE* input = fopen("EMBPOT", "r");
          int npoints;
          int statusvalue=fscanf(input, "%d", &npoints);
          outfile->Printf( "  npoints = %d\n", npoints);
          double x, y, z, w, v;
          double max = 0;
          for(int k=0; k < npoints; k++) {
            statusvalue=fscanf(input, "%lf %lf %lf %lf %lf", &x, &y, &z, &w, &v);
            if(fabs(v) > max) max = fabs(v);

            basisset_->compute_phi(phi_ao, x, y, z);
            // Transform phi_ao to SO basis
            C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1, 0.0, &(phi_so[0]), 1);
            for(int i=0; i < nso; i++)
              for(int j=0; j < nso; j++)
                V_eff[i][j] += w * v * phi_so[i] * phi_so[j];
          } // npoints

          outfile->Printf( "  Max. embpot value = %20.10f\n", max);
          fclose(input);

        } // embpot
        else if(perturb_ == dx) {
          dx_read(V_eff, phi_ao, phi_so, nao, nso, u);

        } // dx file
        else if(perturb_ == sphere) {
          radius_ = options_.get_double("RADIUS");
          thickness_ = options_.get_double("THICKNESS");
          r_points_ = options_.get_int("R_POINTS");
          theta_points_ = options_.get_int("THETA_POINTS");
          phi_points_ = options_.get_int("PHI_POINTS");
          outfile->Printf( "  Hard spherical potential radius         = %3.2f bohr\n", radius_);
          outfile->Printf( "  Spherical potential thickness           = %3.2f bohr\n", thickness_);
          outfile->Printf( "  Number of radial integration points     = %d\n", r_points_);
          outfile->Printf( "  Number of colatitude integration points = %d\n", theta_points_);
          outfile->Printf( "  Number of azimuthal integration points  = %d\n", phi_points_);

          double r_step = thickness_/r_points_; // bohr
          double theta_step = 2*pc_pi/theta_points_; // 1 degree in radians
          double phi_step = 2*pc_pi/phi_points_; // 1 degree in radians
          double weight = r_step * theta_step * phi_step;
          for(double r=radius_; r < radius_+thickness_; r += r_step) {
            for(double theta=0.0; theta < pc_pi; theta += theta_step) {  /* colatitude */
              for(double phi=0.0; phi < 2*pc_pi; phi += phi_step) { /* azimuthal */

                double x = r * sin(theta) * cos(phi);
                double y = r * sin(theta) * sin(phi);
                double z = r * cos(theta);

                double jacobian = weight * r * r * sin(theta);

                basisset_->compute_phi(phi_ao, x, y, z);

                C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1,
                        0.0, &(phi_so[0]), 1);

                for(int i=0; i < nso; i++)
                  for(int j=0; j < nso; j++)
                    V_eff[i][j] += jacobian * (-1.0e6) * phi_so[i] * phi_so[j];
              }
            }
          }
        } // sphere


          outfile->Printf( "  Perturbing H by %f V_eff.\n", lambda_);
          if(options_.get_int("PRINT") > 3) mat_print(V_eff, nso, nso, "outfile");


        if(perturb_ == dx) {
          for(int i=0; i < nso; i++)
            for(int j=0; j < nso; j++)
              V_->set(i, j, V_eff[i][j]); // ignore nuclear potential
        }
        else {
          for(int i=0; i < nso; i++)
            for(int j=0; j < nso; j++)
              V_->set(i, j, (V_eff[i][j] + V_->get(i,j)));
        }

        free(phi_ao);
        free(phi_so);
        free_block(V_eff);
      }  // embpot or sphere
    } // end perturb_h_

    // If an external field exists, add it to the one-electron Hamiltonian
    boost::python::object pyExtern = dynamic_cast<PythonDataType*>(options_["EXTERN"].get())->to_python();
    boost::shared_ptr<ExternalPotential> external = boost::python::extract<boost::shared_ptr<ExternalPotential> >(pyExtern);
    if (external) {
        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY") == false && H_->nirrep() != 1)
            throw PSIEXCEPTION("SCF: External Fields are not consistent with symmetry. Set symmetry c1.");

        SharedMatrix Vprime = external->computePotentialMatrix(basisset_);

        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY")) {
            // Attempt to apply symmetry. No error checking is performed.
            SharedMatrix Vprimesym = factory_->create_shared_matrix("External Potential");
            Vprimesym->apply_symmetry(Vprime, AO2SO_);
            Vprime = Vprimesym;
        }

        if (print_) {
            external->set_print(print_);
            external->print();
        }
        if (print_ > 3)
            Vprime->print();
        V_->add(Vprime);


        // Extra nuclear repulsion
        double enuc2 = external->computeNuclearEnergy(molecule_);
        if (print_) {
               outfile->Printf( "  Old nuclear repulsion        = %20.15f\n", nuclearrep_);
               outfile->Printf( "  Additional nuclear repulsion = %20.15f\n", enuc2);
               outfile->Printf( "  Total nuclear repulsion      = %20.15f\n\n", nuclearrep_ + enuc2);
        }
        nuclearrep_ += enuc2;

    }  // end external

    // Save perturbed V_ for future (e.g. correlated) calcs
    V_->save(psio_, PSIF_OEI);

    H_->copy(T_);
    H_->add(V_);
    // Add PointCharges potential
    SharedMatrix V_pc;
    V_pc = rhf_pc_->compute_V();
    H_->add(V_pc);

    if (print_ > 3)
        H_->print("outfile");
}

double RHF_PC::compute_E()
{
    double one_electron_E = 2.0 * D_->vector_dot(H_);
    double two_electron_E = D_->vector_dot(Fa_) - 0.5 * one_electron_E;

    double pointcharge_nuc_E = rhf_pc_->pointcharge_nuclei_int_get();

    energies_["Nuclear"] = nuclearrep_ + pointcharge_nuc_E;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"] = 0.0;
    energies_["-D"] = 0.0;

    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(Fa_);

    double Etotal = nuclearrep_ + pointcharge_nuc_E + D_->vector_dot(HplusF);
    return Etotal;
}


}}//end namespaces

