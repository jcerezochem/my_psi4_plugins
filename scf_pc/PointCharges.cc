
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

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>
#include <sstream>
#include <string>

#include "PointCharges.h"
#include "psi4-dec.h" //Gives us psi::outfile

#include <physconst.h>

using namespace std;
using namespace boost;

namespace psi { 
extern int str_to_int(const std::string& s);
extern double str_to_double(const std::string& s);

PointCharges::PointCharges(Options &options, boost::shared_ptr<PSIO> psio, int nirrep, boost::shared_ptr<BasisSet> basisset)
{
    outfile->Printf("  **Adding Point Charges to the Hamiltonian**\n");
  
    if(nirrep > 1)
      throw PSIEXCEPTION("You must add\n\n\tsymmetry c1\n\nto the molecule{} block to run use PointCharges.");
  
    basisset_ = basisset; 
  
    boost::shared_ptr<IntegralFactory>
      integrals(new IntegralFactory(basisset, basisset, basisset, basisset));
  
    PetiteList petite(basisset, integrals, true);
    my_aotoso_ = petite.aotoso();
  
    potential_int_ = static_cast<PCMPotentialInt*>(integrals->pcm_potentialint());
  
    std::string pointcharges_file = options.get_str("POINTCHARGES_FILE");
    read_pointcharges(pointcharges_file);
    print_pointcharges();
  
    // Compute interaction with nuclei at initailization
    pointchrg_nuc_int_ = pointcharge_nuclei_int_compute();
  
} // PointCharges()

PointCharges::~PointCharges()
{
    delete [] charges_;
}

double PointCharges::pointcharge_nuclei_int_compute()
{
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    int natom = molecule->natom();

    double ** pcharge_xyz = charge_Zxyz_->pointer();

    double e=0.0;

    for (int atom=0; atom<natom; ++atom) {
        Vector3 xyz_atm = molecule->xyz(atom);
        for (int pointcharge=0; pointcharge<ncharges_; ++pointcharge) {
            Vector3 xyz_pc(pcharge_xyz[pointcharge][1],
                           pcharge_xyz[pointcharge][2],
                           pcharge_xyz[pointcharge][3]); 
            double Zi = molecule->Z(atom);
            double Zj = charges_[pointcharge];
            double distance = xyz_atm.distance(xyz_pc);
            e += Zi * Zj / distance;
        }
    }
    outfile->Printf("\n PointCharges-Nuclei interaction = %17.12f\n\n",e);

    return e;
}

double PointCharges::pointcharge_nuclei_int_get()
{
    return pointchrg_nuc_int_;
}

SharedMatrix PointCharges::compute_V()
{
    // Set the charge field into the PCMPotentialInt object (this could be done only once)
    potential_int_->set_charge_field(charge_Zxyz_);

    // Compute the electrostatic integrals return the contribution to Fock matrix
    SharedMatrix V_pc_cart = SharedMatrix(new Matrix("Point Charges potential cart", basisset_->nao(), basisset_->nao()));
    ContractOverChargesFunctor contract_charges_functor(charges_, V_pc_cart);
    potential_int_->compute(contract_charges_functor);

    // The potential might need to be transformed to the spherical harmonic basis
    SharedMatrix V_pc_pure;
    if(basisset_->has_puream()){
        V_pc_pure = SharedMatrix(new Matrix("Point Charges potential pure", basisset_->nbf(), basisset_->nbf()));
        V_pc_pure->transform(V_pc_cart, my_aotoso_);
    }
    if(basisset_->has_puream()) return V_pc_pure;
    else return V_pc_cart;
}

void PointCharges::read_pointcharges(std::string filename)
{
    // Entire file.
    vector<string> lines;
    // temp
    string text;

    outfile->Printf("    Charges read from: %s\n",filename.c_str());
    // Stream to use (taken from Matrix::load)
    ifstream infile(filename.c_str());
    if (!infile)
        throw PSIEXCEPTION("PointCharges::read_pointcharges: Unable to open file " + filename);

    // reading in entire file
    while (infile.good()) {
        getline(infile, text);
        trim(text);
        if (!text.empty())
            lines.push_back(text);
    }

    ncharges_ = str_to_int(lines[0]);
    //ncharges_=1;
  
    //Allocate arrays
    charge_Zxyz_ = SharedMatrix(new Matrix("Point Charges xyz", ncharges_, 4));
    charges_ = new double[ncharges_];
  
    //double q = 0.5;
    //double x = 0.0;
    //double y = 0.0;
    //double z = 0.800000000000;

    // Go through the file grabbing the data.
    smatch match;
    regex element_line("^\\s*(-*\\d*.\\d*)\\s*(-*\\d*.\\d*)\\s*(-*\\d*.\\d*)\\s*" NUMBER);

    double q, x, y, z;
    for (int elem=0; elem<ncharges_; ++elem) {
        if (regex_match(lines[elem+1], match, element_line)) {
            x = str_to_double(match[1]);
            y = str_to_double(match[2]);
            z = str_to_double(match[3]);
            q = str_to_double(match[4]);
        }
        else
            throw PSIEXCEPTION("PointCharges::read_pointcharges: Unable to match the following line:\n" + lines[elem+1]);

        // Now update the arrays
        charges_[elem] = q;
        charge_Zxyz_->set(elem,0,1.0); // the charge has to be given only once (already in charges_)
        charge_Zxyz_->set(elem,1,x/pc_bohr2angstroms);
        charge_Zxyz_->set(elem,2,y/pc_bohr2angstroms);
        charge_Zxyz_->set(elem,3,z/pc_bohr2angstroms);
    }

    return;

}

void PointCharges:: print_pointcharges()
{           

    double ** Zxyz = charge_Zxyz_->pointer();
              
    outfile -> Printf("\n ==> Point Charges values and coordinates <==\n");
    outfile -> Printf("\n     Charge in |e|; coordinates in Angstrom\n\n");
    outfile->Printf("       Charge              X                  Y                   Z       \n");
    outfile->Printf("    ------------   -----------------  -----------------  -----------------\n");
    for(int i = 0; i < ncharges_; ++i){
        outfile->Printf( " %12.4f    ",charges_[i]);
        for(int j = 1; j < 4; j++)
            outfile->Printf( "  %17.12f", Zxyz[i][j]*pc_bohr2angstroms);
        outfile->Printf("\n");
    }
    outfile->Printf("\n");
  
    return;

}


} // psi namespace

