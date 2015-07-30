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

#ifndef PointCharges_H
#define PointCharges_H

#include <vector>
#include <libmints/mints.h>
#include <libmints/potentialint.h>
//#include "potentialint.h"
#include <libpsio/psio.h>
#include <boost/shared_ptr.hpp>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PointCharges {
  public:
    PointCharges() {};
    PointCharges(Options &options, boost::shared_ptr<PSIO> psio, int nirrep, boost::shared_ptr<BasisSet> basisset);
    ~PointCharges();

    // Functions to get potential
    SharedMatrix compute_V();
    double pointcharge_nuclei_int_compute();
    double pointcharge_nuclei_int_get();

  protected:
    /// The number of point charges.
    int ncharges_;
    /// A matrix to hold the charges and {x,y,z} coordinates of the point charged
    SharedMatrix charge_Zxyz_;
    /// An array to hold the charges
    double * charges_;
    /// An array to hold the charges
    double pointchrg_nuc_int_;

    /// Current basis set (for puream and nao/nso info)
    boost::shared_ptr<BasisSet> basisset_;

    /// The AO->SO transformation matrix, which is used for transforming
    /// matrices between pure and Cartesian representations.
    SharedMatrix my_aotoso_;

    /// Factory for the electrostatic integrals
    PCMPotentialInt* potential_int_;

    /// Read/Print point charges 
    void read_pointcharges(std::string filename);
    void print_pointcharges();

};

typedef boost::shared_ptr<psi::PointCharges> SharedPointCharges;

} // psi
#endif
