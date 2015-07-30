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

#ifndef RHF_PC_H
#define RHF_PC_H

#include <libpsio/psio.hpp>
#include "libscf_solver/rhf.h"
#include "rhf_pc.h"
#include "PointCharges.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class TwoBodySOInt;
class PSIO;
class Chkpt;
class Matrix;
class Vector;

namespace scf {

class RHF_PC : public RHF {
public:
    // PointCharges object
    boost::shared_ptr<PointCharges> rhf_pc_;

protected:
    virtual void form_H();
    virtual double compute_E();

    /// Common initializer
    void common_init();

public:
    RHF_PC(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    RHF_PC(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~RHF_PC();
};

}}

#endif
