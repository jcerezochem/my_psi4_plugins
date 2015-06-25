#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

/* Computes the rotation matrix to minimize the rmsd  *
 * between structures geo1 and geo2 i                 */
psi::Matrix rot_to_fit_rmsd(psi::Matrix &geo1, psi::Matrix &geo2);

/* Apply a 3D rotation to a (3N x M)                  */
void rotate3d_3Nmat(psi::Matrix& T, psi::Matrix& rot); 

