#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

void tdspectra(double DE, psi::Vector &G1_,psi::Vector &G2_,psi::Matrix &GM_,psi::Vector &GK_,
               double* dip0, double* dipm0, double** dipq, double** dipqm,
               psi::Options& options);
