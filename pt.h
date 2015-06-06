#include "gw.h"
#include "mpi.h"

void pt_swap(gw_ensemble * s0, gw_ensemble * s1);
void pt_mpi_temper(gw_ensemble * self, gw_ensemble * partner);