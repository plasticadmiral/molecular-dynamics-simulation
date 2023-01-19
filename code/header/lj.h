
#ifndef __LJ_H
#define __LJ_H

#include <math.h>
#include <Eigen/Dense>
#include "atoms.h"
#include "neighbors.h"

/**
 * @brief           Implements Lennard-Jones potential with per atom potential.          
 * 
 * @param[in,out]   A       Atoms of the system.
 * @param[in,out]   N       Neighbour lists object to operate with.
 * @param[in]       epsilon Measure of how strongly two particles attract.
 * @param[in]       sigma   Distance at which Intermolecular potential is zero.
 * @return          Lennard-Jones potential energy of the system. 
 */
double lj(Atoms &A, NeighborList &N, const double epsilon = 1.0, 
            const double sigma = 1.0);


#endif  // _LJ_H