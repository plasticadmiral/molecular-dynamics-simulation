#ifndef __LJDS_H
#define __LJDS_H

#include <math.h>
#include <Eigen/Dense>
#include "atoms.h"

/**
 * @brief           Implements Lennard-Jones potential without per atom potential.          
 * 
 * @param[in,out]   A       Atoms of the system.
 * @param[in]       epsilon Measure of how strongly two particles attract.
 * @param[in]       sigma   Distance at which Intermolecular potential is zero.
 * @return          Lennard-Jones potential energy of the system. 
 */
double ljDS(Atoms &A, const double epsilon = 1, const double sigma = 1);

#endif  // _LJDS_H