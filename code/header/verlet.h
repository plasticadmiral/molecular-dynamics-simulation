#ifndef __VERLET_H
#define __VERLET_H

#include "atoms.h"

/**
 * @brief           Implements step 1 of velocity verlet algorithm in 
 *                  predictor-corrector scheme.
 * 
 * @param[in,out]   A   Atoms of the system. 
 * @param[in]       TS  Simulation timestep.
 */
void verletStep1( Atoms& A, const double dt);

/**
 * @brief           Implements step 2 of velocity verlet algorithm in 
 *                  predictor-corrector scheme.
 * 
 * @param[in,out]   A   Atoms of the system. 
 * @param[in]       TS  Simulation timestep.
 */
void verletStep2( Atoms& A, const double dt);               

#endif  // __VERLET_H