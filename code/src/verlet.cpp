#include "../header/verlet.h"

/**
 * @brief           Implements step 1 of velocity verlet algorithm in 
 *                  predictor-corrector scheme.
 * 
 * @param[in,out]   A   Atoms of the system. 
 * @param[in]       dt  Simulation timestep.
 */
void verletStep1( Atoms& A, const double dt) {
    A.velocities += 0.5 * A.forces * dt / A.mass;
    A.positions += A.velocities * dt;
}

/**
 * @brief           Implements step 2 of velocity verlet algorithm in 
 *                  predictor-corrector scheme.
 * 
 * @param[in,out]   A   Atoms of the system. 
 * @param[in]       dt  Simulation timestep.
 */
void verletStep2( Atoms& A, const double dt) {
    A.velocities += 0.5 * A.forces * dt / A.mass; 
}