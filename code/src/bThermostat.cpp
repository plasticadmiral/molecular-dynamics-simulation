#include "../header/bThermostat.h"

/**
 * @brief           Implements a berendsen thermostat.
 * 
 * @param[in,out]   A       Atoms of the system.
 * 
 * @param           targetT Temperature of thermal bath. 
 * 
 * @param           dt      Simulation timestep.
 * 
 * @param           tau     Relaxation time determining the strenght of
 *                          coupling between the atomic system and thermal 
 *                          bath.
 */
void bThermostat(Atoms &A, const double targetT, const double dt, 
                            const double tau) {
    A.velocities *= sqrt(1 + ((targetT/getSystemT(A)) - 1) * (dt/tau));
}

/**
 * @brief           Implements a berendsen thermostat without boltzmann constant.
 * 
 * @param[in,out]   A       Atoms of the system.
 * 
 * @param           targetT Temperature of thermal bath. 
 * 
 * @param           dt      Simulation timestep.
 * 
 * @param           tau     Relaxation time determining the strenght of
 *                          coupling between the atomic system and thermal 
 *                          bath.
 */
void altBThermostat(Atoms &A, const double targetT, const double dt, 
                            const double tau) {
    A.velocities *= sqrt(1 + ((targetT/getAltSystemT(A)) - 1) * (dt/tau));
}