#ifndef __BTHERMOSTAT_H
#define __BTHERMOSTAT_H

#include <math.h>
#include "atoms.h"
#include "kineticAndTemp.h"


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
 *                          thermal bath.
 */
void bThermostat(Atoms &A, const double targetT, const double dt, 
                            const double tau);

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
 *                          thermal bath.
 */
void altBThermostat(Atoms &A, const double targetT, const double dt, 
                                const double tau);
#endif  // _BTHERMOSTAT_H