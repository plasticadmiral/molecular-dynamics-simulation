
#ifndef __KINETICANDTEMP_H
#define __KINETICANDTEMP_H

#include <math.h>
#include "atoms.h"

/**
 * @brief           Calculates kinetic energy in the system.
 * 
 * @param[in,out]   A   Atoms of the system.
 * @return          Kinetic energy of the system.
 */
double getKineticE(Atoms &A);

/**
 * @brief           calculates the temperature of the system using velocity of 
 *                  it's atoms.
 * 
 * @param[in]       A   Atoms of the system.
 * @return          Temperature of the system. 
 */
double getSystemT(Atoms &A);

/**
 * @brief           Calculates the temperature of the system using velocity of 
 *                  it's atoms.
 * @param               KE  Kinetic energy of the system.
 * @param               nb  number of atoms in the system.  
 * @return          Temperature of the system 
 */
double getSystemT(const double KE, const unsigned int nb);   

/**
 * @brief           calculates the temperature of the system using velocity of 
 *                  it's atoms without boltzmann constant.
 * 
 * @param[in]       A   Atoms of the system.
 * @return          Temperature of the system. 
 */
double getAltSystemT(Atoms &A);

#endif  // __KINETICANDTEMP_H