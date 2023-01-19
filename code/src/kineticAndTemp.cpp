#include "../header/kineticAndTemp.h"

/**
 * @brief           Calculates kinetic energy of the system.
 * 
 * @param[in,out]   A   Atoms of the system.
 * @return          Kinetic energy of the system.
 */
double getKineticE(Atoms &A) {
    A.kin_energies.setZero();
    for (int i{0}; i < A.velocities.cols(); ++i) {
        for (int j{0}; j < 3; ++j) 
            A.kin_energies(i) +=  pow(A.velocities(j,i), 2);
        A.kin_energies(i) *= A.mass * 0.5;
    }
    return A.kin_energies.sum();
}

/**
 * @brief           Calculates the temperature of the system using velocity of 
 *                  it's atoms.
 * 
 * @param[in,out]       A   Atoms of the system.
 * @return          Temperature of the system. 
 */
double getSystemT(Atoms &A) {    
    return (getKineticE(A) / (A.nb_atoms() * 8.617333262e-5 * 3/2));
}

/**
 * @brief           Calculates the temperature of the system using velocity of 
 *                  it's atoms.
 * @param               KE  Kinetic energy of the system.
 * @param               nb  number of atoms in the system.  
 * @return          Temperature of the system 
 */
double getSystemT(const double KE, const unsigned int nb) {    
    return (KE / (nb * 8.617333262e-5 * 3/2));
}

/**
 * @brief           Calculates the temperature of the system using velocity of 
 *                  it's atoms without boltzmann constant.
 * 
 * @param[in]       A   Atoms of the system.
 * @return          Temperature of the system. 
 */
double getAltSystemT(Atoms &A) {    
    return (getKineticE(A) / (A.nb_atoms() * 1 * 3/2));
}