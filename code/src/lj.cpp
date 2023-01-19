#include "../header/lj.h" 

/**
 * @brief           Implements Lennard-Jones potential with neighbour lists          
 * 
 * @param[in,out]   A       Atoms of the system.
 * @param[in,out]   N       Neighbour lists object to operate with.
 * @param[in]       epsilon Measure of how strongly two particles attract.
 * @param[in]       sigma   Distance at which Intermolecular potential is zero.
 * @return          Lennard-Jones potential energy of the system.
 */
double lj(Atoms &A, NeighborList &N, const double epsilon, const double sigma) {
    A.pot_energies.setZero();
    A.forces.setZero();
    
    for (auto[i, j]: N) {           
        Eigen::Array3d dist_vectr{A.positions.col(i) - A.positions.col(j)};
        double dist{sqrt((dist_vectr * dist_vectr).sum())};

        double derivative{4 * epsilon * (12 * pow(sigma, 12) / pow(dist, 13) - 
                        (6 * pow(sigma, 6) / pow(dist, 7)))};
        
        for(int k{0}; k < 3; ++k)    
            A.forces(k, i) +=  derivative * dist_vectr(k) / dist;

        A.pot_energies(i) += 2 * epsilon * (pow((sigma / dist), 12) - 
                        pow((sigma / dist), 6));
    }
    return A.pot_energies.sum();    
}