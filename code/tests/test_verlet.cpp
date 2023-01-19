#include <cstdlib>
#include <gtest/gtest.h>
#include "../header/atoms.h"
#include "../header/verlet.h"

TEST(VerletTest, Test) {
// Sets dt and a random steps in integer range [200, 300)
    double dt = 0.1;
    int steps = rand() % 300 + 200; 

// Creates in on random within integer range [1000, 1200) and [50, 70) 
// number of atoms and masses of atoms respectively. 
    
// Generates Forces, Positions and Velocities in floating point 
// range [-1,1] and creates the atoms struct. 
    int nb_atoms = rand() % 200 + 1000; 
    double M = rand() % 20 + 50;
    Forces_t F = Forces_t::Random(3, nb_atoms); 
    Positions_t P = Positions_t::Random(3,nb_atoms); 
    Velocities_t V = Velocities_t::Random(3, nb_atoms); 
    Atoms atoms(P, V, F, M);
    
// Runs atoms through the velocity-verlet equation.
    for (int i = 0; i < steps; ++i) {
    verletStep1(atoms, dt);
    verletStep2(atoms, dt); 
    }
    
// Compares with positions and velocities with their respective equations 
    double T = dt * steps;
    for (int i = 0; i < nb_atoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(atoms.positions(j, i), 
                        P(j, i)+V(j, i)*T+ 0.5*(1/M)*F(j, i)*T*T, 1e-9);

            EXPECT_NEAR(atoms.velocities(j, i), 
                        V(j, i)+(1/M)*F(j, i)*T, 1e-9);
        }
    }
}
