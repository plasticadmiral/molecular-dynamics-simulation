
#include <gtest/gtest.h>
#include <vector>
#include "../header/atoms.h"
#include "../header/bThermostat.h"
#include "../header/ljDS.h"
#include "../header/verlet.h"
#include "../header/kineticAndTemp.h"
#include <iostream>

// Temperature changes to targetT more rapidly for lower values of tau
TEST(brendsenTest, CoolingTest) {

    double epsilon = 1;
    double sigma = 1;
    std::vector<double>systemT;
    double targetT = 0;
    double taus[] = {1, 3, 6};
    double dt = 0.01;
    int steps = 200;
    
    Positions_t positions(3, 1);
    positions << 0,
                0,
                0;
    Velocities_t velocities(3, 1);
    velocities << 100,
                0,
                0;
    
    int mass = 1;
    
    Atoms atoms(positions, velocities, mass);
    for (int i = 0; i < steps; ++i) {
        verletStep1(atoms, dt);
        ljDS(atoms, epsilon, sigma);
        verletStep2(atoms, dt);
    }
    systemT.push_back(getAltSystemT(atoms));
    for (int tau : taus) {
        Atoms atoms(positions, velocities, mass);
        for (int i = 0; i < steps; ++i) {
            verletStep1(atoms, dt);
            ljDS(atoms, epsilon, sigma);
            verletStep2(atoms, dt);
            altBThermostat(atoms,targetT, dt, tau);
        }
        
        systemT.push_back(getAltSystemT(atoms));
    }
    std::cout << systemT[0] << " " << systemT[1] << " " <<systemT[2] << " " << systemT[3] << std::endl;
    EXPECT_LT(systemT[1], systemT[2]);
    EXPECT_LT(systemT[2], systemT[3]);
    EXPECT_LT(systemT[1], systemT[0]);

}

TEST(brendsenTest, HeatingTest) {
    double epsilon = 1;
    double sigma = 1;
    std::vector<double>systemT;
    double targetT = 10000;
    double taus[] = {1, 3, 6};
    double dt = 0.01;
    int steps = 200;
    Positions_t positions(3, 1);
    positions << 0,
                0,
                0;
    Velocities_t velocities(3, 1);
    velocities << 100,
                0,
                0;
    int mass = 1;
    Atoms atoms(positions, velocities, mass);
    for (int i = 0; i < steps; ++i) {
        verletStep1(atoms, dt);
        ljDS(atoms, epsilon, sigma);
        verletStep2(atoms, dt);
    }
    systemT.push_back(getAltSystemT(atoms));
                
    for (int tau : taus) {
        Atoms atoms(positions, velocities, mass);
        for (int i = 0; i < steps; ++i) {

            verletStep1(atoms, dt);
            ljDS(atoms, epsilon, sigma);
            verletStep2(atoms, dt);
            altBThermostat(atoms, targetT, dt, tau);

        }
        systemT.push_back(getAltSystemT(atoms));
    }

    std::cout << systemT[0] << " " << systemT[1] << " " <<systemT[2] << " " << systemT[3] << std::endl;
    EXPECT_GT(systemT[1], systemT[2]);
    EXPECT_GT(systemT[2], systemT[3]);
    EXPECT_GT(systemT[1], systemT[0]);
}


