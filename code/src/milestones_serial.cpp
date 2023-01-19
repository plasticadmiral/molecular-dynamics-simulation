#include "../header/milestones_serial.h"

using namespace std;

//4a lennard jones interaction data for plot make.
//5a cooling down and heating up data for gif make.
//7a gupta interaction data for plot make.
//7b make a more exaggerated molten transitioning data for gif. if 7 not sufficient.



/**
 * @brief       Simulates the lj54.xyz file with Lennard Jones potential for 
 *              different dt. 
 *              
 *              "total energy as a function of loops" for different dt 
 *              is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m4.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m4(const vars4& V) {
    cout << "Starting m4!\n\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " for " << V.dtList->size() << " different timesteps.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    unsigned long long int seconds{0};

    for(double dt : *V.dtList) {
        name = V.expName + "_" + to_str(dt);
        ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");
        
        cout << "\nRunning experiment with " << to_str(dt) << " timestep.\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

        auto startTime = chrono::high_resolution_clock::now();
        for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
            if(loop % (V.loops / V.progress_cout_divisions) == 0) 
                cout << "Iteration: " << loop << "/" << V.loops << endl;

            verletStep1(*V.atoms, dt);
            pot_energy = ljDS(*V.atoms, V.epsilon, V.sigma);
            verletStep2(*V.atoms, dt);

            if(loop % V.csv_loops == 0) {
                    kin_energy = getKineticE(*V.atoms);
                    total_energy = pot_energy + kin_energy;
                    temperature = getSystemT(kin_energy, V.atoms->nb_atoms());

                    csv_file << loop << "\t" << total_energy << "\t"
                        << pot_energy << "\t" << kin_energy << "\t"
                        << temperature << "\n";

                    write_xyz(xyz_file, *V.atoms);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();

        elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                    (endTime - startTime));
        
        cout << "Execution time: " << elapsedTime.back().count() 
                                    << " seconds\n";

        csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" << V.atoms->mass    << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << dt               << "\t0\t0\t0\n";
        csv_file << "eps:"      << "\t" << V.epsilon        << "\t0\t0\t0\n";
        csv_file << "sig:"      << "\t" << V.sigma          << "\t0\t0\t0\n";
        
        xyz_file.close();
        csv_file.close();

    }
    seconds = chrono_time_sum(elapsedTime);
    
    cout << "Completed m4!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds <<"\n\n";
    
    return seconds;
}

/**
 * @brief       Moves two atoms away from each other at constant velocity 
 *              while calculating the lennard jones potential between them.
 *              
 *              "Potential energy vs distance" is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m4a.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m4a(const vars4a& V) {
    cout << "Starting m4a!\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " to plot the lennard jones potential.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    if(V.epsilonList->size() != V.sigmaList->size() && 
                V.atoms->size() != V.sigmaList->size()) return 0;

    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    double sigma{0}, epsilon{0};
    unsigned long long int seconds{0};

    for(unsigned int i{0}; i < V.sigmaList->size(); ++i) {
        sigma = (*V.sigmaList)[i];
        epsilon = (*V.epsilonList)[i];
        Atoms atoms{(*V.atoms)[i]};

        name = to_str(V.expName) + "_E_" + to_str(epsilon) + "_S_" 
                                 + to_str(sigma) + "_dt_" + to_str(V.dt);
        ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv"); 

        cout << "Performing simulation with experiment name: " << name << "\n";

        
        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

        auto startTime = chrono::high_resolution_clock::now();
        for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
            if(loop % (V.loops / V.progress_cout_divisions) == 0) 
                cout << "Iteration: " << loop << "/" << V.loops << endl;

            verletStep1(atoms, V.dt);
            pot_energy = ljDS(atoms, epsilon, sigma);
            atoms.forces.setZero();
            verletStep2(atoms, V.dt);

            if(loop % 1 == 0) {
                    kin_energy = getKineticE(atoms);
                    total_energy = pot_energy + kin_energy;
                    temperature = getSystemT(kin_energy, atoms.nb_atoms());

                    csv_file << loop << "\t" << total_energy << "\t"
                        << pot_energy << "\t" << kin_energy << "\t"
                        << temperature << "\n";

                    write_xyz(xyz_file, atoms);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();

        elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                    (endTime - startTime));
        
        cout << "Execution time: " << elapsedTime.back().count() 
                                    << " seconds \n";

        csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" << atoms.mass       << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
        csv_file << "eps:"      << "\t" << epsilon          << "\t0\t0\t0\n";
        csv_file << "sig:"      << "\t" << sigma            << "\t0\t0\t0\n";
        
        xyz_file.close();
        csv_file.close();
    }
    seconds = chrono_time_sum(elapsedTime);

    cout << "Completed m4a!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds << "\n\n";
    
    return seconds;
}

/**
 * @brief       Simulates atoms arranged in cubic lattices of different sizes 
 *              with the Lennard Jones potential to measure computation time. 
 *              
 *              "simulation time vs number of atoms in system" 
 *              is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m5.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m5(const vars5& V) {
    cout << "Starting m5!\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " for " << V.atomsList->size() 
            << " different cubic lattice atoms sets.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";
    
    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    double sigma{0}, epsilon{0};
    unsigned long long int seconds{0};

    for(Atoms atoms : *V.atomsList) {
        name = V.expName + "_" + to_string(atoms.nb_atoms());
        ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");

        cout << "\nRunning experiment with " << to_string(atoms.nb_atoms()) 
                                                << " atoms\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

        auto startTime = chrono::high_resolution_clock::now();
        for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
            if(loop % (V.loops / V.progress_cout_divisions) == 0) 
                cout << "Iteration: " << loop << "/" << V.loops << endl;

            verletStep1(atoms, V.dt);
            pot_energy = ljDS(atoms, V.epsilon, V.sigma);
            verletStep2(atoms, V.dt);

            if(loop % V.csv_loops == 0) {
                kin_energy = getKineticE(atoms);
                total_energy = pot_energy + kin_energy;
                temperature = getSystemT(kin_energy, atoms.nb_atoms());

                csv_file << loop << "\t" << total_energy << "\t"
                    << pot_energy << "\t" << kin_energy << "\t"
                    << temperature << "\n";

                write_xyz(xyz_file, atoms);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();

        elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                    (endTime - startTime));
        
        cout << "Execution time: " << elapsedTime.back().count() 
                                    << " seconds\n";
        
        csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" << atoms.mass       << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
        csv_file << "eps:"      << "\t" << V.epsilon        << "\t0\t0\t0\n";
        csv_file << "sig:"      << "\t" << V.sigma          << "\t0\t0\t0\n";

        xyz_file.close();
        csv_file.close();
    }

    seconds = chrono_time_sum(elapsedTime);
    
    cout << "Completed m5!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds << "\n\n";
    
    return seconds;
}

/**
 * @brief       Simulates a single cubic lattice file at three different 
 *              temperature ranges namely, a natural state, vapourizing state  
 *              and cooled down state temperatures. 
 *              
 *              gifs are made using this data.              
 * 
 * @param V     A struct holding arguments specifically for m5.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m5a(const vars5a& V) {
    cout << "Starting  m5!\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " on cubic lattice with " << V.atoms->nb_atoms() << " atoms.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    vector<int> targetTs{0, V.coolT, V.heatT};
    unsigned long long int seconds{0};

    for(double targetT : targetTs) {
        name = V.expName + "_" + to_str(targetT,0) + "_" 
                         + to_string(V.atoms->nb_atoms());

        ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");

        cout << "\nRunning experiment with " << to_string(V.atoms->nb_atoms()) 
                                                << " atoms\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

        auto startTime = chrono::high_resolution_clock::now();
        for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
            if(loop % (V.loops / V.progress_cout_divisions) == 0)
                cout << "Iteration: " << loop << "/" << V.loops << endl;

            verletStep1(*V.atoms, V.dt);
            pot_energy = ljDS(*V.atoms, V.epsilon, V.sigma);
            verletStep2(*V.atoms, V.dt);
            if(targetT) 
                bThermostat(*V.atoms, targetT, V.dt, V.tau);

            if(loop % V.csv_loops == 0) {
                    kin_energy = getKineticE(*V.atoms);
                    total_energy = pot_energy + kin_energy;
                    temperature = getSystemT(kin_energy, V.atoms->nb_atoms());

                    csv_file << loop << "\t" << total_energy << "\t"
                        << pot_energy << "\t" << kin_energy << "\t"
                        << temperature << "\n";

                    write_xyz(xyz_file, *V.atoms);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();

        elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                    (endTime - startTime));

        cout << "Execution time: " << elapsedTime.back().count() 
                                    << " seconds \n";

        csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" << V.atoms->mass    << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
        csv_file << "eps:"      << "\t" << V.epsilon        << "\t0\t0\t0\n";
        csv_file << "sig:"      << "\t" << V.sigma          << "\t0\t0\t0\n";
        csv_file << "heatT:"    << "\t" << V.heatT          << "\t0\t0\t0\n";
        csv_file << "coolT:"    << "\t" << V.coolT          << "\t0\t0\t0\n";

        xyz_file.close();
        csv_file.close();
    }

    seconds = chrono_time_sum(elapsedTime);
    cout << "Completed m5!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds << "\n\n";
    
    return seconds;
}

/**
 * @brief       Simulates atoms arranged in cubic lattices of different sizes 
 *              with the Lennard Jones potential using Neighbour Lists.
 * 
 *              "simulation time vs number of atoms in system" 
 *              is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m6.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m6(const vars6& V) {
    cout << "Starting m6!\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " for " << V.atomsList->size() 
                << " different cubic lattice atoms sets.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";
    
    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    unsigned long long int seconds{0};

    for(Atoms atoms : *V.atomsList) {
        name = V.expName + "_" + to_string(atoms.nb_atoms());
        ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");

        cout << "Running experiment with " << to_string(atoms.nb_atoms())  
                                    << " atoms\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";
        
        NeighborList neighbors(V.cutoff_range);

        auto startTime = chrono::high_resolution_clock::now();
        for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
            if(loop % (V.loops / V.progress_cout_divisions) == 0) 
                cout << "Iteration: " << loop << "/" << V.loops << endl;

            verletStep1(atoms, V.dt);
            neighbors.update(atoms);
            double pot_energy{lj(atoms, neighbors, V.epsilon, V.sigma)};
            verletStep2(atoms, V.dt);

            if(loop % V.csv_loops == 0) {
                    kin_energy = getKineticE(atoms);
                    total_energy = pot_energy + kin_energy;
                    temperature = getSystemT(kin_energy, atoms.nb_atoms());

                    csv_file << loop << "\t" << total_energy << "\t"
                        << pot_energy << "\t" << kin_energy << "\t"
                        << temperature << "\n";
                    
                    write_xyz(xyz_file, atoms);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();

        elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                    (endTime - startTime));
        
        cout << "Execution time: " << elapsedTime.back().count() 
                                    << " seconds \n";

        csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" <<  atoms.mass      << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
        csv_file << "eps:"      << "\t" << V.epsilon        << "\t0\t0\t0\n";
        csv_file << "sig:"      << "\t" << V.sigma          << "\t0\t0\t0\n";
        csv_file << "cutoff:"   << "\t" << V.cutoff_range   << "\t0\t0\t0\n";
        
        xyz_file.close();
        csv_file.close();
    }


    seconds = chrono_time_sum(elapsedTime);

    cout << "Completed m6!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds << "\n\n";

    return seconds;
}

/**
 * @brief       Simulates the transition from the solid to molten state using 
 *              velocity scaling as means to introduce heat into a gold Mackay 
 *              Icosahedron structure of different cluster sizes. The potential 
 *              energy between atoms is produced using gupta potential.
 * 
 *              "total energy vs temperature",
 *              "melting point vs cluster size", 
 *              "heat capacity vs cluster size" and 
 *              "latent heat vs cluster size" 
 *              are plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m7.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m7(const vars7& V) {
    cout << "Starting m7!\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " for with " << V.atomsList->size() 
            << " different cluster atom sets.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    unsigned long long int seconds{0};
    int i{-1};
    for(Atoms atoms : *V.atomsList) {
        i++;
        NeighborList neighbors(V.cutoff_range);
        name = V.expName + "_" + to_string(atoms.nb_atoms());
        // ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");

        cout << "Running experiment with " << to_string(atoms.nb_atoms()) 
                                            << " atoms cluster.\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";


        auto startTime = chrono::high_resolution_clock::now();
        for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
            if(loop % (V.loops / V.progress_cout_divisions) == 0) 
                cout << "Iteration: " << loop << "/" << V.loops << endl;

            verletStep1(atoms, V.dt);
            neighbors.update(atoms);
            pot_energy = gupta(atoms, neighbors);
            verletStep2(atoms, V.dt);

            if(loop != 0 && loop % V.rel_loops == 0) 
                atoms.velocities *= (((*V.deltaQ)[i]/getKineticE(atoms)) + 1);
            

            if(loop % V.csv_loops == 0) {
                    kin_energy = getKineticE(atoms);
                    total_energy = pot_energy + kin_energy;
                    temperature = getSystemT(kin_energy, atoms.nb_atoms());

                    csv_file << loop << "\t" << total_energy << "\t"
                        << pot_energy << "\t" << kin_energy << "\t" 
                        << temperature <<"\n";

                    // write_xyz(xyz_file, atoms);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();

        elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                    (endTime - startTime));
        
        cout << "Execution time: " << elapsedTime.back().count() 
                                                            << " seconds\n";

        csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" <<  atoms.mass      << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
        csv_file << "cutoff:"   << "\t" << V.cutoff_range   << "\t0\t0\t0\n";
        csv_file << "rsteps:"   << "\t" << V.rel_loops      << "\t0\t0\t0\n";
        csv_file << "deltaq:"   << "\t" << (*V.deltaQ)[i]   << "\t0\t0\t0\n";

        // xyz_file.close();
        csv_file.close();
    }
    // if possible, identify the melting pt, heat cap and latent heat here.

    seconds = chrono_time_sum(elapsedTime);

    cout << "Completed m7!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds <<"\n\n";

    return seconds;
}

/**
 * @brief       Moves two atoms away from each other at constant velocity 
 *              while calculating the gupta potential between them.
 * 
 *              "Potential energy vs distance" is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m7a.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m7a(const vars7a& V) {
    cout << "Starting m7a!\n";
    cout << "Performing simulation with experiment name " << V.expName 
            << " to plot the gupta potential.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    unsigned long long int seconds{0};

    NeighborList neighbors(V.cutoff_range);
    
    name = V.expName + "_" + to_string(V.atoms->nb_atoms());
    ofstream xyz_file(V.outLoc + name + ".xyz"); 
    ofstream csv_file(V.outLoc + name + ".csv");

    cout << "Running experiment with " << to_string(V.atoms->nb_atoms()) 
                                        << " atoms cluster.\n";

    csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

    auto startTime = chrono::high_resolution_clock::now();
    for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
        if(loop % (V.loops / V.progress_cout_divisions) == 0) 
            cout << "Iteration: " << loop << "/" << V.loops << endl;

        verletStep1(*V.atoms, V.dt);
        neighbors.update(*V.atoms);
        pot_energy = gupta(*V.atoms, neighbors);
        V.atoms->forces.setZero();
        verletStep2(*V.atoms, V.dt);

        if(loop % V.csv_loops == 0) {
                kin_energy = getKineticE(*V.atoms);
                total_energy = pot_energy + kin_energy;
                temperature = getSystemT(kin_energy, V.atoms->nb_atoms());
             
                csv_file << loop << "\t" << total_energy << "\t"
                    << pot_energy << "\t" << kin_energy << "\t" 
                    << temperature <<"\n";
             
                write_xyz(xyz_file, *V.atoms);
        }
    }
    auto endTime = chrono::high_resolution_clock::now();

    elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                (endTime - startTime));
    
    cout << "Execution time: " << elapsedTime.back().count() 
                                                        << " seconds\n";

    csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                        << "\ts\t0\t0\n";
    csv_file << "mass:"     << "\t" << V.atoms->mass    << "\t0\t0\t0\n";
    csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
    csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
    csv_file << "cutoff:"   << "\t" << V.cutoff_range   << "\t0\t0\t0\n";

    xyz_file.close();
    csv_file.close();

    seconds = chrono_time_sum(elapsedTime);

    cout << "Completed m7a!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds <<".\n\n";

    return seconds;
}

/**
 * @brief       Performs the simulation from m7 on a single cluster but with 
 *              extra heat introduction and longer relaxation periods to 
 *              make the process more visible.
 *
 *              A gif is made using this data. 
 * 
 * @param V     A struct holding arguments specifically for m7.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m7b(const vars7b& V) {
    cout << "Starting m7b!\n";
    cout << "Performing simulation with experiment name" << V.expName << 
        " on cubic lattice with " << V.atoms->nb_atoms() << " atoms.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsedTime;
    string name;
    double total_energy{0}, pot_energy{0}, kin_energy{0}, temperature{0};
    unsigned long long int seconds{0};

    NeighborList neighbors(V.cutoff_range);

    name = V.expName + "_" + to_string(V.atoms->nb_atoms());
    ofstream xyz_file(V.outLoc + name + ".xyz"); 
    ofstream csv_file(V.outLoc + name + ".csv");

    cout << "Running experiment with " << to_string(V.atoms->nb_atoms()) 
                                        << " atoms\n";

    csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

    auto startTime = chrono::high_resolution_clock::now();
    for (unsigned long long int loop{0}; loop < V.loops; ++loop) {
        if(loop % (V.loops / V.progress_cout_divisions) == 0) 
            cout << "Iteration: " << loop << "/" << V.loops << endl;

        verletStep1(*V.atoms, V.dt);
        neighbors.update(*V.atoms);
        pot_energy = gupta(*V.atoms, neighbors);
        verletStep2(*V.atoms, V.dt);

        // if(loop != 0 && loop % V.rel_loops == 0) 
            // V.atoms->velocities *= ((V.deltaQ/getKineticE(*V.atoms)) + 1);
        

        if(loop % V.csv_loops == 0) {
            kin_energy = getKineticE(*V.atoms);
            total_energy = pot_energy + kin_energy;
            temperature = getSystemT(kin_energy,V.atoms->nb_atoms());

            csv_file << loop << "\t" << total_energy << "\t"
                << pot_energy << "\t" << kin_energy << "\t" 
                << temperature <<"\n";

            write_xyz(xyz_file, *V.atoms);
        }
    }
    auto endTime = chrono::high_resolution_clock::now();

    elapsedTime.push_back(chrono::duration_cast<chrono::seconds>
                                                (endTime - startTime));
    
    cout << "Execution time: " << elapsedTime.back().count() 
                                                        << " seconds\n";

    csv_file << "etime:"    << "\t" << elapsedTime.back().count() 
                                                        << "\ts\t0\t0\n";
    csv_file << "mass:"     << "\t" << V.atoms->mass    << "\t0\t0\t0\n";
    csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
    csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
    csv_file << "cutoff:"   << "\t" << V.cutoff_range   << "\t0\t0\t0\n";
    // csv_file << "rsteps:"   << "\t" << V.rel_loops      << "\t0\t0\t0\n";
    // csv_file << "deltaq:"   << "\t" << V.deltaQ         << "\t0\t0\t0\n";

    xyz_file.close();
    csv_file.close();

    seconds = chrono_time_sum(elapsedTime);

    cout << "Completed m7b!\n";
    cout << "Total execution time in seconds for milestone: " 
            << seconds <<"\n\n";

    return seconds;
}