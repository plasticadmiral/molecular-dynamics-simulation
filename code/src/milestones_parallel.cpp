#include "../header/milestones_parallel.h"

using namespace std;

/**
 * @brief       parallelizing m7
 * 
 * @param V     A struct holding arguments specifically for m8.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m8(vars8& V) {
    cout << "Starting m8!\n";
    cout << "Performing simulation with experiment name " << V.expName 
                << " for " << V.atomsList->size() 
                << " different cluster atom sets.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsed;
    string name;
    double total_energy{0}, temperature{0};
    unsigned long long int seconds{0};

    for(Atoms atoms : *V.atomsList) {
        NeighborList neighbors(V.cutoff_range);
        Domain domain(MPI_COMM_WORLD, V.domain_size, V.divisions, V.isEnabled);
        
        name = V.expName + "_parallel_" + to_string(atoms.nb_atoms());
        ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");

        cout << "Running experiment with " << to_string(atoms.nb_atoms()) 
                                            << " atoms cluster.\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\n";

        auto begin = chrono::high_resolution_clock::now();
        domain.enable(atoms);
        for (int loop{0}; loop < V.loops; ++loop) {
            double lPotE{0}, lKinE{0}, gPotE{0}, gKinE{0};
            
            if(domain.rank() == 0 && loop % 
                (V.loops / V.progress_cout_divisions) == 0) {
                cout << "Iteration: " << loop << "/" << V.loops << endl;
            }

            verletStep1(atoms, V.dt);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, 2 * V.cutoff_range);
            neighbors.update(atoms);
            gupta(atoms, neighbors);
            verletStep2(atoms, V.dt);

            bool csv_condition{loop % V.csv_loops == 0};
            if(csv_condition) {
                getKineticE(atoms);
                for (int i{0}; i < domain.nb_local(); ++i) {
                    lPotE += atoms.pot_energies(i);
                    lKinE += atoms.kin_energies(i);
                }

                gPotE = MPI::allreduce(lPotE, MPI_SUM, MPI_COMM_WORLD);
                gKinE = MPI::allreduce(lKinE, MPI_SUM, MPI_COMM_WORLD);
            }

            if(domain.rank() == 0 && csv_condition) {
                    total_energy = gPotE + gKinE;
                    temperature = getSystemT(gKinE, atoms.nb_atoms());

                    csv_file << loop << "\t" << total_energy << "\t"
                        << gPotE << "\t" << gKinE << "\t" 
                        << temperature <<"\n";

                    write_xyz(xyz_file, atoms);
            }
        }
        domain.disable(atoms);
        auto end = chrono::high_resolution_clock::now();

        elapsed.push_back(chrono::duration_cast<chrono::seconds>(end - begin));
        
        cout << "Execution time: " << elapsed.back().count() << " seconds\n";

        csv_file << "etime:"    << "\t" << elapsed.back().count() 
                                                            << "\ts\t0\t0\n";
        csv_file << "mass:"     << "\t" << atoms.mass       << "\t0\t0\t0\n";
        csv_file << "loops:"    << "\t" << V.loops          << "\t0\t0\t0\n";
        csv_file << "dt:"       << "\t" << V.dt             << "\t0\t0\t0\n";
        csv_file << "cutoff:"   << "\t" << V.cutoff_range   << "\t0\t0\t0\n";

        xyz_file.close();
        csv_file.close();
    }
    seconds = chrono_time_sum(elapsed);
    cout << "Completed m8!\n";
    cout << "Total execution time in sec for milestone: " << seconds <<"\n\n";
    return seconds;
}

/**
 * @brief       stretching gold wire
 * 
 * @param V     A struct holding arguements specifically for m9.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m9(vars9& V) {

    cout << "Starting m9!\n";
    cout << "Performing simulation with experiment name " << V.expName 
                << " for " << V.atomsList->size() 
                << " different cluster atom sets.\n";
    cout << "Storing the data at: " << V.outLoc << "\n";

    vector<chrono::seconds> elapsed;
    string name;
    Eigen::Array3d scaled_domain, domain_size;
    double total_energy{0}, temperature{0}, length{0}, strain{0};
    unsigned long long int seconds{0};
    int sizes{0};
    string padding = "\t0\t0\t0\t0\t0\t0\n";

    for(Atoms atoms : *V.atomsList) {
        NeighborList neighbors(V.cutoff_range);
        domain_size = (*V.domain_size)[sizes]; sizes++;
        length = domain_size[2];
        
        Domain domain(MPI_COMM_WORLD, domain_size, V.divisions, V.isEnabled);

        name = V.expName + "_parallel_" + to_string(atoms.nb_atoms()) 
                            + "_DLM_" + to_str(V.deltaL/length, 3) 
                            + "_T_"  + to_string(V.bt_targetT);
        // ofstream xyz_file(V.outLoc + name + ".xyz"); 
        ofstream csv_file(V.outLoc + name + ".csv");

        cout << "Running experiment with " << to_string(atoms.nb_atoms()) << 
                                        " atoms.\n";

        csv_file << "loop\ttotalEnergy\tpotEnergy\tkinEnergy\ttemperature\t"
                 << "strain\tleftGhostForce\trightGhostForce\n";

        auto begin = chrono::high_resolution_clock::now();
        domain.enable(atoms);
        for (int loop{0}; loop < V.loops; ++loop) {
            double lPotE{0}, lKinE{0}, gPotE{0}, gKinE{0}; 
            double llForce{0}, lrForce{0}, glForce{0}, grForce{0};
            
            if(domain.rank() == 0 && loop % 
                (V.loops / V.progress_cout_divisions) == 0) {
                cout << "Iteration: " << loop << "/" << V.loops << endl;
            }

            if(V.doStretch && loop != 0 && loop % V.stretch_loops == 0) {
                length += V.deltaL;
                scaled_domain << domain_size[0], domain_size[1], length;
                domain.scale(atoms, scaled_domain);
            }

            verletStep1(atoms, V.dt);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, 2 * V.cutoff_range);
            neighbors.update(atoms);
            gupta(atoms, neighbors);
            verletStep2(atoms, V.dt);

            if(V.bt_enable) bThermostat(atoms, V.bt_targetT, V.dt, V.bt_tau);

            // if(true) {
            //     int ghosts{0};
            //     for(int p{domain.nb_local()}; p < atoms.positions.cols(); ++p)
            //         ghosts += 1;
            //     if(loops % 20) {
            //         cout << ghosts << " ghosts in rank " << domain.rank() 
            //                 << "." << atoms.positions.cols() 
            //                 << "atoms present in this rank.\n";
            //     }
            // }

            bool csv_condition{loop % V.csv_loops == 0};
            if(csv_condition) {
                getKineticE(atoms);
                for (int i{0}; i < domain.nb_local(); ++i) {
                    lPotE += atoms.pot_energies(i);
                    lKinE += atoms.kin_energies(i);
                }
                gPotE = MPI::allreduce(lPotE, MPI_SUM, MPI_COMM_WORLD);
                gKinE = MPI::allreduce(lKinE, MPI_SUM, MPI_COMM_WORLD);
            }

            //if-statements made to consider ghost atoms of mid sub-region
            //  with using a rank of 3.
            if(domain.rank() == 0) 
                for(int left{domain.nb_local()}; 
                        left < atoms.positions.cols(); ++left)
                    if(atoms.positions(2,left) < 0)
                        llForce += atoms.forces(2,left);
                
            if(domain.rank() == domain.size() - 1) 
                for(int right{domain.nb_local()}; 
                        right < atoms.positions.cols(); ++right)
                    if(atoms.positions(2,right) > length)
                        lrForce += atoms.forces(2,right);

            glForce = MPI::allreduce(llForce, MPI_SUM, MPI_COMM_WORLD);
            grForce = MPI::allreduce(lrForce, MPI_SUM, MPI_COMM_WORLD);

            if(domain.rank() == 0 && csv_condition) {
                    total_energy = gPotE + gKinE;
                    temperature = getSystemT(gKinE, atoms.nb_atoms());
                    strain = (length - domain_size[2]) / domain_size[2];
                    csv_file << loop << "\t" << total_energy << "\t"
                             << gPotE << "\t" << gKinE << "\t" 
                             << temperature <<"\t" << strain << "\t" 
                             << glForce << "\t" << grForce << "\n";
                    // write_xyz(xyz_file, atoms);
            }
        }
        domain.disable(atoms);
        auto end = chrono::high_resolution_clock::now();

        elapsed.push_back(chrono::duration_cast<chrono::seconds>(end - begin));
        
        cout << "Execution time: " << elapsed.back().count() << " seconds\n";
        
        csv_file << "etime:"    << "\t" << elapsed.back().count() 
                                                            << padding;
        csv_file << "mass:"     << "\t" << atoms.mass       << padding;
        csv_file << "loops:"    << "\t" << V.loops          << padding;
        csv_file << "dt:"       << "\t" << V.dt             << padding;
        csv_file << "cutoff:"   << "\t" << V.cutoff_range   << padding;
        
        // xyz_file.close();
        csv_file.close();
    }
    seconds = chrono_time_sum(elapsed);
    cout << "Completed m9!\n";
    cout << "Total execution time in sec for milestone: " << seconds <<"\n\n";
    return seconds;
}