#include "header/milestones_parallel.h"

using namespace std;

int main(int argc, char *argv[]) {
    cout << "Running m8 and m9...\n";

    unsigned long long int totalT{0};
    double Au_mass = 196.966657 * 103.6; //g/mol
    tuple<Names_t, Positions_t> positions;
    // string clustersPath{"../data/clusters/"}; 
    string clustersPath{"../data/defaults/"}; 
    string whiskersPath{"../data/whiskers/"};

    // vector<Atoms> clusters;
    // for (int i{1}; i < 13; ++i) {
    //     positions = read_xyz(clustersPath+"cluster_"+to_string(i)+"_500.xyz");
    //     Atoms atoms(get<0>(positions), get<1>(positions));
    //     clusters.push_back(atoms);
    // }

    positions = read_xyz(clustersPath + "cluster_923.xyz");
    Atoms cluster = Atoms(get<0>(positions), get<1>(positions), Au_mass);

    vector<Atoms> whiskers;
    for (int i{3}; i < 4; ++i) {
        positions = read_xyz(whiskersPath+"whisker"+to_string(i)+".xyz");
        Atoms atoms(get<0>(positions), get<1>(positions), Au_mass);
        whiskers.push_back(atoms);
    }
    MPI_Init(&argc, &argv);
    /*
    // m8: Parallelized version of m7.
    vars8 V8;
    V8.expName = "cluster";
    V8.outLoc = "../results/m8/";

    vector<Atoms> clusterVect{cluster};
    V8.atomsList = &clusterVect;

    V8.loops = 20000;
    V8.dt = 2;  // femtoseconds
    V8.cutoff_range = 5;
    V8.domain_size << 35, 35, 35;
    V8.divisions << 2, 2, 2;
    V8.isEnabled << 1, 1, 1;


    totalT += m8(V8);
    */
    
    //m9: Stretching gold whiskers.
    vars9 V9;
    V9.progress_cout_divisions = 50;
    V9.csv_loops = 8;
    V9.expName = "whisker";
    V9.outLoc = "../results/m9/";
    V9.loops = 500000;
    V9.stretch_loops = 20000;
    V9.dt = 2; // femtoseconds
    V9.cutoff_range = 5;
    V9.atomsList = &whiskers;
    V9.divisions << 1, 1, 3;
    V9.isEnabled << 0, 0, 1;
    V9.bt_enable = true;
    V9.bt_tau = 3;

    Eigen::Array3d domain;
    vector<Eigen::Array3d> domainVect;
    domain << 16, 22, 57.70;
    domainVect.push_back(domain);
    domain << 25, 34, 57.70;
    domainVect.push_back(domain);
    domain << 31, 42, 57.70;
    domainVect.push_back(domain);
    V9.domain_size = &domainVect;
        
    
    int temperatures[3] = {0, 3, 5}; // kelvin
    double deltaLs[3] = {57.70 * 0.010, 57.70 * 0.015, 57.70 * 0.018};
    // double deltaLs[3] = {57.70 * 0.010, 57.70 * 0.0125, 57.70 * 0.015};
    // int temperatures[3] = {0, 3, 5}; // in kelvin
    // double deltaLs[3] = {57.70 * 0.005, 57.70 * 0.008, 57.70 * 0.010};
    for(int temperature : temperatures) {
        for(double deltaL : deltaLs) {
            V9.deltaL = deltaL; 
            V9.bt_targetT = temperature;
            totalT += m9(V9);
        }
    }
    
    MPI_Finalize();
    unsigned long long int minutes{totalT/60};
    unsigned long long int seconds{(totalT/60 - minutes) * 60}; 
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    cout << "Total execution time of code: " << minutes << " minutes and " 
                            << seconds << " seconds\n"; 
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";

}