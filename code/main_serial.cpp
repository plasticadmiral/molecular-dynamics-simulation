#include <math.h>
#include "header/milestones_serial.h"

int main(int argc, char** argv) {
    cout << "Running m4 - m7 ...\n\n";

    unsigned long long int totalT{0};
    // double mass_SF = 103.6; //  g/mol
    double H_mass = 1.00794; // g/mol
    double Au_mass = 196.966657 * 103.6; // g/mol
    tuple<Names_t, Positions_t> positions;
    tuple<Names_t, Positions_t, Velocities_t> trajectory;
    string lj54Path{"../data/defaults/"}; 
    string cubesPath{"../data/cubes/"}; 
    string clustersPath{"../data/clusters/"}; 
    string defaultsPath{"../data/defaults/"};
    
    trajectory = read_xyz_with_velocities(lj54Path+"lj54.xyz");
    Atoms LJ54(get<0>(trajectory), get<1>(trajectory), get<2>(trajectory), 
                                                        H_mass);

    vector<Atoms> cubes;
    for (int i{11}; i < 19; ++i) {
        positions = read_xyz(cubesPath+"cube_"+to_string(i)+".xyz");
        cubes.push_back(Atoms(get<0>(positions), get<1>(positions), H_mass));
    }

    vector<Atoms> clusters;
    for (int i{13}; i < 14; ++i) {
        positions = read_xyz(clustersPath+"cluster_2.885_"+to_string(i)+".xyz");
        clusters.push_back(Atoms(get<0>(positions), get<1>(positions),
                                                                Au_mass));
    }
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    /*
    // m4: Energy conservation of lj54 structure simulated with LJ potential.
    vars4 V4;
    V4.expName = "lj54";
    V4.outLoc = "../results/m4/";
    V4.loops = 5000;
    V4.epsilon = 1;
    V4.sigma = 1;
    vector<double> dtVect{0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035};
    V4.dtList = &dtVect;

    Atoms v4LJ54 = LJ54;
    V4.atoms = &v4LJ54;
    
    totalT += m4(V4);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    

    // m4a: Potential energy graph between two atoms with lennard jones.
    vars4a V4a;
    V4a.expName = "lj2";
    V4a.outLoc = "../results/m4/";
    V4a.loops = 6000;
    V4a.dt = 0.001;
    vector<double> sigmaVect{2.5, 2, 1.5, 1};
    vector<double> epsilonVect{2.5, 2, 1.5, 1};
    V4a.sigmaList = &sigmaVect;
    V4a.epsilonList = &epsilonVect;
    Positions_t pos(3,2);
    Velocities_t vel(3,2);
    pos.setZero();
    vel.setZero();
    vector<Atoms> lj2atoms;
    for(double i{2}; i < 6; ++i) {
        cout << i/2 << endl;
        pos(0,1) = i/2 - 0.05;
        vel(0,1) = 0.5;
        lj2atoms.push_back(Atoms(pos, vel));
    }
    V4a.atoms = &lj2atoms;

    totalT += m4a(V4a);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";

    // m6: Performance with neighbour lists.
    vars6 V6;
    V6.expName = "ljCubeNL";
    V6.outLoc = "../results/m6/";
    V6.loops = 20;
    V6.epsilon = 1;
    V6.sigma = 1;
    V6.dt = 0.1;
    V6.cutoff_range = 5;

    vector<Atoms>m6Cubes = cubes;
    V6.atomsList = &m6Cubes;
    
    totalT += m6(V6);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    
    // m5: Performance without neighbour lists.
    vars5 V5;
    V5.expName = "ljCube";
    V5.outLoc = "../results/m5/";
    V5.loops = 20;
    V5.epsilon = 1;
    V5.sigma = 1;
    V5.dt = 0.1;

    vector<Atoms>m5Cubes = cubes;
    V5.atomsList = &m5Cubes;

    totalT += m5(V5);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";

    // m5a: Cooling down and heating up a cube with berendsen thermostat 
    //      to make a gif.
    vars5a V5a;
    V5a.expName = "ljCubegif";
    V5a.outLoc = "../results/m5/";
    V5a.loops = 2000;
    V5a.heatT = 5000;
    V5a.coolT = 2500;
    V5a.dt = 0.01;
    V5a.tau = 3;

    positions = read_xyz(cubesPath+"cube_05.xyz");
    Atoms v5aCube = Atoms(get<0>(positions), get<1>(positions), H_mass);
    V5a.atoms = &v5aCube;

    totalT += m5a(V5a);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    */

    // m7: Heating up gold cluster with velocity rescaling to molten state 
    //     simulated with gupta potential.
    vars7 V7;
    V7.expName = "gptaCluster";
    V7.outLoc = "../results/m7/";
    V7.loops = 100000;
    V7.csv_loops = 8;
    V7.rel_loops = 2000;
    V7.dt = 2;
    V7.cutoff_range = 5.0;

    vector<int> deltaQVect;
    vector<Atoms>m7Clusters;
    for (Atoms cluster : clusters) {
        m7Clusters.push_back(cluster);
        deltaQVect.push_back(int(cluster.nb_atoms()/100));
    }

    V7.deltaQ = &deltaQVect;
    V7.atomsList = &m7Clusters;
    
    totalT += m7(V7);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    
    /* 
    // m7a: Potential energy graph between two atoms with gupta.
    vars7a V7a;
    V7.expName = "gpta2";
    V7.outLoc = "../results/m7a/";
    V7a.loops = 2000;
    V7a.dt = 0.01;

    Positions_t posgpta(3,2);
    Velocities_t velgpta(3,2);
    posgpta << 0, 1e-5, 
               0, 0, 
               0, 0; 

    velgpta << 0, 0.1, 
               0, 0, 
               0, 0; 

    Atoms v7aAtoms(posgpta, velgpta, Au_mass);
    V7a.atoms = &v7aAtoms;
    
    totalT += m7a(V7a);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    
    // m7b: Heating up gold cluster with velocity rescaling to molten state 
    //      to make a gif.
    vars7b V7b;
    // V7b.expName = "gptaCluster-gif";
    V7b.expName = "gptaCluster";
    V7b.outLoc = "../results/m7/";
    V7b.loops = 7000;
    // V7b.deltaQ = 19;
    // V7b.rel_loops = 2000;
    V7b.dt = 2;
    V7b.cutoff_range = 5.0;
    positions = read_xyz(defaultsPath+"cluster_923.xyz");
    Atoms m7bCluster(get<0>(positions), get<1>(positions), Au_mass);
    V7b.atoms = &m7bCluster;

    totalT += m7b(V7b);
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    */ 

    unsigned long long int minutes{totalT/60};
    unsigned long long int seconds{(totalT/60 - minutes) * 60}; 
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
    cout << "Total execution time of code: " << minutes << " minutes and " 
                            << seconds << " seconds\n"; 
    cout << "-----x-----x-----x-----x-----x-----x-----x-----x-----x-----x\n";
}