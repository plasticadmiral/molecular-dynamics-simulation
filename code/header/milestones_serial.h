#ifndef __MILESTONES_SERIAL_H
#define __MILESTONES_SERIAL_H

#include <vector>
#include <chrono>
#include <string>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include "atoms.h"
#include "verlet.h"
#include "xyz.h"
#include "ljDS.h"
#include "kineticAndTemp.h"
#include "lj.h"
#include "neighbors.h"
#include "bThermostat.h"
#include "gupta.h"

using namespace std;

template <typename T>
std::string to_str(const T a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

template<typename T>
inline void print_vectr(vector<T> const &vectr) {
    for(T val : vectr) {
        cout << val << " ";
    }
    cout << endl;
}

inline unsigned int chrono_time_sum(vector<chrono::seconds> vectr) {
    unsigned int seconds{0};
    for(chrono::seconds val : vectr) {
        seconds += val.count();
    }
    return seconds;
}


/**
 * @brief A struct containing common arguments.
 */
struct commonVars {
    string expName{"default"};
    string outLoc{"../results/"};
    unsigned long long int loops{2000}; // Number of for-loop iterations.
    unsigned int csv_loops{1}; // Data storage timestep multiplier
    unsigned int progress_cout_divisions{20}; // if(step%(step/??)==0) "cout step"
};

/**
 * @brief A struct used to pass compact arguments to m4 fn.
 */
struct vars4 : commonVars {
    Atoms* atoms;
    double epsilon{1}, sigma{1};
    vector<double>* dtList;
};

/**
 * @brief A struct used to pass compact arguments to m4a fn.
 */
struct vars4a : commonVars { 
    vector<Atoms>* atoms;
    vector<double>* epsilonList;
    vector<double>* sigmaList;
    double dt{0.1};
};

/**
 * @brief A struct used to pass compact arguments to m5 fn.
 */
struct vars5 : commonVars {
    vector<Atoms>* atomsList;
    double epsilon{1}, sigma{1};
    double dt{0.01};
};

/**
 * @brief A struct used to pass compact arguments to m5a fn.
 */
struct vars5a : commonVars {
    Atoms* atoms;
    double epsilon{1}, sigma{1};
    double dt{0.01};
    double tau{3};
    int heatT, coolT;
};

/**
 * @brief A struct used to pass compact arguments to m6 fn.
 */
struct vars6 : commonVars { //cubicLatticeNL
    vector<Atoms>* atomsList;
    double epsilon{1}, sigma{1};
    double dt{0.01};
    double cutoff_range{5.0};
};

/**
 * @brief A struct used to pass compact arguments to m7 fn.
 */
struct vars7 : commonVars { //cluster
    vector<Atoms>* atomsList;
    vector<int>* deltaQ;
    int rel_loops;
    double dt{2};
    double cutoff_range{5.0};
};

/**
 * @brief A struct used to pass compact arguments to m7a fn.
 */
struct vars7a : commonVars { 
    string expName{"gupta_GraphData"};
    string outLoc{"../results/m7a/"};
    Atoms* atoms;
    double dt{0.01};
    double cutoff_range{5.0};
};

/**
 * @brief A struct used to pass compact arguments to m7b fn.
 */
struct vars7b : commonVars { 
    Atoms* atoms;
    double deltaQ;
    int rel_loops;
    double dt{2};
    double cutoff_range{5.0};
};

/**
 * @brief       Simulates the lj54.xyz file with Lennard Jones potential for 
 *              different dt. 
 *              
 *              "total energy as a function of steps" for different dt 
 *              is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m4.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m4(const vars4& V);

/**
 * @brief       Moves two atoms away from each other at constant velocity 
 *              while calculating the lennard jones potential between them.
 *              
 *              "Potential energy vs distance" is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m4a.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m4a(const vars4a& V);


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
unsigned long long int m5(const vars5& V);

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
unsigned long long int m5a(const vars5a& V);


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
unsigned long long int m6(const vars6& V);


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
unsigned long long int m7(const vars7& V);

/**
 * @brief       Moves two atoms away from each other at constant velocity 
 *              while calculating the gupta potential between them.
 * 
 *              "Potential energy vs distance" is plotted with this data.
 * 
 * @param V     A struct holding arguments specifically for m7a.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m7a(const vars7a& V);

/**
 * @brief       Performs the simulation from m7 on a single cluster but with 
 *              extra heat introduction and longer relaxation periods to 
 *              make the process more visible.
 *
 *              A gif is made using this data. 
 * 
 * @param V     A struct holding arguments specifically for m7b.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m7b(const vars7b& V);


#endif  // _MILESTONES_SERIAL_H