#ifndef __MILESTONES_PARALLEL_H
#define __MILESTONES_PARALLEL_H

#include <vector>
#include <chrono>
#include <string>
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
#include "mpi_support.h"
#include "domain.h"

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
    vector<Atoms>* atomsList;
    string expName{"default"};
    string outLoc{"../results/"};
    unsigned int loops;
    unsigned int csv_loops{1};
    unsigned int progress_cout_divisions{20};
    Eigen::Array3i divisions; // number of subdomains in each axis
    Eigen::Array3i isEnabled; // 1 or 0
    double dt{2}; // femtoseconds
    double cutoff_range{5.0}; // Armstrong 

};

/**
 * @brief A struct used to pass compact arguments to m8 fn.
 */
struct vars8 : commonVars { //cluster
    Eigen::Array3d domain_size; // in Armstrong
};

/**
 * @brief A struct used to pass compact arguments to m9 fn.
 */
struct vars9 : commonVars { //cluster
    vector<Eigen::Array3d>* domain_size; // in Armstrong
    bool doStretch; // stretching using domain fn.
    double deltaL; // lenght to change every after every stretch_steps
    int stretch_loops; // time interval after which to strech with domain fn.
    bool bt_enable; // heat using velocity scaling.
    int bt_targetT; // temperature to be sustained by berendsen thermostat.
    int bt_tau; // tau value used by berendsen thermostat.
};

/**
 * @brief       parallelizing m7
 * 
 * @param V     A struct holding arguments specifically for m8.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m8(vars8& V);

/**
 * @brief       stretching gold wire
 * 
 * @param V     A struct holding arguements specifically for m9.
 * @return      Time taken to complete in sec.
 */
unsigned long long int m9(vars9& V);

#endif  // _MILESTONES_PARALLEL_H