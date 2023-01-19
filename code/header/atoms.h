#ifndef __ATOMS_H
#define __ATOMS_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>

using Names_t = std::vector<std::string>;
using Positions_t = Eigen::Array3Xd; 
using Velocities_t = Eigen::Array3Xd; 
using Forces_t = Eigen::Array3Xd; 
using Energy_t = Eigen::ArrayXd; 

/**
 * @brief   A struct to manage Atoms in a system.
 */
struct Atoms { 
    double mass;
    Names_t names;
    Positions_t positions; 
    Velocities_t velocities; 
    Forces_t forces; 
    Energy_t pot_energies; 
    Energy_t kin_energies; 
    
    Atoms(const int nb_atoms) : 
            positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, 
            pot_energies{nb_atoms}, kin_energies{nb_atoms} {
        
        for(int i=0; i<nb_atoms; ++i) names.push_back("H");
        mass = 1.00784;
        velocities.setZero();
        forces.setZero();
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Positions_t &p) : 
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        for(int i=0; i<p.cols(); ++i) names.push_back("H");
        mass = 1.00784;
        velocities.setZero(); 
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Positions_t &p, const Velocities_t &v) : 
            positions{p}, velocities{v}, forces{3, p.cols()}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        assert(p.cols() == v.cols());
        for(int i=0; i<p.cols(); ++i) names.push_back("H");
        mass = 1.00784;
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Positions_t &p, const double m) : 
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        for(int i=0; i<p.cols(); ++i) names.push_back("X");
        mass = m;
        velocities.setZero(); 
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Positions_t &p, const Velocities_t &v, const double m) : 
            positions{p}, velocities{v}, forces{3, p.cols()}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        assert(p.cols() == v.cols());
        for(int i=0; i<p.cols(); ++i) names.push_back("X");
        mass = m;
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Positions_t &p, const Velocities_t &v, const Forces_t &f, 
            const double m) : positions{p}, velocities{v}, forces{f}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        assert(p.cols() == v.cols());
        assert(v.cols() == f.cols());
        for(int i=0; i<p.cols(); ++i) names.push_back("X");
        mass = m; 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Names_t &n, const Positions_t &p) : 
            names{n}, positions{p}, velocities{3, p.cols()}, 
            forces{3, p.cols()}, pot_energies{p.cols()}, 
            kin_energies{p.cols()} { 

        mass = 1;
        velocities.setZero(); 
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v) : 
            names{n}, positions{p}, velocities{v}, forces{3, p.cols()}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        mass = 1;
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Names_t &n, const Positions_t &p, const double m) : 
            names{n}, positions{p}, velocities{3, p.cols()}, 
            forces{3, p.cols()}, pot_energies{p.cols()}, 
            kin_energies{p.cols()} { 

        mass = m;
        velocities.setZero(); 
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v, 
            const double m) : 
            names{n}, positions{p}, velocities{v}, forces{3, p.cols()}, 
            pot_energies{p.cols()}, kin_energies{p.cols()} { 

        mass = m;
        forces.setZero(); 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v, 
            const Forces_t &f, const double m) : names{n}, positions{p}, 
            velocities{v}, forces{f}, pot_energies{p.cols()}, 
            kin_energies{p.cols()} { 

        assert(p.cols() == v.cols());
        assert(v.cols() == f.cols());
        mass = m; 
        pot_energies.setZero();
        kin_energies.setZero();
    } 
    void resize(const int size){
        positions.conservativeResize(3, size);
        velocities.conservativeResize(3, size);
        forces.conservativeResize(3, size);
        pot_energies.conservativeResize(size);
        kin_energies.conservativeResize(size);
    }
    long int nb_atoms() const { 
        return positions.cols();         
    }
};


#endif  // __ATOMS_H
