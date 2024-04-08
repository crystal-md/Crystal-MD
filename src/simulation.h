//
// Created by baihe in 2016-12-22.
//

#ifndef MISA_MD_SIMULATION_H
#define MISA_MD_SIMULATION_H

#include <mpi.h>
#include <cstring>
#include <io/io_writer.h>
#include <eam.h>

#include "newton_motion.h"
#include "world_builder.h"

struct RuntimeStatus {
    /**
     * this flag indicates whether to calculate the system potential energy.
     * In current implementation, it is reset to true/false in function beforeStep at each step.
     */
    bool flag_calc_system_potential_energy = false;
};

class simulation {
public:

    simulation();

    virtual ~simulation();

    /**
     * Denote N as the count of all processors.
     * {@memberof domainDecomposition} will divide the simulation box into N parts,
     * we call each part as a sub-box.
     * And each sub-box will bind to a processor.
     * @param phase_space the lattice count in each dimension.
     * @param lattice_const lattice constance
     * @param cutoff_radius_factor cutoff radius factor = cutoff/lattice_const
     */
    void createDomain(const int64_t phase_space[DIMENSION],
                      const double lattice_const, const double cutoff_radius_factor);

    /**
     * initialize atoms position and velocity by given arguments.
     * @param phase_space the lattice count in each dimension.
     * @param lattice_const lattice const.
     * @param init_step_len time step length as initialized value.
     * @param create_mode whether create atom randomly.
     * @param create_seed seed to create atoms.
     * @param t_set initial temperature.
     * @param types types defines ratio and relative atomic mass of alloy for each types of material.
     * @param read_inp_path the path for reading atoms from.
     */
    void createAtoms(const int64_t phase_space[DIMENSION], const double lattice_const, const double init_step_len,
                     const bool create_mode, const double t_set, const unsigned long create_seed,
                     const std::vector<tp_atom_type_weight> &types_weight,
                     const std::string read_inp_path);

    /**
     * initialize potential function and perform the first simulation step.
     * @param potentialType type of potentail calculation.
     * @param pot_file_path path of potential file.
     */
    void prepareForStart(const unsigned short potentialType,const std::string pot_file_path);

    /**
     * set collision energy.
     * @param potentialType type of potentail calculation.
     * @param coll_step step to perform collision.
     * @param coll_lat lattice position of PKA.
     * @param coll_dir direction of collision.
     * @param coll_pka_energy pka energy, unit eV.
     */
    void collisionStep(const unsigned short potentialType,unsigned long coll_step, const _type_lattice_coord coll_lat[DIMENSION + 1],
                  const double coll_dir[DIMENSION], const double coll_pka_energy);

    /**
     * do time steps loop simulation.
     * @param potentialType type of potentail calculation.
     * @param steps total simulation steps.
     * @param init_step the initial step for step iterating.
     */
    void simulate(const unsigned short potentialType,const unsigned long steps, const unsigned long init_step);

    void finalize();

    /**
     * this function will be called before simulation loop.
     * @param init_step the initial step before simulation.
     * This is usually used in restart simulation mode (otherwise it is 0).
     */
    virtual void onSimulationStarted(const unsigned long init_step) {};

    /**
     * this function will be called after simulation loop finished.
     * @param step total simulation step.
     */
    virtual void onSimulationDone(const unsigned long step) {};

    /**
     * this callback function will be called before a simulation step.
     * @param potentialType type of potentail calculation.
     * @param step current simulation step, starting from 0.
     */
    virtual void beforeStep(const unsigned short potentialType,const unsigned long step) {};

    /**
     * this callback function will be called after a simulation step.
     * @param step current simulation step, starting from 0.
     */
    virtual void postStep(const unsigned long step) {};

    /**
     * this callback function will be called after atoms' force are computed.
     * @param step current simulation step, starting from 0.
     */
    virtual void onForceSolved(const unsigned long step) {};

    void abort(int exitcode);

protected:
    /**
     * the time steps the program have simulated.
     */
    unsigned long _simulation_time_step;

    /**
     * the runtime status.
     * It can be modified by plugin or in BeforeStep/PostStep callback.
     */
    RuntimeStatus runtime_status;

    comm::BccDomain *_p_domain;
    // GlobalDomain *p_domain;  //仅rank==0的进程有效 // todo ??
    atom *_atom;
    NewtonMotion *_newton_motion;

    eam *_pot; // eam potential

    /**
     * set velocity for atoms in a given region
     * @param global_region the given region
     * @param velocity_value velocity value at x,y,z dimension. unit A/ps
     */
    void velocitySetStep(const comm::Region<long> global_region, const double velocity_value[DIMENSION]);
};

#endif //MISA_MD_SIMULATION_H
