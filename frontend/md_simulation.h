//
// Created by genshen on 2019-06-10.
//

#ifndef MISA_MD_MD_SIMULATION_H
#define MISA_MD_MD_SIMULATION_H


#include <simulation.h>
#include "io/output_base_interface.h"

class MDSimulation : public simulation {
public:
    explicit MDSimulation(ConfigValues *p_config_values);

    /**
    * this function will be called before simulation loop.
    */
    void onSimulationStarted(const unsigned long init_step) override;

    /**
     * @note: this function must be called before simulation starting.
     */
    void moveToStageAtStep(const unsigned long step);

    /**
     * this function will be called after simulation loop finished.
     */
    void onSimulationDone(const unsigned long step) override;

    /**
     * this callback function will be called before a simulation step.
     * @param potentialType type of potentail calculation.
     * @param step current simulation step, starting from 0.
     */
    void beforeStep(const unsigned short potentialType, const unsigned long step) override;

    /**
     * this callback function will be called after a simulation step.
     * @param step current simulation step, starting from 0.
     */
    void postStep(const unsigned long step) override;

    /**
     * this callback function will be called after atoms' force are computed.
     * @param step current simulation step, starting from 0.
     */
    void onForceSolved(const unsigned long step) override;

private:
    /**
     * pointer to current stage in config.
     */
    Stage current_stage;
    std::size_t next_stage_index;

    /**
     * steps in current stage.
     */
    unsigned long cur_stage_steps;

    /**
     * physical time
     */
    double phy_time = 0.0;
    /**
     * pointer to config data.
     */
    const ConfigValues *pConfigVal = nullptr;

    // map of dump instance, indexed by dump preset name
    std::map<std::string, OutputBaseInterface *> dump_instances;

    /**
     * log force of all atoms in simulation system.
     * @param filename filename to save atoms' force
     * @param step current simulation step, starting from 0.
     */
    void print_force(const std::string filename, int step);

    /**
     * set the current stage status to the stage specified by @param stage
     * @param stage the stage
     * @param stage_step step offset in the new stage.
     */
    void setCurStageStatus(const Stage stage, const unsigned long stage_step);

#ifdef MD_RUNTIME_CHECKING

    /**
     * check force, only works in dev mode.
     * if the sum of all atoms not 0, program will exists with error message.
     */
    void forceChecking();

#endif

    /**
     * log thermodynamics information to standard output if it is the step
     * @param stage current stage
     * @param cur_stage_step current stage step
     * @param global_step current global step
     * @param phy_time physics time
     */
    void log_thermodynamics(const Stage stage, const unsigned long cur_stage_step,
                            const unsigned long global_step, const double phy_time);
};


#endif //MISA_MD_MD_SIMULATION_H
