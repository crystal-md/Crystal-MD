//
// Created by genshen on 2019-06-10.
//

#include <logs/logs.h>
#include <iostream>
#include <utils/mpi_domain.h>

#include "md_simulation.h"
#include "system_configuration.h"
#include "io/output_dump.h"
#include "io/output_copy.h"

MDSimulation::MDSimulation(ConfigValues *p_config_values)
        : simulation(), cur_stage_steps(0), pConfigVal(p_config_values), dump_instances() {
    // initialize stage, we make a stage with zero steps.
    // then, in first time step simulation, it will move to next stage.
    Stage s;
    s.steps = 0;
    current_stage = s;
    next_stage_index = 0;
}

void MDSimulation::onSimulationStarted(const unsigned long init_step) {
    for (const DumpConfig &preset: pConfigVal->output.presets) {
        switch (preset.mode) {
            case OutputMode::DEBUG:
                this->dump_instances[preset.name] = new OutputDump(preset, *_p_domain);
                break;
            case OutputMode::COPY:
                this->dump_instances[preset.name] = new OutputCopy(preset, *_p_domain);
                break;
        }
    }
    if (init_step != 0) {
        moveToStageAtStep(init_step);
    }
}

void MDSimulation::moveToStageAtStep(const unsigned long step) {
    // test from the first stage.
    unsigned long step_stage_begin = 0;
    unsigned long step_stage_end = 0;
    std::size_t i = 0;
    bool flag = false;
    for (const auto &stage: pConfigVal->stages) {
        step_stage_end += stage.steps;
        if (step >= step_stage_begin && step < step_stage_end) {
            // set stage status
            setCurStageStatus(stage, step - step_stage_begin);
            next_stage_index = i + 1; // prepare for moving to next stage.
            flag = true;
            break;
        }
        step_stage_begin = step_stage_end;
        i++;
    }
    if (!flag) {
        kiwi::logs::e("simulation", "step offset too large when moving to new stage.\n");
        abort(1);
    }
}

void MDSimulation::onSimulationDone(const unsigned long step) {
    for (auto ins_pair :this->dump_instances) {
        ins_pair.second->onAllOut(step);
        delete ins_pair.second;
    }
}

void MDSimulation::beforeStep(const unsigned short potentialType,const unsigned long step) {
    // set new stage
    if (cur_stage_steps >= current_stage.steps) {
        // if we are out of stage steps, then we can move to next stage.
        if (next_stage_index < pConfigVal->stages.size()) {
            // move to next stage.
            setCurStageStatus(pConfigVal->stages[next_stage_index], 0);
            next_stage_index++;
        } else {
            kiwi::logs::e("simulation", "no more stages.\n");
            abort(1);
        }
    }

    kiwi::logs::s(MASTER_PROCESSOR, "simulation", "simulating steps: {}/{}\r",
                  step + 1, pConfigVal->timeSteps);

    // perform collision.
    if (current_stage.collision_set && cur_stage_steps + 1 == current_stage.collisionStep) {
        // just output atoms in preview step of the collision step.
        collisionStep(potentialType, step, current_stage.collisionLat, current_stage.direction, current_stage.pkaEnergy);
    }

    // perform "velocity"
    if (current_stage.velocity_set && cur_stage_steps + 1 == current_stage.velocity_step) {
        // set atoms' velocity in a region
        const comm::Region<long> region(current_stage.velocity_region[0], current_stage.velocity_region[1],
                                        current_stage.velocity_region[2], current_stage.velocity_region[3],
                                        current_stage.velocity_region[4], current_stage.velocity_region[5]);
        velocitySetStep(region, current_stage.velocity_value);
    }

    // perform rescale
    if (current_stage.rescales_set && cur_stage_steps % current_stage.rescale_every == 0) {
        const _type_atom_count n_global_atoms = 2 * _p_domain->phase_space[0] *
                                                _p_domain->phase_space[1] *
                                                _p_domain->phase_space[2];
        configuration::rescale(current_stage.rescale_t, n_global_atoms, _atom->atom_list, _atom->inter_atom_list);
    }
    // log thermodynamic
    // todo: use the same (cur_stage_steps +1)
    if (current_stage.thermo_logs_set && (cur_stage_steps + 1) % current_stage.thermo_logs_every_steps == 0) {
        runtime_status.flag_calc_system_potential_energy = true;
    } else {
        runtime_status.flag_calc_system_potential_energy = false;
    }
}

void MDSimulation::postStep(const unsigned long step) {
    phy_time += current_stage.step_length; // the physical time when this iteration is finished.
    cur_stage_steps++; // add current stage steps.
    kiwi::logs::i(MASTER_PROCESSOR, "simulation", "simulated physical time: {} ps.\n", phy_time);

    // output atoms information if it is the dumping step
    if (current_stage.dump_set && cur_stage_steps % current_stage.dump_every_steps == 0) {
        if (this->dump_instances.find(current_stage.dump_preset_use) == this->dump_instances.end()) {
            kiwi::logs::e("output", "bad output dump `use`: {}.\n", current_stage.dump_preset_use);
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
            OutputBaseInterface *dump_ptr = this->dump_instances[current_stage.dump_preset_use];
            dump_ptr->onOutputStep(step + 1, _atom->getAtomList(), _atom->getInterList());
        }
    }

    // output thermodynamics information if it is the step
    log_thermodynamics(current_stage, cur_stage_steps, step, phy_time);

#ifdef MD_RUNTIME_CHECKING
    {
        unsigned long count[2] = {0, 0};
        count[0] = _atom->realAtoms();
        count[1] = _atom->getInterList()->nlocalinter;
        unsigned long global_count[2] = {0, 0};
        unsigned long &real_atoms = global_count[0], &inter_atoms = global_count[1];
        MPI_Reduce(count, global_count, 2, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_PROCESSOR, MPI_COMM_WORLD);

        kiwi::logs::d("count", "real:{}--inter: {}\n", count[0], count[1]);
        kiwi::logs::d(MASTER_PROCESSOR, "count", "global_real:{}--global_inter: {}\n", real_atoms, inter_atoms);
    }
#endif
}

void MDSimulation::log_thermodynamics(const Stage stage, const unsigned long cur_stage_step,
                                      const unsigned long global_step, const double phy_time) {
    if (current_stage.thermo_logs_set && cur_stage_steps % current_stage.thermo_logs_every_steps == 0) {
        // find by `use`
        md_thermodynamic::OutputThermodynamic preset;
        for (md_thermodynamic::OutputThermodynamic p: pConfigVal->output.thermo_presets) {
            if (p.name == current_stage.thermo_logs_preset_use) {
                preset = p;
                break;
            }
        }

        std::stringstream ss;
        // time and step:
        if ((preset.flags & md_thermodynamic::WithStepMask) != 0) {
            ss << " step: " << global_step;
        }
        if ((preset.flags & md_thermodynamic::WithTimeMask) != 0) {
            ss << " phy_time: " << phy_time;
        }

        if ((preset.flags & md_thermodynamic::WithKineticEnergyMask) != 0 ||
            (preset.flags & md_thermodynamic::WithTemperatureMask) != 0) {
            const double ke = configuration::kineticEnergy(_atom->getAtomList(), _atom->getInterList(),
                                                           configuration::ReturnMod::Root, MASTER_PROCESSOR);
            if ((preset.flags & md_thermodynamic::WithKineticEnergyMask) != 0) {
                ss << " kinetic energy: " << ke;
            }
            if ((preset.flags & md_thermodynamic::WithTemperatureMask) != 0) {
                const _type_atom_count n = 2 * pConfigVal->phaseSpace[0] *
                                           pConfigVal->phaseSpace[1] * pConfigVal->phaseSpace[2];
                const double T = configuration::temperature(ke, n);
                ss << " temperature: " << T;
            }
        }
        if ((preset.flags & md_thermodynamic::WithPotentialEnergyMask) != 0) {
            double pot_energy = _atom->get_system_pot_energy();
            ss << " potential energy: " << pot_energy;
        }
        // write buffer
        kiwi::logs::i(MASTER_PROCESSOR, "thermo", "{}\n", ss.str());
    }
}

void MDSimulation::onForceSolved(const unsigned long step) {
#ifdef MD_RUNTIME_CHECKING
    forceChecking();
#endif
}

void MDSimulation::print_force(const std::string filename, int step) {
    std::ofstream outfile;
    outfile.open(filename);

    _atom->atom_list->foreachSubBoxAtom(
            [&](const _type_atom_index gid) {
                MD_LOAD_ATOM_VAR(_atom_ref, _atom->atom_list, gid);
                outfile << MD_GET_ATOM_ID(_atom_ref, gid) << "\t"
                        << MD_GET_ATOM_X(_atom_ref, gid, 0) << "\t"
                        << MD_GET_ATOM_X(_atom_ref, gid, 1) << "\t"
                        << MD_GET_ATOM_X(_atom_ref, gid, 2) << "\t"
                        << MD_GET_ATOM_V(_atom_ref, gid, 0) << "\t"
                        << MD_GET_ATOM_V(_atom_ref, gid, 1) << "\t"
                        << MD_GET_ATOM_V(_atom_ref, gid, 2) << "\t"
                        << MD_GET_ATOM_F(_atom_ref, gid, 0) << "\t"
                        << MD_GET_ATOM_F(_atom_ref, gid, 1) << "\t"
                        << MD_GET_ATOM_F(_atom_ref, gid, 2)
                        << std::endl;
            }
    );
    outfile << "inter atoms" << std::endl;
    for (_type_inter_list::iterator inter_it = _atom->inter_atom_list->inter_list.begin();
         inter_it != _atom->inter_atom_list->inter_list.end(); inter_it++) {
        AtomElement &_atom_ref = *inter_it;
        outfile << std::endl << _atom_ref.id << "\t"
                << _atom_ref.x[0] << "\t" << _atom_ref.x[1] << "\t" << _atom_ref.x[2] << "\t"
                << _atom_ref.v[0] << "\t" << _atom_ref.v[1] << "\t" << _atom_ref.v[2] << "\t"
                << _atom_ref.f[0] << "\t" << _atom_ref.f[1] << "\t" << _atom_ref.f[2] << std::endl;
    }
    outfile.close();
}

void MDSimulation::setCurStageStatus(const Stage stage, const unsigned long stage_step) {
    current_stage = stage;
    cur_stage_steps = stage_step;
    // set new step length
    if (stage.step_length == 0.0) { // if it is 0, set a default value.
        current_stage.step_length = pConfigVal->timeStepLength;
    }
    // reset time steps length
    _newton_motion->setTimestepLength(current_stage.step_length);
}

#ifdef MD_RUNTIME_CHECKING

void MDSimulation::forceChecking() {
    auto forces = configuration::systemForce(_atom->getAtomList(), _atom->getInterList());
    double fx[3] = {forces[0], forces[1], forces[2]};
    double fx_2[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(fx, fx_2, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (std::abs(fx_2[0]) > 0.00001 || std::abs(fx_2[1]) > 0.00001 || std::abs(fx_2[2]) > 0.00001) {
        kiwi::logs::d("force2", "({}, {}, {})\n\n", forces[0], forces[1], forces[2]);
        kiwi::logs::d(MASTER_PROCESSOR, "sum force:", "({}, {}, {})\n\n", fx_2[0], fx_2[1], fx_2[2]);
        char filename[20];
        sprintf(filename, "force_%d.txt", kiwi::mpiUtils::global_process.own_rank);
        print_force(filename, _simulation_time_step);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
}

#endif //MD_RUNTIME_CHECKING
