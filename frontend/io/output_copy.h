//
// Created by genshen on 2019-07-29.
//

#ifndef MISA_MD_OUTPUT_COPY_H
#define MISA_MD_OUTPUT_COPY_H

#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "../config_values.h"
#include "atom_dump.h"
#include "output_base_interface.h"

/**
 * @brief abstract class of simulation step callback.
 */
class OutputCopy : public OutputBaseInterface {
public:
    explicit OutputCopy(const DumpConfig output, const comm::BccDomain p_domain);

    void prepareOutput(const comm::BccDomain p_domain) override;

    /**
     * this will be call in each output step.
     * @param time_step current time step
     * @param atom_list list of lattice atoms.
     * @param inter_atom_list list of inter atoms.
     */
    void onOutputStep(const unsigned long time_step, AtomList *atom_list, InterAtomList *inter_atom_list) override;

    /**
     * this will be call when all time steps finished.
     * @param time_step current time step
     */
    void onAllOut(const unsigned long time_step) override;

protected:
    /**
     * dump instance if it is not outputting by frame (only one file shared by all processors).
     */
    AtomDump *dumpInstance = nullptr;

    /**
     * time of outputting
     */
    double totalDumpTime = 0;
};

#endif //MISA_MD_OUTPUT_COPY_H
