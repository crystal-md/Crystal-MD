//
// Created by genshen(genshenchu@gmail.com) on 2018-3-9.
//

#ifndef MISA_MD_CONFIG_VALUES_H
#define MISA_MD_CONFIG_VALUES_H

#include <string>
#include <vector>
#include <array>
#include <map>

#include <utils/bundle.h>
#include "io/atom_dump_types.h"
#include "types/pre_define.h"
#include "types/atom_types.h"
#include "config_values_thermodynamic.hpp"

enum OutputMode {
    DEBUG = 0,
    COPY = 1,
};

#define LOGS_MODE_CONSOLE 0
#define LOGS_MODE_FILE 1
#define LOGS_MODE_CONSOLE_STRING "console"
#define LOGS_MODE_FILE_STRING "file"
#define DEFAULT_LOGS_MODE_CONSOLE_STRING LOGS_MODE_CONSOLE_STRING
#define DEFAULT_OUTPUT_DUMP_FILE_PATH "misa_mdl.out"
#define ORIGIN_OUTPUT_DUMP_FILE_PATH "origin_misa_mdl.out"

constexpr atom_dump::type_dump_mask DefaultAtomDumpMask = atom_dump::WithPositionMask;

typedef short _type_logs_mode;

struct Stage {
    unsigned long steps;
    double step_length;

    // dump
    bool dump_set;
    std::string dump_preset_use;
    unsigned int dump_every_steps;

    // collision
    bool collision_set;
    unsigned long collisionStep;
    int collisionLat[4];
    double pkaEnergy;
    double direction[DIMENSION];

    // thermodynamic
    bool thermo_logs_set;
    std::string thermo_logs_preset_use;
    unsigned int thermo_logs_every_steps;

    bool velocity_set;
    unsigned long velocity_step;
    long velocity_region[2 * DIMENSION];
    double velocity_value[DIMENSION];

    // rescale
    bool rescales_set;
    double rescale_t; // rescale to a temperature.
    unsigned long rescale_every; // step to do rescale

    Stage();

    void packdata(kiwi::Bundle &bundle);

    void unnpackdata(int &cursor, kiwi::Bundle &bundle);
};

struct DumpConfig {
    std::string name;
    double region[2 * DIMENSION]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    OutputMode mode;
    // total steps in this dump
    unsigned int steps;
    // path of dumped origin atoms before collision
    std::string file_path;
    // output atoms by frame if true.
    bool by_frame;
    // dump dump_mask
    atom_dump::type_dump_mask dump_mask;
    // dump the whole system, not atoms in the given region. (this field is not in config file)
    bool dump_whole_system;

    DumpConfig() : mode(OutputMode::COPY), steps(0), file_path(DEFAULT_OUTPUT_DUMP_FILE_PATH),
                   by_frame(false), dump_mask(DefaultAtomDumpMask), dump_whole_system(true) {};

    void packdata(kiwi::Bundle &bundle);

    void unnpackdata(int &cursor, kiwi::Bundle &bundle);
};

struct Output {
    // output section
    std::vector<DumpConfig> presets;

    // config fields to output thermodynamics information.
    std::vector<md_thermodynamic::OutputThermodynamic> thermo_presets;
    // logs in output section
    _type_logs_mode logs_mode;
    std::string logs_filename;

    Output() : logs_mode(LOGS_MODE_CONSOLE), logs_filename("") {}
};

struct AtomType {
    std::string name;
    double mass;
    int weight;

    void packdata(kiwi::Bundle &bundle) const;

    void unnpackdata(int &cursor, kiwi::Bundle &bundle);
};

struct ReadPhaseConfig {
    bool enable;
    unsigned int version;
    std::string file_path;
    unsigned int init_step; // initial time step.

    void packdata(kiwi::Bundle &bundle) const;

    void unpackdata(int &cursor, kiwi::Bundle &bundle);

    ReadPhaseConfig() : enable(false), version(0), init_step(0) {};
};

class ConfigValues {
    friend std::ostream &operator<<(std::ostream &os, const ConfigValues &cv);

public:
    // config values start
    // simulation section
    int64_t phaseSpace[DIMENSION];
    double cutoffRadiusFactor;
    double latticeConst;
    unsigned long timeSteps; // total steps is not set in config file, but compute from each stages.
    double timeStepLength; // default step length

    bool createPhaseMode; // enable/disable create mode
    double createTSet; // system temperature for creation.
    int createSeed;
    std::string readPhaseFilename; // for read mode

    // alloy
    int alloyCreateSeed;
    std::vector<AtomType> types;

    // read atoms from file
    ReadPhaseConfig read_phase;

    inline bool createSystemMode() const {
        return createPhaseMode;
    }

    inline bool readSystemMode() const {
        return read_phase.enable;
    }

    // potential config
    std::string potentialFileFormat;
    std::string potentialFilename;
    unsigned short potentialType;

    // simulation section ends
    // output section
    Output output;
    // config values ends
    std::vector<Stage> stages;

    ConfigValues();

    void packdata(kiwi::Bundle &bundle);  // todo override

    void unpackdata(kiwi::Bundle &bundle);

};


#endif //MISA_MD_CONFIG_VALUES_H
