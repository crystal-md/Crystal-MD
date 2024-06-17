//
// Created by genshen on 2019-08-28.
//

#include <fstream>
#include <gtest/gtest.h>
#include <utils/mpi_utils.h>
#include "config_parser.h"
#include "../test_config.h"

const std::string vsl_test_config_header = R"(
simulation:
  phasespace: [50, 50, 50]
  cutoff_radius_factor: 1.96125
  lattice_const: 2.85532
  def_timesteps_length: 0.001

potential:
  format: "setfl"
  type: "eam/alloy"
  file_path: "FeCuNi.eam.alloy"

creation:
  create_phase: true
  create_t_set: 600.0
  create_seed: 466953
  alloy:
    create_seed: 1024
    types:
      - name: Fe
        mass: 55.845
        weight: 97
      - name: Cu
        mass: 63.546
        weight: 2
      - name: Ni
        mass: 58.6934
        weight: 1
output:
  thermo:
    presets:
      - name: my_thermo
        with: [step, time, temp, pe, ke, etotal]
  logs:
    logs_mode: "console"
    logs_filename: ""
)";


// test variable step length in config parsing level.
// @MPI
TEST(config_parsing_vsl_test, config_parsing_test) {
    const std::string config_filePath = TEST_TEMP_FILES_DIR "/test_config_parsing_vsl.yaml";
    std::fstream fs(config_filePath, std::ios::out);
    if (!fs.good()) {
        GTEST_FATAL_FAILURE_("open file fail.");
    }
    fs << vsl_test_config_header << R"(
stages:
  - name: init
    step_length: 0.001
    steps: 4

  - name: collision
    step_length: 0.0001
    steps: 8

  - name: run
    step_length: 0.005
    steps: 6
)";
    fs.close();

    ConfigParser *p_config;
    if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
        // initial config Obj, then read and resolve config file.
        p_config = ConfigParser::newInstance(config_filePath); // todo config file from argv.
        if (p_config->hasError) {
            GTEST_FATAL_FAILURE_(p_config->errorMessage.c_str());
        }
    } else {
        // just initial a empty config Obj.
        p_config = ConfigParser::getInstance();
    }
    p_config->sync(); // sync config data to other processors from master processor.

    // test
    const int expected_vsl_size = 3;
    const unsigned long expected_breaks[] = {4, 8, 6};
    const double expected_vsl_length[] = {0.001, 0.0001, 0.005};

    EXPECT_EQ(p_config->configValues.stages.size(), expected_vsl_size);
    for (int i = 0; i < expected_vsl_size; i++) {
        EXPECT_EQ(p_config->configValues.stages[i].steps, expected_breaks[i]);
        EXPECT_EQ(p_config->configValues.stages[i].step_length, expected_vsl_length[i]);
    }
}
