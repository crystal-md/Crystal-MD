//
// Created by genshen(genshenchu@gmail.com) on 2017/4/16.
//

#include <yaml-cpp/yaml.h>
#include "config_values.h"
#include "config/config.h"

#ifndef MISA_MD_CONFIG_PARSER_H
#define MISA_MD_CONFIG_PARSER_H


const unsigned short EAM_TYPE_ALLOY = 0; // keep the same as `EAM_STYLE_ALLOY` in "eam.h"
const unsigned short EAM_TYPE_FS = 1; // keep the same as `EAM_STYLE_FS` in "eam.h"

/**
 * in our design, rank MASTER read config file,
 * and sync the config data to other processors.
 */
class ConfigParser : public kiwi::config {
public:
    ConfigParser();

    // config values will be stored here.
    ConfigValues configValues;

    /**
     * Create a new config instance first,
     * then it will try to resolve the toml config file.
     * and config data from resolving will be stored in {@var configValues}
     * @param configureFilePath path of toml config file.
     * @return an pointer to {@class ConfigParser}.
     */
    static ConfigParser *newInstance(const std::string &configureFilePath);

    /**
     * Just create an empty {@class ConfigParser} class (not resolve toml config file).
     * @return
     */
    static ConfigParser *getInstance();

    bool configureCheck();

private:

    // pointer of this {@class ConfigParser} for single mode.
    static ConfigParser *m_pInstance; // stored in static area.

    /**
     * resolve toml format config using lib:  https://github.com/jbeder/yaml-cpp
     * visit the url above to get more information.
     * @param config_file the path of config file.
     */
    void parseConfig(const std::string config_file);

    // resolve "simulation" section in yaml config file.
    bool parseConfigSimulation(const YAML::Node &yaml_sim);

    // resolve "create" section in yaml config file.
    bool parseConfigCreation(const YAML::Node &yaml_creation);

    // resolve "read_phase" section in in yaml config file.
    bool parseConfigReadPhase(const YAML::Node &yaml_read);

    // resolve "potential" section in yaml config file.
    bool parseConfigPotential(const YAML::Node &yaml_pot);

    // resolve "creation.alloy" section in yaml config file.
    bool parseConfigAlloy(const YAML::Node &yaml_alloy);

    // resolve "collision" in stages section in yaml config file.
    bool resolveConfigCollision(Stage *stage, const YAML::Node &yaml_collision);

    // resolve "dump" in stages section in yaml config file.
    bool resolveStageDump(Stage *stage, const YAML::Node &yaml_stage_dump);

    // resolve "thermo_logs" field in stage section in yaml config file.
    bool resolveStageThermoLogs(Stage *stage, const YAML::Node &yaml_thermo_logs);

    // resolve "velocity" in stages section in yaml config file.
    bool resolveConfigVelocity(Stage *stage, const YAML::Node &yaml_velocity);

    // resolve "output" section in yaml config file.
    bool parseConfigOutput(const YAML::Node &yaml_output);

    // resolve "output.atom_dump.presets" section in yaml config file.
    bool parseDumpPresets(const YAML::Node &yaml_atom_dump);

    // resolve "output.thermo.presets" section in yaml config file.
    bool parseThermoPresets(const YAML::Node &yaml_atom_dump);

    // resolve "rescale" in stages section in yaml config file.
    bool resolveConfigRescale(Stage *stage, const YAML::Node &yaml_rescale);

    // resolve "stages" section in yaml config file.
    bool parseStages(const YAML::Node &yaml_stages);

    /**
     * [master] put data into bundle, in which bundle is used to buffer config data.
     * later, the bundle will be sent to the other processors (sync config).
     * @param bundle used for buffering config data.
     */
    void putConfigData(kiwi::Bundle &bundle) override;

    /**
     * [non-master] get config data from bundle after finishing sync.
     * @param bundle
     */
    void getConfigData(kiwi::Bundle &bundle) override;

};


#endif //MISA_MD_CONFIG_PARSER_H
