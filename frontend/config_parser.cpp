//
// Created by gensh(genshenchu@gmail.com) on 2017/4/16.
//

#include <iostream>
#include <fstream>
#include <utils/mpi_utils.h>
#include <utils/bundle.h>
#include <climits>
#include <cassert>

#include "config_parser.h"
#include "utils/rpcc.hpp"
#include "def_config_values.h"


ConfigParser *ConfigParser::m_pInstance = nullptr;

ConfigParser::ConfigParser() : kiwi::config::config() {
    configValues = ConfigValues(); // todo ConfigValues cap
}

// a simple single mode.
ConfigParser *ConfigParser::getInstance() {
    if (m_pInstance == nullptr) {
        m_pInstance = new ConfigParser();
    }
    return m_pInstance; // make sure there is a configure instance.
}

ConfigParser *ConfigParser::newInstance(const std::string &configureFilePath) {
    if (m_pInstance == nullptr) {
        m_pInstance = new ConfigParser();  // todo delete
        m_pInstance->parseConfig(configureFilePath);
    }
    return m_pInstance;
}

bool ConfigParser::configureCheck() {
    //todo
    return true;
}

// @override only for master processor.
// @see https://github.com/jbeder/yaml-cpp/wiki for more details.
void ConfigParser::parseConfig(const std::string config_file) {
    std::ifstream ifs(config_file);
    if (!ifs.good()) { // todo fixme important, if file not exist, this branch would not enter.
        setError("can not access the configure file: " + config_file);
        return;
    }
    const YAML::Node config = YAML::Load(ifs);
    ifs.close();

    // parse simulation
    if (!parseConfigSimulation(config["simulation"])) {
        return;
    }

    // parse creation
    if (!parseConfigCreation(config["creation"])) {
        return;
    }

    // parse read_phase field
    if (!parseConfigReadPhase(config["read_phase"])) {
        return;
    }

    // parse potential
    if (!parseConfigPotential(config["potential"])) {
        return;
    }

    // parse output
    if (!parseConfigOutput(config["output"])) {
        return;
    }

    // parse stages
    if (!parseStages(config["stages"])) {
        return;
    }
}

// @override
void ConfigParser::putConfigData(kiwi::Bundle &bundle) {
    configValues.packdata(bundle);
}

// @override
void ConfigParser::getConfigData(kiwi::Bundle &bundle) {
    configValues.unpackdata(bundle);
}

bool ConfigParser::parseConfigSimulation(const YAML::Node &yaml_sim) {
    if (!yaml_sim) {
        setError("\"simulation\" in config is not specified.");
        return false;
    }
    // resolve simulation.phase_space
    if (!yaml_sim["phasespace"].IsSequence() || yaml_sim["phasespace"].size() != DIMENSION) {
        setError("\"simulation.phasespace\" must be a array with length 3.");
        return false;
    } else {
        for (int i = 0; i < DIMENSION; i++) {
            configValues.phaseSpace[i] = yaml_sim["phasespace"][i].as<unsigned long>(0); // todo uint64
        }
    }

    //resolve simulation.cutoff_radius
    if (!yaml_sim["cutoff_radius_factor"]) {
        setError("cutoff_radius_factor must be set in config.");
        return false;
    } else {
        configValues.cutoffRadiusFactor = yaml_sim["cutoff_radius_factor"].as<double>(0.0);
    }

    //resolve simulation.latticeconst
    if (!yaml_sim["lattice_const"]) {
        setError("lattice const must be set in config.");
        return false;
    } else {
        configValues.latticeConst = yaml_sim["lattice_const"].as<double>(0.0);
    }

    configValues.timeStepLength = yaml_sim["def_timesteps_length"].as<double>(default_time_length);
    return true;
}

bool ConfigParser::parseConfigCreation(const YAML::Node &yaml_creation) {
    if (!yaml_creation) {
        setError("\"creation\" in config is not specified.");
        return false;
    }

    configValues.createPhaseMode = yaml_creation["create_phase"].as<bool>(default_create_phase);// default value is true
    if (configValues.createPhaseMode) { //create mode
        if (yaml_creation["create_t_set"]) {
            configValues.createTSet = yaml_creation["create_t_set"].as<double>(0.0);
        } else {
            setError("creation `create_t_set` must be specified.");
            return false;
        }
        configValues.createSeed = yaml_creation["create_seed"].as<int>(default_random_seek);
    } else {
        // read mode.
    }

    // resolve simulation.alloy
    return parseConfigAlloy(yaml_creation["alloy"]);
}

#define PARSE_FILED_ELSE_RETURN_FALSE(value, parent_node, field_name, field_message, type, fallback) \
const YAML::Node yaml_##field_name = (parent_node)[#field_name]; \
if (yaml_##field_name) { \
  (value) = yaml_##field_name.as<type>(fallback); \
} else { \
  setError(field_message " must be specified."); \
  return false; \
}


bool ConfigParser::parseConfigReadPhase(const YAML::Node &yaml_read) {
    if (!yaml_read) {
        return true; // if it is not specified, it's also ok.
    }

    PARSE_FILED_ELSE_RETURN_FALSE(configValues.read_phase.enable, yaml_read, enable, "`enable` in read_phase", bool,
                                  false);
    PARSE_FILED_ELSE_RETURN_FALSE(configValues.read_phase.version, yaml_read, version, "`enable` in read_phase",
                                  unsigned int, 0);
    PARSE_FILED_ELSE_RETURN_FALSE(configValues.read_phase.file_path, yaml_read, file_path, "`file_path` in read_phase",
                                  std::string, "");
    PARSE_FILED_ELSE_RETURN_FALSE(configValues.read_phase.init_step, yaml_read, init_step, "`init_step` in read_phase",
                                  unsigned int, 0);
    return true;
}

bool ConfigParser::parseConfigPotential(const YAML::Node &yaml_pot) {
    //potential_file
    if (!yaml_pot) {
        setError("\"potential\" in config is not specified.");
        return false;
    }
    const YAML::Node yaml_pot_format = yaml_pot["format"];
    const YAML::Node yaml_pot_file = yaml_pot["file_path"];
    const YAML::Node yaml_pot_type = yaml_pot["type"];

    if (yaml_pot_format) {
        configValues.potentialFileFormat = yaml_pot_format.as<std::string>();
    } else {
        setError("potential file format must be specified.");
        return false;
    }

    if (yaml_pot_file) {
        configValues.potentialFilename = yaml_pot_file.as<std::string>();;
    } else {
        setError("potential filename must be specified.");
        return false;
    }
    if (yaml_pot_type) {
        std::string type = yaml_pot_type.as<std::string>();;
        if (type.compare("eam/fs") == 0) {
            configValues.potentialType = EAM_TYPE_FS;
        } else if (type.compare("eam/alloy") == 0) {
            configValues.potentialType = EAM_TYPE_ALLOY;
        } else {
            setError("potential type must be specified111.");
            return false;
        }
    } else {
        setError("potential type must be specified.");
        return false;
    }
    return true;
}

bool ConfigParser::parseThermoPresets(const YAML::Node &yaml_thermo_dump) {
    if (!yaml_thermo_dump) {
        // skip the field "thermo" is ok.
        return true;
    }

    const YAML::Node yaml_thermo_preset = yaml_thermo_dump["presets"];
    if (!yaml_thermo_preset) {
        // skip "thermo.presets"
        return true;
    }

    if (!yaml_thermo_preset.IsSequence()) {
        setError("`output.thermo.presets` is not correctly set in config file.");
        return false;
    }
    const std::size_t presets_size = yaml_thermo_preset.size();

    for (std::size_t i = 0; i < presets_size; i++) {
        const YAML::Node yaml_preset = yaml_thermo_preset[i];
        md_thermodynamic::OutputThermodynamic thermodynamic_config;
        // parse a preset
        const YAML::Node yaml_with = yaml_preset["with"];
        if (!yaml_with) {
            thermodynamic_config.flags = 0x00;
        } else {
            if (!yaml_with.IsSequence()) {
                setError("\"output.thermo.presets.with\" must be a sequence.");
                return false;
            } else {
                atom_dump::type_dump_mask local_mask = 0;
                for (auto ele: yaml_with) {
                    const std::string ele_str = ele.as<std::string>();
                    if (ele_str == "time") { // current global simulation time
                        local_mask |= md_thermodynamic::WithTimeMask;
                    } else if (ele_str == "step") { // current global step
                        local_mask |= md_thermodynamic::WithStepMask;
                    } else if (ele_str == "temp") { // temperature
                        local_mask |= md_thermodynamic::WithTemperatureMask;
                    } else if (ele_str == "pe") { // potential energy
                        local_mask |= md_thermodynamic::WithPotentialEnergyMask;
                    } else if (ele_str == "ke") { // kinetic energy
                        local_mask |= md_thermodynamic::WithKineticEnergyMask;
                    } else if (ele_str == "etotal") { // total energy
                        local_mask |= md_thermodynamic::WithKineticEnergyMask;
                    } else {
                        setError("unrecognized value `" + ele_str + "` for \"output.thermo.presets.with\".");
                        return false;
                    }
                }
                thermodynamic_config.flags = local_mask;
            }
        }

        thermodynamic_config.name = yaml_preset["name"].as<std::string>("default");
        configValues.output.thermo_presets.emplace_back(thermodynamic_config);
    }
    return true;
}

bool ConfigParser::parseDumpPresets(const YAML::Node &yaml_atom_dump) {
    if (!yaml_atom_dump) {
        // skip "atom_dump"
        return true;
    }
    const YAML::Node yaml_dump_preset = yaml_atom_dump["presets"];
    if (!yaml_dump_preset) {
        // skip "atom_dump.presets"
        return true;
    }

    if (!yaml_dump_preset.IsSequence()) {
        setError("`output.atom_dump.presets` is not correctly set in config file.");
        return false;
    }
    const std::size_t presets_size = yaml_dump_preset.size();

    for (std::size_t i = 0; i < presets_size; i++) {
        const YAML::Node yaml_preset = yaml_dump_preset[i];
        DumpConfig dump_config;
        // dump mode
        if (yaml_preset["mode"]) {
            if ("copy" == yaml_preset["mode"].as<std::string>("copy")) { // todo equal?
                dump_config.mode = OutputMode::COPY;
            } else {
                dump_config.mode = OutputMode::DEBUG;
            }
        }
        // dump region
        const YAML::Node yaml_region = yaml_preset["region"];
        if (yaml_region) {
            dump_config.dump_whole_system = false;
            if (yaml_region.IsSequence() && yaml_region.size() == 2 * DIMENSION) {
                for (std::size_t r = 0; r < 2 * DIMENSION; r++) {
                    dump_config.region[r] = yaml_region[r].as<double>(0.0);
                }
            } else { //the array length must be 6.
                setError("array length of value \"output.atom_dump.presets.region\" must be 6.");
                return false;
            }
        } else {
            dump_config.dump_whole_system = true;
        }

        // parse `with` field
        const YAML::Node yaml_with = yaml_preset["with"];
        if (!yaml_with) {
            dump_config.dump_mask = DefaultAtomDumpMask;
        } else {
            if (yaml_with.IsSequence()) {
                atom_dump::type_dump_mask local_mask = 0;
                for (auto ele : yaml_with) {
                    const std::string ele_str = ele.as<std::string>();
                    if (ele_str == "location") {
                        local_mask |= atom_dump::WithPositionMask;
                    } else if (ele_str == "velocity") {
                        local_mask |= atom_dump::WithVelocityMask;
                    } else if (ele_str == "force") {
                        local_mask |= atom_dump::WithForceMask;
                    } else {
                        setError("unrecognized value `" + ele_str + "` for \"output.atom_dump.presets.with\".");
                        return false;
                    }
                }
                dump_config.dump_mask = local_mask;
            } else {
                setError("\"output.atom_dump.presets.with\" must be a sequence.");
                return false;
            }
        }

        dump_config.by_frame = yaml_preset["by_frame"].as<bool>(false);
        dump_config.name = yaml_preset["name"].as<std::string>("default");
        dump_config.file_path = yaml_preset["file_path"].as<std::string>(DEFAULT_OUTPUT_DUMP_FILE_PATH);
        // todo check if it is a real path.
        if (dump_config.by_frame && dump_config.file_path.find("{}") == std::string::npos) {
            setError("error format of dump file path: " + dump_config.file_path);
        }
        configValues.output.presets.emplace_back(dump_config);
    }
    return true;
}

bool ConfigParser::parseConfigOutput(const YAML::Node &yaml_output) {
    if (!yaml_output) {
        setError("output section in config file is not specified.");
        return false;
    }
    // resole dump config.
    const YAML::Node yaml_atom_dump = yaml_output["atom_dump"];
    if (!parseDumpPresets(yaml_atom_dump)) {
        return false;
    }

    // resole thermodynamics config
    const YAML::Node yaml_thermo_dump = yaml_output["thermo"];
    if (yaml_thermo_dump["interval"]) {
        // error only for compatibility issue in v0.4.0*.
        setError("`output.thermo.interval` is deprecated now, please use `output.thermo.presets` instead.");
        return false;
    }
    if (!parseThermoPresets(yaml_thermo_dump)) {
        return false;
    }

    // resolve logs.
    const YAML::Node yaml_logs = yaml_output["logs"];
    if (!yaml_logs) {
        setError("output logs in config file is not specified.");
        return false;
    }

    std::string logs_mode_str = yaml_logs["logs_mode"].as<std::string>(DEFAULT_LOGS_MODE_CONSOLE_STRING);
    if (LOGS_MODE_CONSOLE_STRING == logs_mode_str) {
        configValues.output.logs_mode = LOGS_MODE_CONSOLE;
    } else {
        configValues.output.logs_mode = LOGS_MODE_FILE;
        configValues.output.logs_filename = yaml_logs["logs_filename"].as<std::string>("");
        // if filename is empty, we will generate a filename.
        if (configValues.output.logs_filename.empty()) {
            std::ostringstream str_stream;
            str_stream << "md.";
            str_stream << rpcc();
            str_stream << ".log";
            configValues.output.logs_filename = str_stream.str();
        }
    }
    return true;
}

bool ConfigParser::parseConfigAlloy(const YAML::Node &yaml_alloy) {
    if (!yaml_alloy) {
        setError("alloy config is not specified in creation mode.");
        return false;
    }
    configValues.alloyCreateSeed = yaml_alloy["create_seed"].as<int>(default_random_seek);

    const YAML::Node yaml_types = yaml_alloy["types"];
    if (!yaml_types || !yaml_types.IsSequence()) {
        setError("alloy types must be an array.");
        return false;
    } else {
        for (const auto &yaml_type : yaml_types) {
            AtomType atom_type;
            atom_type.name = yaml_type["name"].as<std::string>("undefined");
            atom_type.mass = yaml_type["mass"].as<double>(1.0);
            atom_type.weight = yaml_type["weight"].as<int>(1);
            configValues.types.emplace_back(atom_type);
        }
    }
    return true;
}

bool ConfigParser::resolveConfigCollision(Stage *stage, const YAML::Node &yaml_collision) {
    if (!yaml_collision) {
        setError("collision config is not set.");
        return false;
    }
    stage->collisionStep = yaml_collision["collision_step"].as<unsigned long>(0); // todo default 0.

    const YAML::Node yaml_lat = yaml_collision["lat"];
    if (yaml_lat.IsSequence() && yaml_lat.size() == 4) {
        for (std::size_t i = 0; i < 4; i++) {
            stage->collisionLat[i] = yaml_lat[i].as<int64_t>();  // todo conversion.
        }
    } else { //the array length must be 4.
        setError("array length of value \"collision.lat\" must be 4.");
        return false;
    }
    stage->pkaEnergy = yaml_collision["energy"].as<double>(0.0);


    const YAML::Node yaml_dir = yaml_collision["direction"];
    if (yaml_dir.IsSequence() && yaml_dir.size() == DIMENSION) {
        for (std::size_t i = 0; i < DIMENSION; i++) {
            stage->direction[i] = yaml_dir[i].as<double>(0.0);  // todo conversion.
        }
    } else { //the array length must be 3.
        setError("array length of value \"simulation.collision.collision_v\" must be 3.");
        return false;
    }
    return true;
}

bool ConfigParser::resolveStageDump(Stage *stage, const YAML::Node &yaml_stage_dump) {

    stage->dump_preset_use = yaml_stage_dump["use"].as<std::string>("default");
    // todo 0 as default.
    stage->dump_every_steps = yaml_stage_dump["every_steps"].as<unsigned long>(ULONG_MAX);

    // verify (e.g. the preset must exist) and calculate total step in this dump preset.
    bool found = false;
    size_t index = 0;
    for (auto &p : configValues.output.presets) {
        if (p.name == stage->dump_preset_use) {
            found = true;
            break;
        }
        index++;
    }
    if (!found) {
        setError("dump preset not found for `" + stage->dump_preset_use + "`");
        return false;
    }
    configValues.output.presets[index].steps += stage->steps / stage->dump_every_steps;
    return true;
}

bool ConfigParser::resolveStageThermoLogs(Stage *stage, const YAML::Node &yaml_thermo_logs) {
    stage->thermo_logs_preset_use = yaml_thermo_logs["use"].as<std::string>("default");
    stage->thermo_logs_every_steps = yaml_thermo_logs["every_steps"].as<unsigned int>(0);

    // verify if the specified preset exists.
    bool found = false;
    std::size_t index = 0;
    for (auto &p: configValues.output.thermo_presets) {
        if (p.name == stage->thermo_logs_preset_use) {
            found = true;
            break;
        }
        index++;
    }
    if (!found) {
        setError("the thermo_logs preset not found for `" + stage->thermo_logs_preset_use + "`");
        return false;
    }
    configValues.output.thermo_presets[index].steps += stage->steps / stage->thermo_logs_every_steps;
    return true;
}

bool ConfigParser::resolveConfigVelocity(Stage *stage, const YAML::Node &yaml_velocity) {
    if (!yaml_velocity) {
        setError("`velocity` config is not set.");
        return false;
    }
    stage->velocity_step = yaml_velocity["step"].as<unsigned long>(0);

    const YAML::Node yaml_velocity_value = yaml_velocity["v"];
    if (yaml_velocity_value.IsSequence() && yaml_velocity_value.size() == 3) {
        for (std::size_t i = 0; i < 3; i++) {
            stage->velocity_value[i] = yaml_velocity_value[i].as<double>();
        }
    } else { //the array length must be 3.
        setError("array length of value \"velocity.v\" must be 3.");
        return false;
    }

    const YAML::Node yaml_region = yaml_velocity["region"];
    if (yaml_region.IsSequence() && yaml_region.size() == 6) {
        for (std::size_t i = 0; i < 6; i++) {
            stage->velocity_region[i] = yaml_region[i].as<long>();
        }
    } else { //the array length must be 6.
        setError("array length of value \"velocity.region\" must be 6.");
        return false;
    }
    return true;
}

bool ConfigParser::resolveConfigRescale(Stage *stage, const YAML::Node &yaml_rescale) {
    if (!yaml_rescale) {
        setError("rescale config is not set.");
        return false;
    }
    stage->rescale_t = yaml_rescale["t"].as<double>(0.0);
    stage->rescale_every = yaml_rescale["every_steps"].as<unsigned long>(0);
    return true;
}

bool ConfigParser::parseStages(const YAML::Node &yaml_stages) {
    if (!yaml_stages && yaml_stages.IsSequence()) {
        setError("stages is not correctly set in config file.");
        return false;
    }

    //resolve simulation.collision
//    resolveConfigCollision(table->get_table("collision"));
    std::vector<unsigned long> vsl_break_points;
    std::vector<double> vsl_lengths;
    const std::size_t stages_size = yaml_stages.size();
    for (std::size_t i = 0; i < stages_size; i++) {
        const YAML::Node yaml_stage = yaml_stages[i];
        Stage stage;
        stage.steps = yaml_stage["steps"].as<unsigned long>(0);
        stage.step_length = yaml_stage["step_length"].as<double>(0.0);
        // collision
        if (yaml_stage["set_v"]) {
            stage.collision_set = true;
            if (!resolveConfigCollision(&stage, yaml_stage["set_v"])) {
                return false;
            }
        } else {
            stage.collision_set = false;
        }
        // dump
        if (yaml_stage["dump"]) {
            stage.dump_set = true;
            if (!resolveStageDump(&stage, yaml_stage["dump"])) {
                return false;
            }
        } else {
            stage.dump_set = false;
        }
        // thermo_logs
        if (yaml_stage["thermo_logs"]) {
            stage.thermo_logs_set = true;
            if (!resolveStageThermoLogs(&stage, yaml_stage["thermo_logs"])) {
                return false;
            }
        } else {
            stage.thermo_logs_set = false;
        }
        // velocity
        if (yaml_stage["velocity"]) {
            stage.velocity_set = true;
            if (!resolveConfigVelocity(&stage, yaml_stage["velocity"])) {
                return false;
            }
        } else {
            stage.velocity_set = false;
        }
        // rescale
        if (yaml_stage["rescale"]) {
            stage.rescales_set = true;
            if (!resolveConfigRescale(&stage, yaml_stage["rescale"])) {
                return false;
            }
        } else {
            stage.rescales_set = false;
        }
        configValues.stages.emplace_back(stage);
        configValues.timeSteps += stage.steps; // calculate total steps.
    }
    return true;
}
