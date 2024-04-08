//
// Created by genshen on 2023/12/4.
//

#ifndef MISA_MD_CONFIG_VALUES_THERMODYNAMIC_HPP
#define MISA_MD_CONFIG_VALUES_THERMODYNAMIC_HPP

#include <string>
#include <fstream>

#include <utils/bundle.h>

// the namespace for outputting of system thermodynamic, such as energy, temperature.
namespace md_thermodynamic {
    typedef unsigned long FlagOutThermodynamic;

    constexpr FlagOutThermodynamic WithTimeMask = 0x1 << 0;
    constexpr FlagOutThermodynamic WithStepMask = 0x1 << 1;
    constexpr FlagOutThermodynamic WithTemperatureMask = 0x1 << 2;
    constexpr FlagOutThermodynamic WithPotentialEnergyMask = 0x1 << 3;
    constexpr FlagOutThermodynamic WithKineticEnergyMask = 0x1 << 4;

    /**
     * This class is used as thermodynamic output preset.
     * see `output.thermo.presets` in the example input config file.
     */
    struct OutputThermodynamic {
        /**
         * a flag to record which values should be printed on screen or log file.
         */
        FlagOutThermodynamic flags;
        /**
         * the name of current preset config.
         */
        std::string name;

        /**
         * total steps in this preset. 
         * It can be calculated on the fly and does not need sync among processor.
         */
        unsigned int steps = 0;
    public:
        /**
         * pack data for MPI communication.
         * @param bundle data buffer for storing the data
         */
        void packdata(kiwi::Bundle &bundle) {
            bundle.put(flags);
            bundle.put(name);
        }

        /**
         * unpack data from a buffer of MPI communication.
         * @param cursor cursor in the input data buffer
         * @param bundle input data buffer
         */
        void unnpackdata(int &cursor, kiwi::Bundle &bundle) {
            bundle.get(cursor, flags);
            bundle.get(cursor, name);
        }
    };

    inline std::ostream &operator<<(std::ostream &os, const OutputThermodynamic &value) {
        os << value.name << "," << value.flags;
        return os;
    }
}

#endif //MISA_MD_CONFIG_VALUES_THERMODYNAMIC_HPP
