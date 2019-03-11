//
// Created by genshen on 5/31/18.
//

#ifndef CRYSTAL_MD_ATOM_INFO_DUMP_H
#define CRYSTAL_MD_ATOM_INFO_DUMP_H


#include <cstddef>
#include <mpi.h>
#include "pre_define.h"

namespace atom_dump {
    struct AtomInfoDump {
    public:
        _type_atom_id id;
        std::size_t step;
        _type_atom_type type;
        _type_inter_type inter_type;
        _type_atom_location atom_location[DIMENSION]; // atom location
        _type_atom_velocity atom_velocity[DIMENSION]; // atom velocity
    };

    static MPI_Datatype mpi_dump_type;

    void registerAtomDumpMPIDataType();
}

#endif //CRYSTAL_MD_ATOM_INFO_DUMP_H
