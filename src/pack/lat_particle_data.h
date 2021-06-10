//
// Created by baihe back to 2016-06-28.
//

#ifndef MISA_MD_LAT_PARTICLE_DATA_H
#define MISA_MD_LAT_PARTICLE_DATA_H

#include <mpi.h>
#include "../types/atom_types.h"
#include "md_building_config.h"

class LatParticleData {
public:
    // 定义一个数据类型，用于MPI数据传输
    static void setMPIType(MPI_Datatype &sendPartType);

    LatParticleData();

    atom_type::atom_type type;
#ifdef DEV_MD_COMM_INC_ATOM_ID
    _type_atom_id id;
#endif
    double r[3];
};

#endif //MISA_MD_LAT_PARTICLE_DATA_H
