//
// Created by genshen on 5/31/18.
//

#include <logs/logs.h>
#include "atom_dump_types.h"


size_t atom_dump::getAtomDumpDataSize(const atom_dump::type_dump_mask mask) {
    if (((mask & atom_dump::WithPositionMask) | (mask & atom_dump::WithVelocityMask) |
         (mask & atom_dump::WithForceMask)) == 0) {
        return offsetof(AtomInfoDump, inter_type) + sizeof(_type_inter_type);
    }
    size_t base = offsetof(AtomInfoDump, atom_location[0]);
    if ((mask & atom_dump::WithPositionMask) != 0) {
        base += DIMENSION * sizeof(double);
    }
    if ((mask & atom_dump::WithVelocityMask) != 0) {
        base += DIMENSION * sizeof(double);
    }
    if ((mask & atom_dump::WithForceMask) != 0) {
        base += DIMENSION * sizeof(double);
    }
    return base;
}

#define MPI_TYPE_OP_CHECK(OP) \
{                           \
  const int status = OP; \
  if (status != MPI_SUCCESS) { \
    kiwi::logs::e("mpi", "bad MPI type creation\n"); \
    MPI_Abort(MPI_COMM_WORLD, 0); \
    return MPI_DATATYPE_NULL; \
  }                           \
}

MPI_Datatype atom_dump::registerAtomDumpMPIDataType(const type_dump_mask mask) {
    if (mpi_dump_types[mask] != MPI_DATATYPE_NULL) {
        return mpi_dump_types[mask];
    }

    /* create a type for struct AtomInfoDump */
    int block_len[] = {0, 0, 0};
    if ((mask & atom_dump::WithPositionMask) != 0) {
        block_len[0] = DIMENSION;
    }
    if ((mask & atom_dump::WithVelocityMask) != 0) {
        block_len[1] = DIMENSION;
    }
    if ((mask & atom_dump::WithForceMask) != 0) {
        block_len[2] = DIMENSION;
    }
    int disp[] = {0, DIMENSION, 2 * DIMENSION};
    MPI_Datatype atom_prop = MPI_DATATYPE_NULL; // atom's position, velocity and force.
    MPI_TYPE_OP_CHECK(MPI_Type_indexed(3, block_len, disp, MPI_DOUBLE, &atom_prop))
    MPI_TYPE_OP_CHECK(MPI_Type_commit(&atom_prop))

    MPI_Datatype ext_atom_prop = MPI_DATATYPE_NULL;
    MPI_TYPE_OP_CHECK(MPI_Type_create_resized(atom_prop, 0, 3 * DIMENSION * sizeof(double), &ext_atom_prop));
    MPI_TYPE_OP_CHECK(MPI_Type_commit(&ext_atom_prop))

    int block_lengths[] = {1, 1, 1};
    MPI_Datatype types[] = {
            MPI_UNSIGNED_LONG, // id
            MPI_INT,  // type
            ext_atom_prop, // location, velocity or force
    };
    MPI_Aint offsets[] = {
            offsetof(AtomInfoDump, id),
            offsetof(AtomInfoDump, type),
            offsetof(AtomInfoDump, atom_location)
    };

    MPI_Datatype l_mpi_dump_type = MPI_DATATYPE_NULL;
    MPI_TYPE_OP_CHECK(MPI_Type_create_struct(3, block_lengths, offsets, types, &l_mpi_dump_type))
    MPI_TYPE_OP_CHECK(MPI_Type_commit(&l_mpi_dump_type))

    MPI_Datatype ext_l_mpi_dump_type = MPI_DATATYPE_NULL;
    MPI_TYPE_OP_CHECK(MPI_Type_create_resized(l_mpi_dump_type, 0, sizeof(atom_dump::AtomInfoDump), &ext_l_mpi_dump_type));

    mpi_dump_types[mask] = ext_l_mpi_dump_type;
    return ext_l_mpi_dump_type;
}

MPI_Datatype atom_dump::registerFrameMetaMPIDataType() {
    if (mpi_frame_meta_datatype != MPI_DATATYPE_NULL) {
        return mpi_frame_meta_datatype;
    }

    int block_lengths[] = {3, 1};
    MPI_Datatype types[] = {
            MPI_UNSIGNED_LONG,
            MPI_DOUBLE
    };

    MPI_Aint offsets[] = {
            offsetof(FrameMetaData, atoms_num),
            offsetof(FrameMetaData, time),
    };

    MPI_Datatype l_mpi_frame_meta_type = 0;
    MPI_TYPE_OP_CHECK(MPI_Type_create_struct(2, block_lengths, offsets, types, &l_mpi_frame_meta_type))
    MPI_TYPE_OP_CHECK(MPI_Type_commit(&l_mpi_frame_meta_type))

    mpi_frame_meta_datatype = l_mpi_frame_meta_type;
    return l_mpi_frame_meta_type;
}
