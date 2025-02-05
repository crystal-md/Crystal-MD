//
// Created by genshen on 5/5/18.
//

#include <cmath>
#include <iostream>
#include <logs/logs.h>
#include <random>

#include "utils/random/random.h"
#include "world_builder.h"
#include "system_configuration.h"

WorldBuilder::WorldBuilder() : box_x(0), box_y(0), box_z(0), tset(0), _atoms_ratio() {
    _p_domain = nullptr;
    _p_atom = nullptr;
}

WorldBuilder &WorldBuilder::setDomain(comm::BccDomain *p_domain) {
    this->_p_domain = p_domain;
    return *this;
}

WorldBuilder &WorldBuilder::setAtomsContainer(AtomSet *p_atom) {
    this->_p_atom = p_atom;
    return *this;
}

WorldBuilder &WorldBuilder::setRandomSeed(int seed) {
    if (seed == md_rand::seed_auto) {
        std::random_device rd;
        uint32_t auto_seed = rd();
        kiwi::logs::v("random", "random seed was set to {} in auto mode.\n", auto_seed);
        md_rand::initSeed(auto_seed);
    } else {
        md_rand::initSeed(seed);
    }
    return *this;
}

WorldBuilder &WorldBuilder::setTset(double t_set) {
    this->tset = t_set;
    return *this;
}

WorldBuilder &WorldBuilder::setLatticeConst(double lattice_const) {
    this->_lattice_const = lattice_const;
    return *this;
}

WorldBuilder &WorldBuilder::setAlloyRatio(const std::vector<tp_atom_type_weight> &types_weight) {
    _atoms_ratio = types_weight;
    return *this;
}

WorldBuilder &WorldBuilder::setBoxSize(int64_t _box_x, int64_t _box_y, int64_t _box_z) {
    this->box_x = _box_x;
    this->box_y = _box_y;
    this->box_z = _box_z;
    return *this;
}

void WorldBuilder::build() {
    if (_p_atom == nullptr) {
        throw std::invalid_argument("no atom container");
        // todo return error.
    }
    if (_p_domain == nullptr) {
        throw std::invalid_argument("no domain");
        // todo return error
    }

    createPhaseSpace();

    double p[4] = {0.0, 0.0, 0.0, 0.0}; // index 0-2: mv in 3 d; index 3: mass total at local.
    double _vcm[4]; // same as @var p, but global.
    vcm(p);
    MPI_Allreduce(p, _vcm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const _type_atom_count n_atoms_global =
            2 * (unsigned long) box_x * (unsigned long) box_y * (unsigned long) box_z;  // todo type
    double &mass_total = _vcm[3];

    if (mass_total > 0.0) {
        _vcm[0] /= n_atoms_global; // the momentum to be cut off.
        _vcm[1] /= n_atoms_global;
        _vcm[2] /= n_atoms_global;
    }

    zeroMomentum(_vcm);
#ifdef MD_DEV_MODE
    vcm(p);
    kiwi::logs::d("momentum", "momentum:{0} {1} {2}\n", p[0], p[1], p[2]);
#endif
    // set temperature for all atoms
    if (tset != 0.0) {
        const _type_atom_count n_global_atoms = 2 * _p_domain->phase_space[0] *
                                                _p_domain->phase_space[1] *
                                                _p_domain->phase_space[2];
        configuration::rescale(tset, n_global_atoms, _p_atom->atom_list, _p_atom->inter_atom_list);
    }
}

void WorldBuilder::createPhaseSpace() {
    unsigned long id_pre = (unsigned long) box_x * box_y * _p_domain->dbx_sub_box_lattice_region.z_low * 2
                           + (unsigned long) _p_domain->dbx_sub_box_lattice_region.y_low *
                             box_x * _p_domain->dbx_sub_box_lattice_size[2] * 2
                           + (unsigned long) _p_domain->dbx_sub_box_lattice_region.x_low *
                             _p_domain->dbx_sub_box_lattice_size[1] *
                             _p_domain->dbx_sub_box_lattice_size[2];
    _type_atom_mass mass = 0;
    md_rand::type_rng atom_type_rng;
    atom_type_rng.seed(0);
    for (int k = 0; k < _p_domain->dbx_sub_box_lattice_size[2]; k++) {
        for (int j = 0; j < _p_domain->dbx_sub_box_lattice_size[1]; j++) {
            for (int i = 0; i < _p_domain->dbx_sub_box_lattice_size[0]; i++) {
                const _type_atom_index gid = _p_atom->getAtomList()->_atoms.getAtomIndexInSubBox(i, j, k);
                MD_LOAD_ATOM_VAR(atom_, _p_atom->getAtomList(), gid);
                MD_SET_ATOM_ID(atom_, gid, ++id_pre);
                MD_SET_ATOM_TYPE(atom_, gid, randomAtomsType(atom_type_rng)); // set random atom type.
                mass = atom_type::getAtomMass(MD_GET_ATOM_TYPE(atom_, gid)); // get atom mass of this kind of atom.
                MD_SET_ATOM_X(atom_, gid, 0,
                              (_p_domain->dbx_sub_box_lattice_region.x_low + i) * 0.5 * (_lattice_const));
                MD_SET_ATOM_X(atom_, gid, 1, (_p_domain->dbx_sub_box_lattice_region.y_low + j) * _lattice_const +
                                             (i % 2) * (_lattice_const / 2));
                MD_SET_ATOM_X(atom_, gid, 2, (_p_domain->dbx_sub_box_lattice_region.z_low + k) * _lattice_const +
                                             (i % 2) * (_lattice_const / 2));
                MD_SET_ATOM_V(atom_, gid, 0, (md_rand::random() - 0.5) / mass);
                MD_SET_ATOM_V(atom_, gid, 1, (md_rand::random() - 0.5) / mass);
                MD_SET_ATOM_V(atom_, gid, 2, (md_rand::random() - 0.5) / mass);
            }
        }
    }
}

/**
 * To achieve the goal of zero momentum,
 * we denote the total momentum of the system as vcm, the total mass of all atoms in system as M.
 * And we assume the total count of atoms is N.
 *
 * We cutoff the momentum of each atom by vcm/N, that is $ v_i' = v_i - vcm/(N * m_i) $.
 * in which, v_i is the velocity before cutting of, v_i' is the velocity after cutting of.
 * Because, \sum_{i=1}^{N} v_i * m_i = \sum_{i=1}^{N} vcm/N = vcm.
 * Thus, \sum_{i=1}^{N} v_i' * m_i  =0.
 */
void WorldBuilder::zeroMomentum(double *vcm) {
    _type_atom_mass mass;
    for (int k = 0; k < _p_domain->dbx_sub_box_lattice_size[2]; k++) {
        for (int j = 0; j < _p_domain->dbx_sub_box_lattice_size[1]; j++) {
            for (int i = 0; i < _p_domain->dbx_sub_box_lattice_size[0]; i++) {
                const _type_atom_index gid = _p_atom->getAtomList()->_atoms.getAtomIndexInSubBox(i, j, k);
                MD_LOAD_ATOM_VAR(atom_, _p_atom->getAtomList(), gid);
                mass = atom_type::getAtomMass(MD_GET_ATOM_TYPE(atom_, gid));
                MD_ADD_ATOM_V(atom_, gid, 0, -(vcm[0] / mass));
                MD_ADD_ATOM_V(atom_, gid, 1, -(vcm[1] / mass));
                MD_ADD_ATOM_V(atom_, gid, 2, -(vcm[2] / mass));
            }
        }
    }
}

void WorldBuilder::vcm(double p[DIMENSION + 1]) {
    _type_atom_mass mass_one = 0.0;
    // reset p.
    for (int i = 0; i < DIMENSION; i++) {
        p[i] = 0;
    }
    p[DIMENSION] = 0;

    for (int k = 0; k < _p_domain->dbx_sub_box_lattice_size[2]; k++) {
        for (int j = 0; j < _p_domain->dbx_sub_box_lattice_size[1]; j++) {
            for (int i = 0; i < _p_domain->dbx_sub_box_lattice_size[0]; i++) {
                const _type_atom_index gid = _p_atom->getAtomList()->_atoms.getAtomIndexInSubBox(i, j, k);
                MD_LOAD_ATOM_VAR(atom_, _p_atom->getAtomList(), gid);
                mass_one = atom_type::getAtomMass(MD_GET_ATOM_TYPE(atom_, gid));
                p[0] += MD_GET_ATOM_V(atom_, gid, 0) * mass_one;
                p[1] += MD_GET_ATOM_V(atom_, gid, 1) * mass_one;
                p[2] += MD_GET_ATOM_V(atom_, gid, 2) * mass_one;
                p[3] += mass_one; // all mass.
            }
        }
    }
}

atom_type::atom_type WorldBuilder::randomAtomsType(md_rand::type_rng &rng) {
    unsigned int ratio_total = 0;
    for (int i = 0; i < atom_type::num_atom_types; i++) {
        ratio_total += _atoms_ratio[i];
    }
#ifdef MD_DEV_MODE
    unsigned int rand_ = rand() % ratio_total;
#else
    unsigned int rand_ = rng() % ratio_total; // todo srank. Rand() has limited randomness; use C++ lib instead.
#endif
//    return atom_type::getAtomTypeByOrder();
    unsigned int rank_local = 0;
    for (int i = 0; i < atom_type::num_atom_types; i++) {
        rank_local += _atoms_ratio[i];
        if (rand_ < rank_local) {
            return atom_type::getAtomTypeByNum(i);
        }
    }
//    return atom_type::Fe;
}
