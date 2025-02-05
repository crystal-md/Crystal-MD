//
// Created by genshen on 2019-04-25.
//

#include "system_configuration.h"
#include "atom/atom_element.h"

std::array<_type_atom_force, DIMENSION> configuration::systemForce(
        AtomList *atom_list, InterAtomList *inter_atom_list) {
    _type_atom_force force_x = 0.0, force_y = 0.0, force_z = 0.0;
    atom_list->foreachSubBoxAtom([atom_list, &force_x, &force_y, &force_z](const _type_atom_index gid) {
        MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
        if (MD_GET_ATOM_TYPE(_atom_ref, gid) != atom_type::INVALID) {
            force_x += MD_GET_ATOM_F(_atom_ref, gid, 0);
            force_y += MD_GET_ATOM_F(_atom_ref, gid, 1);
            force_z += MD_GET_ATOM_F(_atom_ref, gid, 2);
        }
    });
    for (AtomElement &atom_ele:inter_atom_list->inter_list) {
        force_x += atom_ele.f[0];
        force_y += atom_ele.f[1];
        force_z += atom_ele.f[2];
    }
    return std::array<_type_atom_force, DIMENSION>{force_x, force_y, force_z};
}

double configuration::kineticEnergy(AtomList *atom_list, InterAtomList *inter_atom_list,
                                    ReturnMod mode, const kiwi::RID root) {
    const double energy = mvv(atom_list, inter_atom_list);
    // sum energy from all processes.
    return reduceEnergy(energy, mode, root);
}

double configuration::kineticEnergy(AtomList *atom_list, ReturnMod mode, const kiwi::RID root) {
    const double energy = mvv(atom_list, nullptr);
    // sum energy from all processes.
    return reduceEnergy(energy, mode, root);
}

double configuration::kineticEnergy(InterAtomList *inter_atom_list, ReturnMod mode, const kiwi::RID root) {
    const double energy = mvv(nullptr, inter_atom_list);
    // sum energy from all processes.
    return reduceEnergy(energy, mode, root);
}

double configuration::temperature(const double ke, const _type_lattice_size n) {
    // ke = 3nkT/2
    const double to_T = mvv2e / ((3 * n - 3) * BOLTZ); // 2 * ke * mvv2e / ((3 * n - 3) * BOLTZ); /
    return 2 * ke * to_T; // todo better times order for precision.
}

double configuration::temperature(const _type_atom_count n_atoms,
                                  AtomList *atom_list, InterAtomList *inter_atom_list) {
    const double _mvv = mvv(atom_list, inter_atom_list);
    double t_global;
    MPI_Allreduce(&_mvv, &t_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // The factor 3(n-1) appears because the center of mass (COM) is fixed in space.
    const _type_atom_count dof = 3 * n_atoms - 3;
    return t_global * mvv2e / (dof * BOLTZ);
}

double configuration::mvv(AtomList *atom_list, InterAtomList *inter_atom_list) {
    double energy = 0.0;
    if (atom_list) {
        atom_list->foreachSubBoxAtom([atom_list, &energy](const _type_atom_index gid) {
            MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
            if (MD_GET_ATOM_TYPE(_atom_ref, gid) != atom_type::INVALID) {
                energy += (MD_GET_ATOM_V(_atom_ref, gid, 0) * MD_GET_ATOM_V(_atom_ref, gid, 0) +
                           MD_GET_ATOM_V(_atom_ref, gid, 1) * MD_GET_ATOM_V(_atom_ref, gid, 1) +
                           MD_GET_ATOM_V(_atom_ref, gid, 2) * MD_GET_ATOM_V(_atom_ref, gid, 2)) *
                          atom_type::getAtomMass(MD_GET_ATOM_TYPE(_atom_ref, gid));
            }
        });
    }
    if (inter_atom_list) {
        for (_type_inter_list::iterator itl = inter_atom_list->inter_list.begin();
             itl != inter_atom_list->inter_list.end(); ++itl) {
            AtomElement &_atom_ref = *itl;
            energy += (_atom_ref.v[0] * _atom_ref.v[0] +
                       _atom_ref.v[1] * _atom_ref.v[1] +
                       _atom_ref.v[2] * _atom_ref.v[2]) *
                      atom_type::getAtomMass(_atom_ref.type);
        }
    }
    return energy;
}

void configuration::rescale(const double T, const _type_atom_count n_atoms_global,
                            AtomList *atom_list, InterAtomList *inter_atom_list) {
    const double scalar = temperature(n_atoms_global, atom_list, inter_atom_list);

    /**
     * \sum { m_i(v_i)^2 }= 3*nkT  => scale = T = \sum { m_i(v_i)^2 / 3nk }
     * thus: T / T_set =  \sum { m_i(v_i)^2 } / \sum { m_i(v'_i)^2 }
     * then: \sum { m_i(v'_i)^2 } = \sum{ m_i(v_i)^2 }* (T_set / T) = \sum{ m_i(v_i * rescale_factor)^2 }
     * so, v'_i = v_i * rescale_factor
     */
    const double rescale_factor = sqrt(T / scalar);

    // perform resale
    atom_list->foreachSubBoxAtom([atom_list, rescale_factor](const _type_atom_index gid) {
        MD_LOAD_ATOM_VAR(_atom_ref, atom_list, gid);
        MD_SET_ATOM_V(_atom_ref, gid, 0, rescale_factor * MD_GET_ATOM_V(_atom_ref, gid, 0));
        MD_SET_ATOM_V(_atom_ref, gid, 1, rescale_factor * MD_GET_ATOM_V(_atom_ref, gid, 1));
        MD_SET_ATOM_V(_atom_ref, gid, 2, rescale_factor * MD_GET_ATOM_V(_atom_ref, gid, 2));
    });
    for (_type_inter_list::iterator itl = inter_atom_list->inter_list.begin();
         itl != inter_atom_list->inter_list.end(); ++itl) {
        AtomElement &_atom_ref = *itl;
        _atom_ref.v[0] *= rescale_factor;
        _atom_ref.v[1] *= rescale_factor;
        _atom_ref.v[2] *= rescale_factor;
    }
}

double configuration::reduceEnergy(const double mvv, const ReturnMod mode, const kiwi::RID root) {
    double e_global = 0;
    switch (mode) {
        case Local:
            return mvv / 2;
        case Root:
            MPI_Reduce(&mvv, &e_global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            break;
        case All:
            MPI_Allreduce(&mvv, &e_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            break;
    }
    return e_global / 2;
}
