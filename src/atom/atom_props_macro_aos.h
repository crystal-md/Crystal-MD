//
// Created by chu genshen on 2021/8/24.
//

#ifndef MISA_MD_ATOM_PROPS_MACRO_AOS_H
#define MISA_MD_ATOM_PROPS_MACRO_AOS_H

#include "types/atom_types.h"
#include "arch/arch_atom_list_collection.h"

#define MD_HASH_LIST_DECLARE() \
AtomPropList<AtomElement> _atoms; // atoms in 3d.

#define MD_HASH_LIST_INIT(lattice) \
_atoms(AtomPropList<AtomElement>(lattice))

#define MD_HASH_LIST_DESTROY(lattice) \
_atoms.destroyPropList();

/**
 * id: global id
 */
#define MD_LOAD_ATOM_VAR(name, list, id)  \
AtomElement &name = list->_atoms.getAtomEleByLinearIndex(id)

#define MD_GET_ATOM_ID(name, gid) \
name.id

#define MD_SET_ATOM_ID(name, gid, _id) \
name.id = _id

#define MD_GET_ATOM_TYPE(name, gid) \
name.type

#define MD_SET_ATOM_TYPE(name, gid, _tp) \
name.type = _tp

#define MD_IS_ATOM_TYPE_INTER(name, gid) \
(name.type == atom_type::INVALID)

#define MD_GET_ATOM_X(name, gid, i) \
name.x[i]

#define MD_SET_ATOM_X(name, gid, i, _x) \
name.x[i] = _x

#define MD_GET_ATOM_X_ALL(name, gid) \
name.x

#define MD_GET_ATOM_V(name, gid, i) \
name.v[i]

#define MD_ADD_ATOM_X(name, gid, i, _x) \
name.x[i] += _x

#define MD_SET_ATOM_V(name, gid, i, _v) \
name.v[i] = _v

#define MD_ADD_ATOM_V(name, gid, i, _v) \
name.v[i] += _v

#define MD_GET_ATOM_F(name, gid, i) \
name.f[i]

#define MD_SET_ATOM_F(name, gid, i, _f) \
name.f[i] = _f

#define MD_ADD_ATOM_F(name, gid, i, _f) \
name.f[i] += _f

#define MD_GET_ATOM_RHO(name, gid) \
name.rho

#define MD_SET_ATOM_RHO(name, gid, _rho) \
name.rho = _rho

#define MD_ADD_ATOM_RHO(name, gid, _rho) \
name.rho += _rho

#define MD_GET_ATOM_DF(name, gid) \
name.df

#define MD_SET_ATOM_DF(name, gid, _df) \
name.df = _df

inline _type_atom_list_collection to_atom_list_coll(AtomElement *atoms_prop) {
    _type_atom_list_collection coll{.atoms = atoms_prop};
    return coll;
}

#define TO_ATOM_LIST_COLL(p_list) to_atom_list_coll(p_list->_atoms._data())

#define MD_TO_ATOM_ELEMENT(name, gid) \
name

#endif //MISA_MD_ATOM_PROPS_MACRO_AOS_H
