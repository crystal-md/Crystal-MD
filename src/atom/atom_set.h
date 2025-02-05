//
// Created by genshen on 2018-12-31.
//

#ifndef MISA_MD_ATOM_SET_H
#define MISA_MD_ATOM_SET_H

#include <vector>

#include <comm/domain/domain.h>

#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "atom/neighbour_index.h"
#include "../types/atom_types.h"

/**
 * the container of latticed atom list and inter atom list (all atoms in this sub-box).
 *
 */
class AtomSet {
public:
    typedef NeighbourIndex<_type_neighbour_index_ele> AtomNei;
public:
    /**
     * create atoms list(include lattice atoms list and inter atoms list) and neighbour relative index.
     * @param cutoff_radius cutoff radius.
     * @param extended_lattice_size lattice size with ghost extended.
     * @param sub_box_lattice_size lattice size of current sub box.
     * @param ghost_lattice_size ghost lattice size in each dimension.
     */
    AtomSet(const double cutoff_radius,
            const _type_lattice_size extended_lattice_size[DIMENSION],
            const _type_lattice_size sub_box_lattice_size[DIMENSION],
            const _type_lattice_size ghost_lattice_size[DIMENSION]);

    ~AtomSet();

    inline AtomList *getAtomList() {
        return atom_list;
    }

    inline AtomList &getAtomListRef() {
        return *atom_list;
    }

    inline InterAtomList *getInterList() {
        return inter_atom_list;
    }

    _type_atom_count getnlocalatom(comm::Domain *p_domain);

    /**
     * compute the index offset of neighbour atoms.
     */
    void calcNeighbourIndices(const double cutoff_radius_factor, const _type_lattice_size cut_lattice);

    /**
     * used in read creating mode.
     * @return true for atom is added into the sub-box of current process. false for otherwise.
     */
    bool addAtom(comm::BccDomain *p_domain, const AtomElement atom);

#ifdef MD_RUNTIME_CHECKING

    /**
     * @return the count of atoms that are not vacancy (not type of INVALID).
     */
    _type_atom_count realAtoms() const;

#endif

    inline NeighbourIndex<_type_neighbour_index_ele> *getNeiOffsets() {
        return neighbours;
    }

public:
    AtomList *atom_list;
    InterAtomList *inter_atom_list;

protected:

    double _cutoffRadius;
//    int _cutlattice;
    //   double _latticeconst;
    NeighbourIndex<_type_neighbour_index_ele> *neighbours;
};


#endif //MISA_MD_ATOM_SET_H
