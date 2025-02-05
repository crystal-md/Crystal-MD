//
// Created by genshen on 5/8/18.
//

#include <comm/preset/comm_forwarding_region.h>

#include "arch/hardware_accelerate.hpp"
#include "atom_list.h"

AtomList::AtomList(_type_atom_count size_x, _type_atom_count size_y, _type_atom_count size_z,
                   _type_atom_count size_sub_box_x, _type_atom_count size_sub_box_y, _type_atom_count size_sub_box_z,
                   _type_atom_count ghost_count_x, _type_atom_count ghost_count_y, _type_atom_count ghost_count_z) :
        lattice(size_x, size_y, size_z, size_sub_box_x, size_sub_box_y, size_sub_box_z,
                ghost_count_x, ghost_count_y, ghost_count_z), MD_HASH_LIST_INIT(lattice) {
}

AtomList::~AtomList() {
    MD_HASH_LIST_DESTROY();
}

bool AtomList::isBadList(comm::Domain domain) {
    const double delta = 10;
    for (long z = lattice.purge_ghost_count_z; z < lattice._size_sub_box_z + lattice.purge_ghost_count_z; z++) {
        for (long y = lattice.purge_ghost_count_y; y < lattice._size_sub_box_y + lattice.purge_ghost_count_y; y++) {
            for (long x = lattice.purge_ghost_count_x; x < lattice._size_sub_box_x + lattice.purge_ghost_count_x; x++) {
                const _type_atom_index gid = _atoms.getAtomIndex(x, y, z); // todo long type
                MD_LOAD_ATOM_VAR(_atom, this, gid);
                if (MD_GET_ATOM_X(_atom, gid, 0) < domain.meas_global_region.x_low - delta ||
                    MD_GET_ATOM_X(_atom, gid, 1) < domain.meas_global_region.y_low - delta ||
                    MD_GET_ATOM_X(_atom, gid, 2) < domain.meas_global_region.z_low - delta ||
                    MD_GET_ATOM_X(_atom, gid, 0) > domain.meas_global_region.x_high + delta ||
                    MD_GET_ATOM_X(_atom, gid, 1) > domain.meas_global_region.y_high + delta ||
                    MD_GET_ATOM_X(_atom, gid, 2) > domain.meas_global_region.z_high + delta) {
                    return true;
                }
            }
        }
    }
    return false;
}

/**
 *******************************************
 * following (to end of this file) are the
 * implementations about atom list iterator.
 *******************************************
 */
AtomList::iterator AtomList::begin() {
    return AtomList::iterator(nullptr);
}

AtomList::iterator AtomList::end() {
    return AtomList::iterator(nullptr);
}

AtomList::const_iterator AtomList::begin() const {
    return AtomList::const_iterator(nullptr);
}

AtomList::const_iterator AtomList::end() const {
    return AtomList::const_iterator(nullptr);
}

// iterator
AtomList::iterator::iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::iterator::self_type &AtomList::iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::iterator::self_type AtomList::iterator::operator++(int) {
    self_type i = *this;
    ptr_++;
    return i;
}

AtomElement &AtomList::iterator::operator*() {
    return *ptr_;
}

AtomElement *AtomList::iterator::operator->() {
    return ptr_;
}

bool AtomList::iterator::operator==(const AtomList::iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::iterator::operator!=(const AtomList::iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}

void AtomList::iterator::next() {

}

// const iterator
AtomList::const_iterator::const_iterator(AtomElement *ptr) : ptr_(ptr) {}

AtomList::const_iterator::self_type &AtomList::const_iterator::operator++() {
    ptr_++;
    return *this;
}

AtomList::const_iterator::self_type AtomList::const_iterator::operator++(int junk) {
    self_type i = *this;
    ptr_++;
    return i;
}

const AtomElement &AtomList::const_iterator::operator*() {
    return *ptr_;
}

bool AtomList::const_iterator::operator==(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ == rhs.ptr_;
}

bool AtomList::const_iterator::operator!=(const AtomList::const_iterator::self_type &rhs) {
    return ptr_ != rhs.ptr_;
}
