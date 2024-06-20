//
// Created by genshen on 2019/10/4.
//

#include <gtest/gtest.h>
#include <atom.h>
#include <world_builder.h>

#include "fixtures/world_builder_test_fixture.h"

// fixme in mpi, may crash. index out of range
TEST_F(WorldBuilderTestFixture, atom_set_v_test) {
    // test case 1: change atom velocity and check it.
    const _type_lattice_coord pos[3]{25, 25, 25};
    const double dir[3]{1, 3, 5};
    const _type_atom_index gid1 = _atom->atom_list->_atoms.getAtomIndexInSubBox(2 * pos[0], pos[1], pos[2]);
    MD_LOAD_ATOM_VAR(atom_1, _atom->atom_list, gid1);
    AtomElement my_atom_before = MD_TO_ATOM_ELEMENT(atom_1, gid1);

    const _type_lattice_coord pos2[4]{25, 25, 25, 0};
    _atom->setv(pos2, dir, 5000);

    const _type_atom_index gid2 = _atom->atom_list->_atoms.getAtomIndexInSubBox(2 * pos[0], pos[1], pos[2]);
    MD_LOAD_ATOM_VAR(atom_2, _atom->atom_list, gid2);
    AtomElement my_atom_post = MD_TO_ATOM_ELEMENT(atom_2, gid2);

    EXPECT_FLOAT_EQ(my_atom_post.v[0], my_atom_before.v[0] + 222.17975);
    EXPECT_FLOAT_EQ(my_atom_post.v[1], my_atom_before.v[1] + 666.53928);
    EXPECT_FLOAT_EQ(my_atom_post.v[2], my_atom_before.v[2] + 1110.89879);

    // reset
    const _type_atom_index gid_temp = _atom->atom_list->_atoms.getAtomIndexInSubBox(2 * pos[0], pos[1], pos[2]);
    MD_LOAD_ATOM_VAR(atom_tmp, _atom->atom_list, gid_temp);
    MD_SET_ATOM_V(atom_tmp, gid_temp, 0, 0.0);
    MD_SET_ATOM_V(atom_tmp, gid_temp, 1, 0.0);
    MD_SET_ATOM_V(atom_tmp, gid_temp, 2, 0.0);

    // test case 2
    const double dir2[3]{1, 2, 2};
    const _type_atom_index gid_bf2 = _atom->atom_list->_atoms.getAtomIndexInSubBox(2 * pos[0], pos[1], pos[2]);
    MD_LOAD_ATOM_VAR(my_atom_bf_2, _atom->atom_list, gid_bf2);
    AtomElement my_atom_before_2 = MD_TO_ATOM_ELEMENT(my_atom_bf_2, gid_bf2);
    _atom->setv(pos2, dir2, 20000);
    const _type_atom_index gid_pt_2 = _atom->atom_list->_atoms.getAtomIndexInSubBox(2 * pos[0], pos[1], pos[2]);
    MD_LOAD_ATOM_VAR(my_atom_post_2, _atom->atom_list, gid_pt_2);

    EXPECT_FLOAT_EQ(MD_GET_ATOM_V(my_atom_post_2, gid_pt_2, 0), my_atom_before_2.v[0] + 876.28882);
    EXPECT_FLOAT_EQ(MD_GET_ATOM_V(my_atom_post_2, gid_pt_2, 1), my_atom_before_2.v[1] + 1752.5776);
    EXPECT_FLOAT_EQ(MD_GET_ATOM_V(my_atom_post_2, gid_pt_2, 2), my_atom_before_2.v[2] + 1752.5776);
}
