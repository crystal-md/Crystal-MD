//
// Created by wangxinwei on 2019-07-14.
//

#include <cmath>
#include <random>
#include <map>
#include <gtest/gtest.h>
#include <types/pre_define.h>
#include <lattice/ws_utils.h>
#include <atom/neighbour_index.h>
#include <atom/neighbour_iterator.h>


// used for reading protected member.
class NeiIndexTests : public NeighbourIndex<_type_neighbour_index_ele> {
public:
    explicit NeiIndexTests(AtomList &atom_list)
            : NeighbourIndex<_type_neighbour_index_ele>(atom_list._atoms._data(), atom_list.lattice) {}

    FRIEND_TEST(atom_cutoff_radius_included_test, lattice_included_test);
};

// used to save atoms pair which are in cutoff radius range.
struct CutoffAtom {
    std::array<_type_atom_index, DIMENSION> coords;
    std::array<double, DIMENSION> position;
};

bool verify_offset_pairs(std::vector<std::pair<CutoffAtom, CutoffAtom>> &offset_pairs,
                         std::vector<NeiOffset> &nei_offsets, BccLattice lat, const double lattice_const) {
    //assert whether 2 atoms is in nei_index
    bool is_included = false;

    for (auto atomIt : offset_pairs) {
        for (auto it : nei_offsets) {
            if (lat.IndexOf3DIndex(atomIt.first.coords[0], atomIt.first.coords[1], atomIt.first.coords[2])
                + it == lat.IndexOf3DIndex(atomIt.second.coords[0], atomIt.second.coords[1], atomIt.second.coords[2])) {
                is_included = true;
                break;
            }
        }

        // log failed pair
        if (!is_included) {
            std::cout << "failed pair offset : " << std::endl;
            std::cout << atomIt.first.coords[0] << " "
                      << atomIt.first.coords[1] << " "
                      << atomIt.first.coords[2] << std::endl;
            std::cout << atomIt.second.coords[0] << " "
                      << atomIt.second.coords[1] << " "
                      << atomIt.second.coords[2] << std::endl;

            std::cout << "failed pair position : " << std::endl;
            std::cout << atomIt.first.position[0] / lattice_const << " "
                      << atomIt.first.position[1] / lattice_const << " "
                      << atomIt.first.position[2] / lattice_const << std::endl;
            std::cout << atomIt.second.position[0] / lattice_const << " "
                      << atomIt.second.position[1] / lattice_const << " "
                      << atomIt.second.position[2] / lattice_const << std::endl;
            return false;
        }
    }
    return true;
}

// we test that:
// we generate 2 atoms A and B with random position, if A B are in cutoff radius range,
// they will be recorded.
// On the other hand, A's lattice coordinate(integer value) can be verified by B's coordinate plus a offset,
// where the offset is generated by function NeighbourIndex::make.
TEST(atom_cutoff_radius_included_test, lattice_included_test) {
    //initial
    const int grid_size[3] = {2, 2, 2};
    const int grid_coord[3] = {0, 0, 0};
    const int64_t space[3] = {25 * grid_size[0], 25 * grid_size[1], 25 * grid_size[2]};
    const double lattice_const = 2.85532;
    const double cutoff_radius_factor = 1.96125;
    const double cutoff_radius = cutoff_radius_factor * lattice_const;

    comm::Domain *p_domain = comm::Domain::Builder()
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius_factor)
            .setGhostSize(static_cast<_type_atom_index>(ceil(cutoff_radius_factor)) + 1)
            .setLatticeConst(lattice_const)
            .localBuild(grid_size, grid_coord);
    AtomList atom_list(p_domain->ghost_extended_lattice_size[0] * 2,
                       p_domain->ghost_extended_lattice_size[1],
                       p_domain->ghost_extended_lattice_size[2],
                       p_domain->sub_box_lattice_size[0] * 2,
                       p_domain->sub_box_lattice_size[1],
                       p_domain->sub_box_lattice_size[2],
                       p_domain->lattice_size_ghost[0] * 2,
                       p_domain->lattice_size_ghost[1],
                       p_domain->lattice_size_ghost[2]);

    //initial offsets
    NeiIndexTests nei_index(atom_list);
    nei_index.make(p_domain->cut_lattice, p_domain->cutoff_radius_factor);

    //random generate
    const _type_lattice_size max = 5000, min = -5000;
    std::default_random_engine e(time(NULL));
    std::uniform_int_distribution<_type_lattice_size> u(min, max);

    for (int l = 0; l < 512; l++) {
        //random atom pair whose distance is short than cutt off
        std::vector<std::pair<CutoffAtom, CutoffAtom>> offset_pairs_odd;
        std::vector<std::pair<CutoffAtom, CutoffAtom>> offset_pairs_even;

        for (int i = 0; i < 2000; i++) {
            //原子位置
            AtomElement src_atom1;
            AtomElement src_atom2;
            src_atom1.x[0] = u(e) / 1000.0;
            src_atom1.x[1] = u(e) / 1000.0;
            src_atom1.x[2] = u(e) / 1000.0;

            src_atom2.x[0] = u(e) / 1000.0;
            src_atom2.x[1] = u(e) / 1000.0;
            src_atom2.x[2] = u(e) / 1000.0;

            // two atoms is in the cutoff radius range.
            if (((src_atom1.x[0] - src_atom2.x[0]) * (src_atom1.x[0] - src_atom2.x[0]) +
                 (src_atom1.x[1] - src_atom2.x[1]) * (src_atom1.x[1] - src_atom2.x[1]) +
                 (src_atom1.x[2] - src_atom2.x[2]) * (src_atom1.x[2] - src_atom2.x[2])) <
                cutoff_radius * cutoff_radius) {

                _type_atom_index coords1[DIMENSION];
                _type_atom_index coords2[DIMENSION];
                ws::getNearLatSubBoxCoord(src_atom1, p_domain, coords1);
                ws::getNearLatSubBoxCoord(src_atom2, p_domain, coords2);

                //to save offset of 2 atoms
                const CutoffAtom atom_meta1{
                        {coords1[0],     coords1[1],     coords1[2]},
                        {src_atom1.x[0], src_atom1.x[1], src_atom1.x[2]},
                };
                const CutoffAtom atom_meta2{
                        {coords2[0],     coords2[1],     coords2[2]},
                        {src_atom2.x[0], src_atom2.x[1], src_atom2.x[2]},
                };

                // filter atoms with the same coords.
                if (!(atom_meta1.coords[0] == atom_meta2.coords[0] &&
                      atom_meta1.coords[1] == atom_meta2.coords[1] &&
                      atom_meta1.coords[2] == atom_meta2.coords[2])) {
                    if (coords1[0] % 2 == 0) { // even
                        offset_pairs_even.emplace_back(std::pair<CutoffAtom, CutoffAtom>(atom_meta1, atom_meta2));
                    } else { // odd
                        offset_pairs_odd.emplace_back(std::pair<CutoffAtom, CutoffAtom>(atom_meta1, atom_meta2));
                    }
                }
            }
        }
        static unsigned int loops = 0;
        loops++;
        std::cout << loops << ": " << offset_pairs_odd.size() << ", " << offset_pairs_even.size() << std::endl;

        const bool passed_odd = verify_offset_pairs(offset_pairs_odd, nei_index.nei_odd_offsets,
                                                    atom_list.lattice, lattice_const);
        EXPECT_TRUE(passed_odd);
        if (!passed_odd) {
            GTEST_FAIL();
        }

        const bool passed_even = verify_offset_pairs(offset_pairs_even, nei_index.nei_even_offsets,
                                                     atom_list.lattice, lattice_const);
        EXPECT_TRUE(passed_even);
        if (!passed_even) {
            GTEST_FAIL();
        }
    }
}
