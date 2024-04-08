#include <cmath>

#include <utils/mpi_domain.h>
#include <eam.h>
#include <comm/comm.hpp>

#include "atom.h"
#include "lattice/ws_utils.h"
#include "atom/atom_props_macro_wrapper.h"
#include "md_building_config.h"
#include "pack/rho_packer.h"
#include "pack/force_packer.h"
#include "pack/df_embed_packer.h"
#include "arch/hardware_accelerate.hpp" // use hardware(eg.GPU, MIC,Sunway slave cores.) to achieve calculate accelerating.

atom::atom(comm::BccDomain *domain)
        : AtomSet(domain->lattice_const * domain->cutoff_radius_factor,
                  domain->ghost_extended_lattice_size,
                  domain->sub_box_lattice_size,
                  domain->lattice_size_ghost),
          p_domain(domain) {
    p_send_recv_list = new SendRecvLists(*(this->atom_list), *(this->inter_atom_list));
}

atom::~atom() {
    delete p_send_recv_list;
}

int atom::decide() {
    inter_atom_list->clearGhost();
    int nflag = 0;
    double dist;
    double xtemp, ytemp, ztemp;

    //对本地晶格点原子进行判断，看是否运动为间隙原子
    for (_type_atom_index k = 0; k < p_domain->dbx_sub_box_lattice_size[2]; k++) {
        for (_type_atom_index j = 0; j < p_domain->dbx_sub_box_lattice_size[1]; j++) {
            for (_type_atom_index i = 0; i < p_domain->dbx_sub_box_lattice_size[0]; i++) {
//                kk = atom_list->IndexOf3DIndex(i, j, k);
                const _type_atom_index gid = atom_list->_atoms.getAtomIndexInSubBox(i, j, k); // todo long type
                MD_LOAD_ATOM_VAR(atom_, atom_list, gid);
                if (!MD_IS_ATOM_TYPE_INTER(atom_, gid)) {
                    xtemp = (i + p_domain->dbx_sub_box_lattice_region.x_low) * 0.5 * p_domain->lattice_const;
                    ytemp = (j + p_domain->dbx_sub_box_lattice_region.y_low + (i % 2) * 0.5) *
                            p_domain->lattice_const;
                    ztemp = (k + p_domain->dbx_sub_box_lattice_region.z_low + (i % 2) * 0.5) *
                            p_domain->lattice_const;
                    dist = (MD_GET_ATOM_X(atom_, gid, 0) - xtemp) * (MD_GET_ATOM_X(atom_, gid, 0) - xtemp);
                    dist += (MD_GET_ATOM_X(atom_, gid, 1) - ytemp) * (MD_GET_ATOM_X(atom_, gid, 1) - ytemp);
                    dist += (MD_GET_ATOM_X(atom_, gid, 2) - ztemp) * (MD_GET_ATOM_X(atom_, gid, 2) - ztemp);
                    if (dist > (pow(0.2 * p_domain->lattice_const, 2.0))) { /**超过距离则判断为间隙原子*/
                        inter_atom_list->addInterAtom(MD_TO_ATOM_ELEMENT(atom_, gid));
                        MD_SET_ATOM_TYPE(atom_, gid, atom_type::INVALID);
                        MD_SET_ATOM_V(atom_, gid, 0, 0.0);
                        MD_SET_ATOM_V(atom_, gid, 1, 0.0);
                        MD_SET_ATOM_V(atom_, gid, 2, 0.0);
                        nflag = 1;
                    }
                }
            }
        }
    }

    // 如果间隙原子跑入晶格点,且晶格点为空位, 则空位-间隙发生复合.
    // note: the out-of-lattice atoms can be out of sub-box.
    _type_inter_list::iterator inter_it;
    for (inter_it = inter_atom_list->inter_list.begin(); inter_it != inter_atom_list->inter_list.end();) {
        AtomElement &inter_ref = *inter_it;
        // get the near atom of the inter atom
        // If we use func findNearLatAtom to find a near atom of an inter atom in atoms list
        // the near atom can be an ghost atom (but the position of that ghost atom may still be in sub-box).
        // we should find near atom only in lattice atoms(exclude ghost atoms), so we use func finNearLatAtomInSubBox.
        const _type_atom_index near_atom_inx = ws::findNearLatIndexInSubBox(atom_list->lattice, inter_ref, p_domain);

        // the near atom must be in sub-box, and it is in the lattice atom lists.
        if (near_atom_inx != box::IndexNotExists) {
            MD_LOAD_ATOM_VAR(near_atom, atom_list, near_atom_inx);
            if (MD_IS_ATOM_TYPE_INTER(near_atom, near_atom_inx) &&
                ws::isOutBox(MD_GET_ATOM_X_ALL(near_atom, near_atom_inx), p_domain) == box::IN_BOX) {
                MD_SET_ATOM_ID(near_atom, near_atom_inx, inter_ref.id);
                MD_SET_ATOM_TYPE(near_atom, near_atom_inx, inter_ref.type); // set type to valid.
                MD_SET_ATOM_X(near_atom, near_atom_inx, 0, inter_ref.x[0]);
                MD_SET_ATOM_X(near_atom, near_atom_inx, 1, inter_ref.x[1]);
                MD_SET_ATOM_X(near_atom, near_atom_inx, 2, inter_ref.x[2]);
                MD_SET_ATOM_V(near_atom, near_atom_inx, 0, inter_ref.v[0]);
                MD_SET_ATOM_V(near_atom, near_atom_inx, 1, inter_ref.v[1]);
                MD_SET_ATOM_V(near_atom, near_atom_inx, 2, inter_ref.v[2]);

                // remove this atom from inter list.
                inter_it = inter_atom_list->removeInter(inter_it);
            } else {
                inter_it++;
            }
        } else {
            inter_it++;
        }
    }
    return nflag;
}

void atom::clearForce() {
#ifdef MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS
    for (_type_atom_index i = 0; i < atom_list->cap(); i++) {
        MD_LOAD_ATOM_VAR(atom_, atom_list, i);
        MD_SET_ATOM_F(atom_, i, 0, 0.0);
        MD_SET_ATOM_F(atom_, i, 1, 0.0);
        MD_SET_ATOM_F(atom_, i, 2, 0.0);
        MD_SET_ATOM_RHO(atom_, i, 0.0);
    }
#endif
#ifdef MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_SOA
    // quick set force and rho for SoA memory layout.
    memset(atom_list->atom_f._data(), 0, sizeof(_type_atom_force[DIMENSION]) * atom_list->cap());
    memset(atom_list->atom_rho._data(), 0, sizeof(_type_atom_rho) * atom_list->cap());
#endif

    for (AtomElement &inter_ref: inter_atom_list->inter_list) {
        inter_ref.f[0] = 0.0;
        inter_ref.f[1] = 0.0;
        inter_ref.f[2] = 0.0;
        inter_ref.rho = 0.0;
    }
}

void atom::computeEamWrapper(const unsigned short pot_type, bool calc_pot_energy, eam *pot, double &comm) {
    if (pot_type == EAM_STYLE_ALLOY) {
        if (calc_pot_energy) {
            computeEam<EAM_STYLE_ALLOY, true>(pot, comm);
        } else {
            computeEam<EAM_STYLE_ALLOY, false>(pot, comm);
        }
    } else if (pot_type == EAM_STYLE_FS) {
        if (calc_pot_energy) {
            computeEam<EAM_STYLE_FS, true>(pot, comm);
        } else {
            computeEam<EAM_STYLE_FS, false>(pot, comm);
        }
    } else {
        kiwi::logs::e("eam", "unsupported potential type\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

template<int POT_TYPE, bool CALC_SYS_POT_ENERGY>
void atom::computeEam(eam *pot, double &comm) {
    double starttime, stoptime;
    inter_atom_list->makeIndex(atom_list, p_domain); // create index for inter atom and inter ghost atoms.

    system_total_embed = 0.0;
    system_total_pair = 0.0;

    latRho<POT_TYPE>(pot);
    interRho<POT_TYPE>(pot);

    {
        // 发送电子云密度
        // 将ghost区域的粒子的电子云密度发送给其所在的进程，得到完整的电子云密度
        starttime = MPI_Wtime();
        RhoPacker rho_packer(getAtomListRef(), p_send_recv_list->sendlist, p_send_recv_list->recvlist);
        comm::neiSendReceive<double, true>(&rho_packer, MPIDomain::toCommProcess(),
                                           MPI_DOUBLE, p_domain->rank_id_neighbours);
        stoptime = MPI_Wtime();
        comm = stoptime - starttime;
    }

    //本地晶格点计算嵌入能导数
    latDf<CALC_SYS_POT_ENERGY>(pot);

    {
        // 发送嵌入能导数
        // 将本地box属于邻居进程ghost区域的粒子的嵌入能导数发送给邻居进程
        starttime = MPI_Wtime();
        DfEmbedPacker packer(getAtomListRef(), p_send_recv_list->sendlist, p_send_recv_list->recvlist,
                             p_send_recv_list->intersendlist,
                             p_send_recv_list->interrecvlist);
        comm::neiSendReceive<double>(&packer, MPIDomain::toCommProcess(), MPI_DOUBLE, p_domain->rank_id_neighbours);
        stoptime = MPI_Wtime();
        comm += stoptime - starttime;
    }

    // force for local lattice.
    latForce<CALC_SYS_POT_ENERGY>(pot);

    //间隙原子计算嵌入能和对势带来的力
    interForce<CALC_SYS_POT_ENERGY>(pot);

    // send force
    starttime = MPI_Wtime();
    ForcePacker force_packer(getAtomListRef(), p_send_recv_list->sendlist, p_send_recv_list->recvlist);
    comm::neiSendReceive<double, true>(&force_packer, MPIDomain::toCommProcess(),
                                       MPI_DOUBLE, p_domain->rank_id_neighbours);
    stoptime = MPI_Wtime();
    comm += stoptime - starttime;
}

template void atom::computeEam<EAM_STYLE_ALLOY, true>(eam *pot, double &comm);

template void atom::computeEam<EAM_STYLE_ALLOY, false>(eam *pot, double &comm);

template void atom::computeEam<EAM_STYLE_FS, true>(eam *pot, double &comm);

template void atom::computeEam<EAM_STYLE_FS, false>(eam *pot, double &comm);

template<int POT_TYPE>
void atom::latRho(eam *pot) {
    double delx, dely, delz;
    double dist2;
    const int xstart = p_domain->dbx_lattice_size_ghost[0];
    const int ystart = p_domain->dbx_lattice_size_ghost[1];
    const int zstart = p_domain->dbx_lattice_size_ghost[2];

    // 本地晶格点上的原子计算电子云密度
    if (isArchAccSupport()) {
        archAccEamRhoCalc(pot, TO_ATOM_LIST_COLL(atom_list), _cutoffRadius); // fixme
    } else { // calculate electron density use cpu only.
        for (int k = zstart; k < p_domain->dbx_sub_box_lattice_size[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_sub_box_lattice_size[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_sub_box_lattice_size[0] + xstart; i++) {
                    const _type_atom_index gid = atom_list->_atoms.getAtomIndex(i, j, k);
                    MD_LOAD_ATOM_VAR(atom_central, atom_list, gid);

                    if (MD_IS_ATOM_TYPE_INTER(atom_central, gid)) {
                        continue;
                    }
                    //对晶格点邻居原子遍历
                    // only consider the atoms whose id is bigger than {@var atom_central}, just single side.
                    AtomNei::iterator nei_itl_end = neighbours->end(true, i, j, k);
                    for (AtomNei::iterator nei_itl = neighbours->begin(true, i, j, k);
                         nei_itl != nei_itl_end; ++nei_itl) {
                        const _type_atom_index nei_id = nei_itl.cur_index;
                        MD_LOAD_ATOM_VAR(atom_neighbour, atom_list, nei_id);

                        if (MD_IS_ATOM_TYPE_INTER(atom_neighbour, nei_id)) {
                            continue;
                        }
                        delx = MD_GET_ATOM_X(atom_central, gid, 0) - MD_GET_ATOM_X(atom_neighbour, nei_id, 0);
                        dely = MD_GET_ATOM_X(atom_central, gid, 1) - MD_GET_ATOM_X(atom_neighbour, nei_id, 1);
                        delz = MD_GET_ATOM_X(atom_central, gid, 2) - MD_GET_ATOM_X(atom_neighbour, nei_id, 2);
                        dist2 = delx * delx + dely * dely + delz * delz;
                        if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                            if (POT_TYPE == EAM_STYLE_FS) {
                                MD_ADD_ATOM_RHO(atom_central, gid, pot->chargeDensity(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_central, gid)),
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_neighbour, nei_id)), dist2));
                                MD_ADD_ATOM_RHO(atom_neighbour, nei_id, pot->chargeDensity(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_neighbour, nei_id)),
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_central, gid)), dist2));
                            } else if (POT_TYPE == EAM_STYLE_ALLOY) {
                                MD_ADD_ATOM_RHO(atom_central, gid, pot->chargeDensity(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_neighbour, nei_id)), dist2));
                                MD_ADD_ATOM_RHO(atom_neighbour, nei_id, pot->chargeDensity(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_central, gid)), dist2));
                            }
                        }
                    }
                }
            }
        }
    }
}

template<int POT_TYPE>
void atom::interRho(eam *pot) {
    double delx, dely, delz;
    double dist2;
    double dfEmbed;

    std::vector<_type_atom_index>::iterator neighbourOffsetsIter;

    //间隙原子电子云密度
    _type_atom_index near_atom_index;
    // note that some inter atoms may be not in sub-box, we should exclude them.
    for (_type_inter_list::iterator inter_it = inter_atom_list->inter_list.begin();
         inter_it != inter_atom_list->inter_list.end(); inter_it++) {
        // get index of nearest atom of inter atoms.
        near_atom_index = ws::findNearLatIndexInSubBox(atom_list->lattice, *inter_it, p_domain);
#ifdef MD_RUNTIME_CHECKING
        if (near_atom_index == box::IndexNotExists) {
            assert(false);
            continue; // make sure the inter atoms is in sub box.
            // todo find a good way to filter out-of-box atoms while exchanging inter atoms.
        }
#endif
        MD_LOAD_ATOM_VAR(atom_near, atom_list, near_atom_index);

        delx = (*inter_it).x[0] - MD_GET_ATOM_X(atom_near, near_atom_index, 0);
        dely = (*inter_it).x[1] - MD_GET_ATOM_X(atom_near, near_atom_index, 1);
        delz = (*inter_it).x[2] - MD_GET_ATOM_X(atom_near, near_atom_index, 2);
        dist2 = delx * delx + dely * dely + delz * delz;
        if (POT_TYPE == EAM_STYLE_FS) {
            if (!MD_IS_ATOM_TYPE_INTER(atom_near, near_atom_index) && dist2 < (_cutoffRadius * _cutoffRadius)) {
                inter_it->rho += pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type),
                                                    atom_type::getTypeIdByType(
                                                            MD_GET_ATOM_TYPE(atom_near, near_atom_index)), dist2);
                MD_ADD_ATOM_RHO(atom_near, near_atom_index,
                                pot->chargeDensity(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_near, near_atom_index)),
                                        atom_type::getTypeIdByType((*inter_it).type), dist2));
            }
        } else if (POT_TYPE == EAM_STYLE_ALLOY) {
            if (!MD_IS_ATOM_TYPE_INTER(atom_near, near_atom_index) && dist2 < (_cutoffRadius * _cutoffRadius)) {
                inter_it->rho += pot->chargeDensity(
                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_near, near_atom_index)), dist2);
                MD_ADD_ATOM_RHO(atom_near, near_atom_index,
                                pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type), dist2));
            }
        }

        _type_atom_index x, y, z;
        atom_list->lattice.get3DIndexByLinearIndex(near_atom_index, x, y, z);
        // rho between inter atoms and lattice atoms (use full neighbour index).
        AtomNei::iterator nei_full_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_full_itl_end; ++nei_itl) {
            const _type_atom_index nei_id = nei_itl.cur_index;
            MD_LOAD_ATOM_VAR(lat_nei_atom, atom_list, nei_id); // this is a lattice atom.
            if (!MD_IS_ATOM_TYPE_INTER(lat_nei_atom, nei_id)) {
                delx = (*inter_it).x[0] - MD_GET_ATOM_X(lat_nei_atom, nei_id, 0);
                dely = (*inter_it).x[1] - MD_GET_ATOM_X(lat_nei_atom, nei_id, 1);
                delz = (*inter_it).x[2] - MD_GET_ATOM_X(lat_nei_atom, nei_id, 2);
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    if (POT_TYPE == EAM_STYLE_FS) {
                        (*inter_it).rho += pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type),
                                                              atom_type::getTypeIdByType(
                                                                      MD_GET_ATOM_TYPE(lat_nei_atom, nei_id)), dist2);
                        MD_ADD_ATOM_RHO(lat_nei_atom, nei_id, pot->chargeDensity(
                                atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(lat_nei_atom, nei_id)),
                                atom_type::getTypeIdByType((*inter_it).type), dist2));
                    } else if (POT_TYPE == EAM_STYLE_ALLOY) {
                        (*inter_it).rho += pot->chargeDensity(
                                atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(lat_nei_atom, nei_id)), dist2);
                        MD_ADD_ATOM_RHO(lat_nei_atom, nei_id,
                                        pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type), dist2));
                    }
                }
            }
        }

        // rho contribution of neighbour atoms in the same bucket.
        {
            const inter_map_range inter_map_range = inter_atom_list->inter_map.equal_range(near_atom_index);
            for (inter_map_range_itl bucket_nei_itl = inter_map_range.first;
                 bucket_nei_itl != inter_map_range.second; ++bucket_nei_itl) {
                if (bucket_nei_itl->second->id == inter_it->id) {
                    continue; // can not be itself.
                }
                delx = (*inter_it).x[0] - bucket_nei_itl->second->x[0];
                dely = (*inter_it).x[1] - bucket_nei_itl->second->x[1];
                delz = (*inter_it).x[2] - bucket_nei_itl->second->x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    if (POT_TYPE == EAM_STYLE_FS) {
                        (*inter_it).rho += pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type),
                                                              atom_type::getTypeIdByType(bucket_nei_itl->second->type),
                                                              dist2);
                    } else if (POT_TYPE == EAM_STYLE_ALLOY) {
                        (*inter_it).rho += pot->chargeDensity(
                                atom_type::getTypeIdByType(bucket_nei_itl->second->type), dist2);
                    }
                }
            }
        }

        // rho between inter atoms and inter atoms (use full neighbour index).
        AtomNei::iterator nei_half_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_half_itl_end; ++nei_itl) {
            const _type_atom_index inter_nei_id = nei_itl.cur_index; // get index of the neighbour lattice.
            // get inter atoms on this neighbour lattice and calculate inter-rho.
            inter_map_range inter_map_range = inter_atom_list->inter_map.equal_range(inter_nei_id);
            for (inter_map_range_itl itl = inter_map_range.first; itl != inter_map_range.second; ++itl) {
                delx = (*inter_it).x[0] - itl->second->x[0];
                dely = (*inter_it).x[1] - itl->second->x[1];
                delz = (*inter_it).x[2] - itl->second->x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    if (POT_TYPE == EAM_STYLE_FS) {
                        (*inter_it).rho += pot->chargeDensity(atom_type::getTypeIdByType((*inter_it).type),
                                                              atom_type::getTypeIdByType(itl->second->type), dist2);
                    } else if (POT_TYPE == EAM_STYLE_ALLOY) {
                        (*inter_it).rho += pot->chargeDensity(atom_type::getTypeIdByType(itl->second->type), dist2);
                        //                    itl->second->rho += pot->rhoContribution(atom_type::getTypeIdByType((*inter_it).type), dist2);
                    }
                }
            }
        }
        //计算间隙原子嵌入能导数
        dfEmbed = pot->dEmbedEnergy(atom_type::getTypeIdByType((*inter_it).type), (*inter_it).rho);
        (*inter_it).df = dfEmbed;
    }
}

template<bool WITH_ENERGY>
void atom::latDf(eam *pot) {
    double dfEmbed;
    const int xstart = p_domain->dbx_lattice_size_ghost[0];
    const int ystart = p_domain->dbx_lattice_size_ghost[1];
    const int zstart = p_domain->dbx_lattice_size_ghost[2];

    //本地晶格点计算嵌入能导数
    if (isArchAccSupport()) {
        archAccEamDfCalc(pot, TO_ATOM_LIST_COLL(atom_list), _cutoffRadius);    // fixme
    } else {
        for (int k = zstart; k < p_domain->dbx_sub_box_lattice_size[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_sub_box_lattice_size[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_sub_box_lattice_size[0] + xstart; i++) {
                    const _type_atom_index gid = atom_list->_atoms.getAtomIndex(i, j, k);
                    MD_LOAD_ATOM_VAR(atom_, atom_list, gid);

                    if (MD_IS_ATOM_TYPE_INTER(atom_, gid)) {
                        continue;
                    }
                    dfEmbed = pot->dEmbedEnergy(atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_, gid)),
                                                MD_GET_ATOM_RHO(atom_, gid));
                    MD_SET_ATOM_DF(atom_, gid, dfEmbed);
                    // embed energy
                    if (WITH_ENERGY) {
                        system_total_embed += pot->embedEnergy(atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_, gid)),
                                                               MD_GET_ATOM_RHO(atom_, gid),
                                                               MD_GET_ATOM_RHO(atom_, gid) + 0.1);
                    }
                }
            }
        }
    }
}

template<bool WITH_ENERGY>
void atom::latForce(eam *pot) {
    double delx, dely, delz;
    double dist2;
    double fpair;

    const int xstart = p_domain->dbx_lattice_size_ghost[0];
    const int ystart = p_domain->dbx_lattice_size_ghost[1];
    const int zstart = p_domain->dbx_lattice_size_ghost[2];

    if (isArchAccSupport()) {
        archAccEamForceCalc(pot, TO_ATOM_LIST_COLL(atom_list), _cutoffRadius);
    } else {
        for (int k = zstart; k < p_domain->dbx_sub_box_lattice_size[2] + zstart; k++) {
            for (int j = ystart; j < p_domain->dbx_sub_box_lattice_size[1] + ystart; j++) {
                for (int i = xstart; i < p_domain->dbx_sub_box_lattice_size[0] + xstart; i++) {
                    const _type_atom_index gid = atom_list->_atoms.getAtomIndex(i, j, k);
                    MD_LOAD_ATOM_VAR(atom_, atom_list, gid);

                    if (MD_IS_ATOM_TYPE_INTER(atom_, gid)) {
                        continue;
                    }

                    // force between lattice atoms and lattice atoms.
                    AtomNei::iterator nei_itl_end = neighbours->end(true, i, j, k);
                    for (AtomNei::iterator nei_itl = neighbours->begin(true, i, j, k);
                         nei_itl != nei_itl_end; ++nei_itl) {
                        const _type_atom_index nei_id = nei_itl.cur_index;
                        MD_LOAD_ATOM_VAR(atom_n, atom_list, nei_id);

                        delx = MD_GET_ATOM_X(atom_, gid, 0) - MD_GET_ATOM_X(atom_n, nei_id, 0);
                        dely = MD_GET_ATOM_X(atom_, gid, 1) - MD_GET_ATOM_X(atom_n, nei_id, 1);
                        delz = MD_GET_ATOM_X(atom_, gid, 2) - MD_GET_ATOM_X(atom_n, nei_id, 2);
                        dist2 = delx * delx + dely * dely + delz * delz;
                        if (dist2 < (_cutoffRadius * _cutoffRadius) && !MD_IS_ATOM_TYPE_INTER(atom_n, nei_id)) {
                            fpair = pot->toForce(atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_, gid)),
                                                 atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_n, nei_id)),
                                                 dist2, MD_GET_ATOM_DF(atom_, gid), MD_GET_ATOM_DF(atom_n, nei_id));

                            MD_ADD_ATOM_F(atom_, gid, 0, delx * fpair);
                            MD_ADD_ATOM_F(atom_, gid, 1, dely * fpair);
                            MD_ADD_ATOM_F(atom_, gid, 2, delz * fpair);

                            MD_ADD_ATOM_F(atom_n, nei_id, 0, -delx * fpair);
                            MD_ADD_ATOM_F(atom_n, nei_id, 1, -dely * fpair);
                            MD_ADD_ATOM_F(atom_n, nei_id, 2, -delz * fpair);

                            if (WITH_ENERGY) {
                                system_total_pair += pot->pairPotential(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_, gid)),
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_n, nei_id)),
                                        dist2);
                                system_total_pair += pot->pairPotential(
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_n, nei_id)),
                                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_, gid)),
                                        dist2);
                            }
                        }
                    }
                }
            }
        }
    } // end of if-isArchAccSupport.
}

template<bool WITH_ENERGY>
void atom::interForce(eam *pot) {
    double delx, dely, delz;
    double dist2;
    double fpair;

    //间隙原子计算嵌入能和对势带来的力
    _type_atom_index _atom_near_index;
    for (_type_inter_list::iterator inter_it = inter_atom_list->inter_list.begin();
         inter_it != inter_atom_list->inter_list.end(); inter_it++) {
        _atom_near_index = ws::findNearLatIndexInSubBox(atom_list->lattice, *inter_it, p_domain);
#ifdef MD_RUNTIME_CHECKING
        if (_atom_near_index == box::IndexNotExists) {
            assert(false);
            continue; // make sure the inter atoms is in sub box.
            // todo find a good way to filter out-of-box atoms while exchanging inter atoms.
        }
#endif 
        // 间隙原子所在晶格处的原子
        MD_LOAD_ATOM_VAR(atom_central, atom_list, _atom_near_index);

        delx = (*inter_it).x[0] - MD_GET_ATOM_X(atom_central, _atom_near_index, 0);
        dely = (*inter_it).x[1] - MD_GET_ATOM_X(atom_central, _atom_near_index, 1);
        delz = (*inter_it).x[2] - MD_GET_ATOM_X(atom_central, _atom_near_index, 2);
        dist2 = delx * delx + dely * dely + delz * delz;
        if (dist2 < (_cutoffRadius * _cutoffRadius) && !MD_IS_ATOM_TYPE_INTER(atom_central, _atom_near_index)) {
            fpair = pot->toForce(
                    atom_type::getTypeIdByType((*inter_it).type),
                    atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_central, _atom_near_index)),
                    dist2, (*inter_it).df, MD_GET_ATOM_DF(atom_central, _atom_near_index));

            (*inter_it).f[0] += delx * fpair;
            (*inter_it).f[1] += dely * fpair;
            (*inter_it).f[2] += delz * fpair;

            MD_ADD_ATOM_F(atom_central, _atom_near_index, 0, -delx * fpair);
            MD_ADD_ATOM_F(atom_central, _atom_near_index, 1, -dely * fpair);
            MD_ADD_ATOM_F(atom_central, _atom_near_index, 2, -delz * fpair);

            if (WITH_ENERGY) {
                system_total_pair += pot->pairPotential(
                        atom_type::getTypeIdByType((*inter_it).type),
                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_central, _atom_near_index)),
                        dist2);
                system_total_pair += pot->pairPotential(
                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(atom_central, _atom_near_index)),
                        atom_type::getTypeIdByType((*inter_it).type),
                        dist2);
            }
        }

        _type_atom_index x, y, z;
        atom_list->lattice.get3DIndexByLinearIndex(_atom_near_index, x, y, z);
        // force between inter atoms and lattice atoms (use full neighbour index).
        AtomNei::iterator nei_full_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_full_itl_end; ++nei_itl) {
            const _type_atom_index nei_atom_inx = nei_itl.cur_index;
            MD_LOAD_ATOM_VAR(lattice_neighbour, atom_list, nei_atom_inx); // this is a lattice atom.
            delx = (*inter_it).x[0] - MD_GET_ATOM_X(lattice_neighbour, nei_atom_inx, 0);
            dely = (*inter_it).x[1] - MD_GET_ATOM_X(lattice_neighbour, nei_atom_inx, 1);
            delz = (*inter_it).x[2] - MD_GET_ATOM_X(lattice_neighbour, nei_atom_inx, 2);
            dist2 = delx * delx + dely * dely + delz * delz;
            if (dist2 < (_cutoffRadius * _cutoffRadius) && !MD_IS_ATOM_TYPE_INTER(lattice_neighbour, nei_atom_inx)) {
                fpair = pot->toForce(
                        atom_type::getTypeIdByType((*inter_it).type),
                        atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(lattice_neighbour, nei_atom_inx)),
                        dist2, (*inter_it).df, MD_GET_ATOM_DF(lattice_neighbour, nei_atom_inx));

                (*inter_it).f[0] += delx * fpair;
                (*inter_it).f[1] += dely * fpair;
                (*inter_it).f[2] += delz * fpair;

                MD_ADD_ATOM_F(lattice_neighbour, nei_atom_inx, 0, -delx * fpair);
                MD_ADD_ATOM_F(lattice_neighbour, nei_atom_inx, 1, -dely * fpair);
                MD_ADD_ATOM_F(lattice_neighbour, nei_atom_inx, 2, -delz * fpair);

                if (WITH_ENERGY) {
                    system_total_pair += pot->pairPotential(
                            atom_type::getTypeIdByType((*inter_it).type),
                            atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(lattice_neighbour, nei_atom_inx)),
                            dist2);
                    system_total_pair += pot->pairPotential(
                            atom_type::getTypeIdByType(MD_GET_ATOM_TYPE(lattice_neighbour, nei_atom_inx)),
                            atom_type::getTypeIdByType((*inter_it).type),
                            dist2);
                }
            }
        }
        // force contribution of neighbour atoms in the same bucket.
        {
            const inter_map_range inter_map_range = inter_atom_list->inter_map.equal_range(_atom_near_index);
            for (inter_map_range_itl bucket_nei_itl = inter_map_range.first;
                 bucket_nei_itl != inter_map_range.second; ++bucket_nei_itl) {
                if (bucket_nei_itl->second->id == inter_it->id) {
                    continue; // can not be itself.
                }
                delx = (*inter_it).x[0] - bucket_nei_itl->second->x[0];
                dely = (*inter_it).x[1] - bucket_nei_itl->second->x[1];
                delz = (*inter_it).x[2] - bucket_nei_itl->second->x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    fpair = pot->toForce(
                            atom_type::getTypeIdByType((*inter_it).type),
                            atom_type::getTypeIdByType(bucket_nei_itl->second->type),
                            dist2, (*inter_it).df, bucket_nei_itl->second->df);

                    (*inter_it).f[0] += delx * fpair;
                    (*inter_it).f[1] += dely * fpair;
                    (*inter_it).f[2] += delz * fpair;

                    if (WITH_ENERGY) {
                        system_total_pair += pot->pairPotential(
                                atom_type::getTypeIdByType((*inter_it).type),
                                atom_type::getTypeIdByType(bucket_nei_itl->second->type),
                                dist2);
                    }
                }
            }
        }

        // force between inter atoms and inter atoms(including inter ghost atom) (use full neighbour index).
        AtomNei::iterator nei_half_itl_end = neighbours->end(false, x, y, z);
        for (AtomNei::iterator nei_itl = neighbours->begin(false, x, y, z);
             nei_itl != nei_half_itl_end; ++nei_itl) {
            const _type_atom_index inter_nei_id = nei_itl.cur_index; // get index of the neighbour lattice.
            inter_map_range inter_map_range_up = inter_atom_list->inter_map.equal_range(inter_nei_id);
            for (inter_map_range_itl itl_up = inter_map_range_up.first;
                 itl_up != inter_map_range_up.second; ++itl_up) {
                delx = (*inter_it).x[0] - itl_up->second->x[0];
                dely = (*inter_it).x[1] - itl_up->second->x[1];
                delz = (*inter_it).x[2] - itl_up->second->x[2];
                dist2 = delx * delx + dely * dely + delz * delz;
                if (dist2 < (_cutoffRadius * _cutoffRadius)) {
                    fpair = pot->toForce(
                            atom_type::getTypeIdByType((*inter_it).type),
                            atom_type::getTypeIdByType(itl_up->second->type),
                            dist2, (*inter_it).df, itl_up->second->df);

                    (*inter_it).f[0] += delx * fpair;
                    (*inter_it).f[1] += dely * fpair;
                    (*inter_it).f[2] += delz * fpair;

                    if (WITH_ENERGY) {
                        system_total_pair += pot->pairPotential(
                                atom_type::getTypeIdByType((*inter_it).type),
                                atom_type::getTypeIdByType(itl_up->second->type),
                                dist2);
                    }
                }
            }
        }
    }
}

void atom::setv(const _type_lattice_coord lat[4], const double direction[3], const double energy) {
    long kk;
    if ((lat[0] * 2) >= p_domain->dbx_sub_box_lattice_region.x_low &&
        (lat[0] * 2) < (p_domain->dbx_sub_box_lattice_region.x_low + p_domain->dbx_sub_box_lattice_size[0])
        && lat[1] >= p_domain->dbx_sub_box_lattice_region.y_low &&
        lat[1] < (p_domain->dbx_sub_box_lattice_region.y_low + p_domain->dbx_sub_box_lattice_size[1])
        && lat[2] >= p_domain->dbx_sub_box_lattice_region.z_low &&
        lat[2] < (p_domain->dbx_sub_box_lattice_region.z_low + p_domain->dbx_sub_box_lattice_size[2])) {
        kk = (atom_list->lattice.IndexOf3DIndex(lat[0] * 2 - p_domain->dbx_ghost_ext_lattice_region.x_low,
                                                lat[1] - p_domain->dbx_ghost_ext_lattice_region.y_low,
                                                lat[2] - p_domain->dbx_ghost_ext_lattice_region.z_low) + lat[3]);
        // todo verify the position.
        MD_LOAD_ATOM_VAR(atom_, atom_list, kk);
        // the unit of v is A/ps (or 100m/s)
        const double v_ = sqrt(2 * energy / atom_type::getAtomMass(MD_GET_ATOM_TYPE(atom_, kk)) / mvv2e);
        const double d_ = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
        // fixme: to set v, just set, no adding.
        MD_ADD_ATOM_V(atom_, kk, 0, v_ * direction[0] / d_);
        MD_ADD_ATOM_V(atom_, kk, 1, v_ * direction[1] / d_);
        MD_ADD_ATOM_V(atom_, kk, 2, v_ * direction[2] / d_);
    }
}

void atom::setv(const _type_lattice_coord lat_x, const _type_lattice_coord lat_y,
                const _type_lattice_coord lat_z, const double v[DIMENSION]) {
    if ((lat_x * 2) >= p_domain->dbx_sub_box_lattice_region.x_low &&
        (lat_x * 2) < (p_domain->dbx_sub_box_lattice_region.x_low + p_domain->dbx_sub_box_lattice_size[0])
        && lat_y >= p_domain->dbx_sub_box_lattice_region.y_low &&
        lat_y < (p_domain->dbx_sub_box_lattice_region.y_low + p_domain->dbx_sub_box_lattice_size[1])
        && lat_z >= p_domain->dbx_sub_box_lattice_region.z_low &&
        lat_z < (p_domain->dbx_sub_box_lattice_region.z_low + p_domain->dbx_sub_box_lattice_size[2])) {
        const long kk = atom_list->lattice.IndexOf3DIndex(
                lat_x * 2 - p_domain->dbx_ghost_ext_lattice_region.x_low,
                lat_y - p_domain->dbx_ghost_ext_lattice_region.y_low,
                lat_z - p_domain->dbx_ghost_ext_lattice_region.z_low);

        MD_LOAD_ATOM_VAR(atom_, atom_list, kk);
        MD_SET_ATOM_V(atom_, kk, 0, v[0]);
        MD_SET_ATOM_V(atom_, kk, 1, v[1]);
        MD_SET_ATOM_V(atom_, kk, 2, v[2]);
    }
}
