//
// Created by baihe back to 2016-01-06.
//

#ifndef MISA_MD_ATOM_H
#define MISA_MD_ATOM_H

#include <cstdio>
#include <vector>
#include <eam.h>

#include <comm/domain/bcc_domain.h>
#include <logs/logs.h>

#include "atom/atom_set.h"
#include "atom/atom_element.h"
#include "atom/atom_list.h"
#include "atom/inter_atom_list.h"
#include "pack/send_recv_lists.h"
#include "pack/particledata.h"
#include "pack/lat_particle_data.h"

class atom : public AtomSet {
public :
    atom(comm::BccDomain *domain);

    ~ atom();

    /**
     * move atoms to inter-atom list if the atoms is not in its lattice.
     * @return n_flag
     */
    int decide();

    void clearForce();

    /**
     * @tparam POT_TYPE the potential type used for calculation.
     * @tparam CALC_SYS_POT_ENERGY enable/disable the system potential energy calculation.
     * @param pot pointer of eam potential object
     * @param comm timer for recording communication time.
     */
    template<int POT_TYPE, bool CALC_SYS_POT_ENERGY>
    void computeEam(eam *pot, double &comm);

    /**
     * wrapper function to call template function computeEam
     * @param pot_type potential type
     * @param calc_pot_energy enable/disable the system potential energy calculation
     * @param pot pointer of eam potential object
     * @param comm timer for recording communication time.
     */
    void computeEamWrapper(const unsigned short pot_type, bool calc_pot_energy, eam *pot, double &comm);

    /**
     * set velocity of a atom whose position is specified by array @param lat
     * This atom is called PKA (Primary Knock-on Atom).
     * We first compute the the velocity in each direction from the energy, and the set velocity of the atom.
     * The unit of energy is ev, the unit of velocity is 100m/s.
     *
     * @param lat the position of pka.
     * @param direction the vector of velocity to be set.
     * @param energy the energy of PKA.
     */
    void setv(const _type_lattice_coord lat[4], const double direction[3], const double energy);

public:
    SendRecvLists *p_send_recv_list;

    /**
     * set velocity for a atom whose position is specified by array @param lat.
     * @param lat_x, lat_y, lat_z lattice position of the atom (x dimension is not doubled), whose velocity will be changed.
     * @param v the velocity value at x, y, z dimension.
     */
    void setv(const _type_lattice_coord lat_x, const _type_lattice_coord lat_y,
              const _type_lattice_coord lat_z, const double v[DIMENSION]);


    /**
     * Note: only master process can get the correct value of system potential energy.
     * @return system potential energy.
     */
    inline double get_system_pot_energy() const {
        double energy_data_src[2] = {system_total_embed, system_total_pair};
        double energy_data_dest[2] = {0.0, 0.0};
        MPI_Reduce(energy_data_src, energy_data_dest, 2, MPI_DOUBLE, MPI_SUM, MASTER_PROCESSOR, MPI_COMM_WORLD);
        return energy_data_dest[1] / 2 + energy_data_dest[0];
    }

private:
    comm::BccDomain *p_domain;

    double system_total_embed = 0.0;
    double system_total_pair = 0.0;

    /**
     * calculate electron density for all lattice atoms.
     * @tparam POT_TYPE the potential type used for calculation.
     * @param pot pointer of eam potential object.
     */
    template<int POT_TYPE>
    void latRho(eam *pot);

    /**
     * calculate electron density for all interstitial atoms.
     * @tparam POT_TYPE the potential type used for calculation.
     * @param pot pointer of eam potential object.
     */
    template<int POT_TYPE>
    void interRho(eam *pot);

    /**
     * calculate derivative of embedded energy for all lattice atoms.
     * @tparam WITH_ENERGY this template parameter determines whether or not
     * the embedded energy is calculated.
     * @param pot pointer of eam potential object.
     */
    template<bool WITH_ENERGY = false>
    void latDf(eam *pot);

    /**
     * calculate force (or also with pair potential) for all lattice atoms.
     * @tparam WITH_ENERGY this template parameter determines whether or not
     * the pair potential energy is calculated.
     * @param pot pointer of eam potential object.
     */
    template<bool WITH_ENERGY = false>
    void latForce(eam *pot);

    /**
     * calculate force for all interstitial atoms.
     * @tparam WITH_ENERGY this template parameter determines whether or not
     * the pair potential energy is calculated.
     * @param pot pointer of eam potential object.
     */
    template<bool WITH_ENERGY = false>
    void interForce(eam *pot);
};

#endif // MISA_MD_ATOM_H
