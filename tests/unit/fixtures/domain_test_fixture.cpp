//
// Created by genshen on 5/12/18.
//

#include <utils/mpi_utils.h>
#include <types/atom_types.h>
#include "domain_test_fixture.h"

comm::BccDomain *DomainFixture::getDomainInstance(int64_t space[3], double lattice_const, double cutoff_radius) {
    MPI_Comm mpi_comm; // new domain
    comm::mpi_process m_process{};
    m_process.comm = MPI_COMM_WORLD;
    MPI_Comm_size(m_process.comm, &m_process.all_ranks);
    MPI_Comm_rank(m_process.comm, &m_process.own_rank);

    comm::BccDomain *ptr_domain = comm::BccDomain::Builder()
            .setComm(m_process, &mpi_comm)
            .setPhaseSpace(space)
            .setCutoffRadius(cutoff_radius)
            .setLatticeConst(lattice_const)
            .build();
    kiwi::mpiUtils::onGlobalCommChanged(mpi_comm); // notice MPI changed
    return ptr_domain;
}

void DomainFixture::SetUp() {
    this->p_domain = getDomainInstance(const_cast<int64_t *>(space), lattice_const, cutoff_radius_factor);
    // create mass
    std::vector<double> masses = {55.845, 63.546, 58.6934};
    atom_type::setGlobalAtomMasses(masses);
}
