#include <cmath>
#include <string>

#include <utils/mpi_utils.h>
#include <logs/logs.h>
#include <eam.h>
#include <parser/setfl_parser.h>
#include <comm/domain/domain.h>

#include "simulation.h"
#include "input.h"
#include "utils/mpi_domain.h"
#include "arch/hardware_accelerate.hpp"
#include "world_builder.h"

simulation::simulation() :
        _p_domain(nullptr), _atom(nullptr),
        _newton_motion(nullptr), _pot(nullptr) {
    runtime_status.flag_calc_system_potential_energy = false;
//    createDomainDecomposition();
//    collision_step = -1;
}

simulation::~simulation() {
    delete _atom;
    delete _newton_motion;
    delete _pot;
}

void simulation::createDomain(const int64_t phase_space[DIMENSION],
                              const double lattice_const, const double cutoff_radius_factor) {

    //进行区域分解
    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initializing GlobalDomain decomposition.\n");
    MPI_Comm new_comm;
    comm::mpi_process pro = comm::mpi_process{
            MPIDomain::sim_processor.own_rank,
            MPIDomain::sim_processor.all_ranks,
            MPIDomain::sim_processor.comm,
    };
    // In domain creation, we set ghost size as (cut_lattice +1) to avoid neighbor index overflowing.
    _p_domain = comm::BccDomain::Builder()
            .setComm(pro, &new_comm)
            .setPhaseSpace(phase_space)
            .setCutoffRadius(cutoff_radius_factor)
            .setLatticeConst(lattice_const)
            .setGhostSize(static_cast<int>(ceil(cutoff_radius_factor)) + 1)
            .build();
    kiwi::mpiUtils::onGlobalCommChanged(new_comm); // set new domain.
    MPIDomain::sim_processor = kiwi::mpiUtils::global_process;

    // init domain for architectures calculation.
    if (isArchAccSupport()) {
        archAccDomainInit(_p_domain);
    }
    kiwi::logs::v(MASTER_PROCESSOR, "domain", "Initialization done.\n");
}

void simulation::createAtoms(const int64_t phase_space[DIMENSION], const double lattice_const,
                             const double init_step_len, const bool create_mode,
                             const double t_set, const unsigned long create_seed,
                             const std::vector<tp_atom_type_weight> &types_weight,
                             const std::string read_inp_path) {
    _atom = new atom(_p_domain);
    // establish index offset for neighbour.
    _atom->calcNeighbourIndices(_p_domain->cutoff_radius_factor, _p_domain->cut_lattice);

    // init domain and neighbor offset indexes for architectures calculation.
    if (isArchAccSupport()) {
        archAccNeiOffsetInit(_atom->getNeiOffsets());
    }
    if (create_mode) {  //创建原子坐标、速度信息
        WorldBuilder mWorldBuilder;
        mWorldBuilder.setDomain(_p_domain)
                .setAtomsContainer(_atom)
                .setBoxSize(phase_space[0], phase_space[1], phase_space[2])
                .setRandomSeed(create_seed)
                .setLatticeConst(lattice_const)
                .setTset(t_set)
                .setAlloyRatio(types_weight)
                .build();
    } else { //读取原子坐标、速度信息
        input _input;
        _input.readPhaseSpace(read_inp_path, _atom, _p_domain);
    }
    _newton_motion = new NewtonMotion(init_step_len, atom_type::num_atom_types); // time step length.
}

void simulation::prepareForStart(const unsigned short potentialType, const std::string pot_file_path) {
    double starttime, stoptime;
    double commtime, computetime, comm;

    // todo file type funl support. pConfigVal->potentialFileFormat
    // 读取势函数文件
    //atom_type::_type_atom_types eles = 0;
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        char tmp[4096];
        sprintf(tmp, "%s", pot_file_path.c_str());

        FILE *pot_file = fopen(tmp, "r");
        if (pot_file == nullptr) { // todo open too many in md.
            kiwi::logs::e("pot", "open potential file {} failed.\n", pot_file_path);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return;
        }
        // new parser
        SetflParser *parser = new SetflParser(pot_file_path); // todo delete (vector)
        parser->parseHeader(); // elements count got. // todo parsing error.
        // eles = parser->getEles(); // elements values on non-root processors are 0.
        _pot = eam::newInstance(potentialType, parser->getEles(),
                                MASTER_PROCESSOR,
                                MPIDomain::sim_processor.own_rank,
                                MPIDomain::sim_processor.comm);
        // read data
        parser->parseBody(_pot); // todo parsing error.
        parser->done();
/*      BCast type list, the lists size is the same as eam potential elements size.
        parser->type_lists.sync(MASTER_PROCESSOR,
                                MPIDomain::sim_processor.own_rank,
                                MPIDomain::sim_processor.comm,
                                _pot->geEles());
*/
    } else {
        _pot = eam::newInstance(potentialType, 0,
                                MASTER_PROCESSOR,
                                MPIDomain::sim_processor.own_rank,
                                MPIDomain::sim_processor.comm);
    }
    // BCast Potential
    _pot->eamBCast(MASTER_PROCESSOR,
                   MPIDomain::sim_processor.own_rank,
                   MPIDomain::sim_processor.comm);
    _pot->interpolateFile(); // interpolation.

    archAccPotInit(_pot); // it runs after atom and boxes creation, but before simulation running.

    starttime = MPI_Wtime();
    // todo: make _cut_lattice a member of class AtomList
    // borderInter is only for [read atom mode], because in [create atom mode], inter is empty at first step;
    _atom->p_send_recv_list->borderInter(_p_domain);
    _atom->p_send_recv_list->exchangeAtomFirst(_p_domain);
    stoptime = MPI_Wtime();
    commtime = stoptime - starttime;

    _atom->clearForce(); // clear force before running simulation.
    starttime = MPI_Wtime();
    _atom->computeEamWrapper(potentialType, runtime_status.flag_calc_system_potential_energy, _pot, comm);
    stoptime = MPI_Wtime();
    computetime = stoptime - starttime - comm;
    commtime += comm;


    kiwi::logs::i(MASTER_PROCESSOR, "sim", "first step comm time: {}\n", commtime);
    kiwi::logs::i(MASTER_PROCESSOR, "sim", "first step compute time: {}\n", computetime);
}

void simulation::simulate(const unsigned short potentialType, const unsigned long steps, const unsigned long init_step) {
    double starttime, stoptime;
    double commtime = 0.0, computetime = 0.0, comm = 0.0;
    double alltime, allstart, allstop;

    allstart = MPI_Wtime();
    onSimulationStarted(init_step);

    // start simulation
    for (_simulation_time_step = init_step; _simulation_time_step < steps; _simulation_time_step++) {
        beforeStep(potentialType, _simulation_time_step);
        //先进行求解牛顿运动方程第一步
        _newton_motion->firststep(_atom->getAtomList(), _atom->getInterList());

        //判断是否有粒子跑出晶格点
        _atom->decide();

        //通信ghost区域，交换粒子
        starttime = MPI_Wtime();
        _atom->p_send_recv_list->exchangeInter(_p_domain);
        _atom->p_send_recv_list->borderInter(_p_domain);
        _atom->p_send_recv_list->exchangeAtom(_p_domain);
        stoptime = MPI_Wtime();
        commtime += stoptime - starttime;

        //计算力
        _atom->clearForce();
        starttime = MPI_Wtime();
        _atom->computeEamWrapper(potentialType, runtime_status.flag_calc_system_potential_energy, _pot, comm);
        stoptime = MPI_Wtime();
        computetime += stoptime - starttime - comm;
        commtime += comm;

        onForceSolved(_simulation_time_step);

        //求解牛顿运动方程第二步
        _newton_motion->secondstep(_atom->getAtomList(), _atom->getInterList());

        postStep(_simulation_time_step);
    }
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("simulation", "loop comm time: {}.\n", commtime);
        kiwi::logs::i("simulation", "loop compute time: {}.\n", computetime);
    }
    onSimulationDone(_simulation_time_step);

    allstop = MPI_Wtime();
    alltime = allstop - allstart;
    if (MPIDomain::sim_processor.own_rank == MASTER_PROCESSOR) {
        kiwi::logs::i("simulation", "total time: {}.\n", alltime);
    }
}

void simulation::collisionStep(const unsigned short potentialType, unsigned long coll_step, const _type_lattice_coord coll_lat[DIMENSION + 1],
                               const double coll_dir[DIMENSION], const double coll_pka_energy) {
    double comm = 0;
    _atom->setv(coll_lat, coll_dir, coll_pka_energy);
    _atom->p_send_recv_list->exchangeInter(_p_domain);
    _atom->p_send_recv_list->borderInter(_p_domain);
    _atom->p_send_recv_list->exchangeAtom(_p_domain);
    _atom->clearForce();
    _atom->computeEamWrapper(potentialType, runtime_status.flag_calc_system_potential_energy, _pot, comm);
}

void simulation::velocitySetStep(const comm::Region<long> global_region, const double velocity_value[DIMENSION]) {
    for (long z = global_region.z_low; z < global_region.z_high; z++) {
        for (long y = global_region.y_low; y < global_region.y_high; y++) {
            for (long x = global_region.x_low; x < global_region.x_high; x++) {
                _atom->setv(x, y, z, velocity_value);
            }
        }
    }
}

void simulation::finalize() {
    if (_p_domain != nullptr) {
        delete _p_domain;
        _p_domain = nullptr;
    }
}

void simulation::abort(int exitcode) {
    MPI_Abort(MPI_COMM_WORLD, exitcode);
}
