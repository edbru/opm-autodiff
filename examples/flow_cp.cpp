/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <dune/common/version.hh>

#include <opm/core/utility/platform_dependent/disable_warnings.h>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/common/GridAdapter.hpp>

#include <opm/core/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/cornerpoint_grid.h>
#include <opm/autodiff/GridHelpers.hpp>

#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/initStateEquil.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/thresholdPressures.hpp> // Note: the GridHelpers must be included before this (to make overloads available). \TODO: Fix.

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/ExtractParallelGridInformationToISTL.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/autodiff/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/RedistributeDataHandles.hpp>

#include <opm/core/utility/share_obj.hpp>

#include <opm/parser/eclipse/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/OpmLog/StreamLog.hpp>
#include <opm/parser/eclipse/OpmLog/CounterLog.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <memory>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <numeric>

namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    // Must ensure an instance of the helper is created to initialise MPI.
    const Dune::MPIHelper& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    using namespace Opm;

    std::cout << "******************************************************************************************\n";
    std::cout << "*                                                                                        *\n";
    std::cout << "*     This is Flow (version 2015.04, experimental variant using dune-cornerpoint)        *\n";
    std::cout << "*                                                                                        *\n";
    std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow that is part of OPM. *\n";
    std::cout << "* For more information see:                                                              *\n";
    std::cout << "*                             http://opm-project.org                                     *\n";
    std::cout << "*                                                                                        *\n";
    std::cout << "******************************************************************************************\n\n";

    // Read parameters, see if a deck was specified on the command line.
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;
    parameter::ParameterGroup param(argc, argv, false);
    if (!param.unhandledArguments().empty()) {
        if (param.unhandledArguments().size() != 1) {
            OPM_THROW(std::runtime_error, "You can only specify a single input deck on the command line.");
        } else {
            param.insertParameter("deck_filename", param.unhandledArguments()[0]);
        }
    }

    // We must have an input deck. Grid and props will be read from that.
    if (!param.has("deck_filename")) {
        std::cerr << "This program must be run with an input deck.\n"
            "Specify the deck filename either\n"
            "    a) as a command line argument by itself\n"
            "    b) as a command line parameter with the syntax deck_filename=<path to your deck>, or\n"
            "    c) as a parameter in a parameter file (.param or .xml) passed to the program.";
        OPM_THROW(std::runtime_error, "Input deck required.");
    }
    std::shared_ptr<Dune::CpGrid> grid;
    std::shared_ptr<BlackoilPropertiesInterface> props;
    std::shared_ptr<BlackoilPropsAdFromDeck> new_props;
    std::shared_ptr<RockCompressibility> rock_comp;
    BlackoilState state;
    // bool check_well_controls = false;
    // int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    std::string deck_filename = param.get<std::string>("deck_filename");

    // Write parameters used for later reference. (only if rank is zero)
    bool output = ( mpi_helper.rank() == 0 ) && param.getDefault("output", true);
    std::string output_dir;
    if (output) {
        // Create output directory if needed.
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        // Write simulation parameters.
        param.writeParam(output_dir + "/simulation.param");
    }

    std::string logFile = output_dir + "/LOGFILE.txt";
    Opm::ParserPtr parser(new Opm::Parser());
    {
        std::shared_ptr<Opm::StreamLog> streamLog = std::make_shared<Opm::StreamLog>(logFile , Opm::Log::DefaultMessageTypes);
        std::shared_ptr<Opm::CounterLog> counterLog = std::make_shared<Opm::CounterLog>(Opm::Log::DefaultMessageTypes);

        Opm::OpmLog::addBackend( "STREAM" , streamLog );
        Opm::OpmLog::addBackend( "COUNTER" , counterLog );
    }

    Opm::DeckConstPtr deck;
    std::shared_ptr<EclipseState> eclipseState;
    try {
        deck = parser->parseFile(deck_filename);
        Opm::checkDeck(deck);
        eclipseState.reset(new Opm::EclipseState(deck));
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Failed to create valid ECLIPSESTATE object. See logfile: " << logFile << std::endl;
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Grid init
    grid.reset(new Dune::CpGrid);
    std::vector<double> porv = eclipseState->getDoubleGridProperty("PORV")->getData();
    grid->processEclipseFormat(deck, false, false, false, porv);

    const PhaseUsage pu = Opm::phaseUsageFromDeck(deck);
    Opm::BlackoilOutputWriter outputWriter(*grid, param, eclipseState, pu );

    // Rock and fluid init
    props.reset(new BlackoilPropertiesFromDeck(deck, eclipseState,
                                               Opm::UgGridHelpers::numCells(*grid),
                                               Opm::UgGridHelpers::globalCell(*grid),
                                               Opm::UgGridHelpers::cartDims(*grid),
                                               Opm::UgGridHelpers::beginCellCentroids(*grid),
                                               Opm::UgGridHelpers::dimensions(*grid), param));
    new_props.reset(new BlackoilPropsAdFromDeck(deck, eclipseState, *grid));
    // check_well_controls = param.getDefault("check_well_controls", false);
    // max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
    // Rock compressibility.
    rock_comp.reset(new RockCompressibility(deck, eclipseState));

    // Gravity.
    gravity[2] = deck->hasKeyword("NOGRAV") ? 0.0 : unit::gravity;

    // Init state variables (saturation and pressure).
    if (param.has("init_saturation")) {
        initStateBasic(grid->numCells(), &(grid->globalCell())[0],
                       &(grid->logicalCartesianSize()[0]),
                       grid->numFaces(), UgGridHelpers::faceCells(*grid),
                       grid->beginFaceCentroids(),
                       grid->beginCellCentroids(), Dune::CpGrid::dimension,
                       *props, param, gravity[2], state);
        initBlackoilSurfvol(grid->numCells(), *props, state);
        enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour };
        if (pu.phase_used[Oil] && pu.phase_used[Gas]) {
            const int np = props->numPhases();
            const int nc = grid->numCells();
            for (int c = 0; c < nc; ++c) {
                state.gasoilratio()[c] = state.surfacevol()[c*np + pu.phase_pos[Gas]]
                    / state.surfacevol()[c*np + pu.phase_pos[Oil]];
            }
        }
    } else if (deck->hasKeyword("EQUIL") && props->numPhases() == 3) {
        state.init(grid->numCells(), grid->numFaces(), props->numPhases());
        const double grav = param.getDefault("gravity", unit::gravity);
        initStateEquil(*grid, *props, deck, eclipseState, grav, state);
        state.faceflux().resize(grid->numFaces(), 0.0);
    } else {
        initBlackoilStateFromDeck(grid->numCells(), &(grid->globalCell())[0],
                                  grid->numFaces(), UgGridHelpers::faceCells(*grid),
                                  grid->beginFaceCentroids(),
                                  grid->beginCellCentroids(), Dune::CpGrid::dimension,
                                  *props, deck, gravity[2], state);
    }


    // The capillary pressure is scaled in new_props to match the scaled capillary pressure in props.
    if (deck->hasKeyword("SWATINIT")) {
        const int nc = grid->numCells();
        std::vector<int> cells(nc);
        for (int c = 0; c < nc; ++c) { cells[c] = c; }
        std::vector<double> pc = state.saturation();
        props->capPress(nc, state.saturation().data(), cells.data(), pc.data(),NULL);
        new_props->setSwatInitScaling(state.saturation(),pc);
    }

    BlackoilState distributed_state;
    std::shared_ptr<Opm::BlackoilPropsAdFromDeck> distributed_props = new_props;
    Dune::CpGrid distributed_grid = *grid;
    // At this point all properties and state variables are correctly initialized
    // If there are more than one processors involved, we now repartition the grid
    // and initilialize new properties and states for it.
    bool must_distribute = ( grid->comm().size() > 1 );
    if( must_distribute )
    {
        if( param.getDefault("output_matlab", false) || param.getDefault("output_ecl", true) )
        {
            OPM_THROW(std::logic_error, "We only support vtk output during parallel runs. "
                      <<"Please use \"output_matlab=false output_ecl=false\" to deactivate the "
                      <<"other outputs!");
        }
        grid->loadBalance();
        Dune::CpGrid global_grid      = *grid;
        distributed_grid = *grid;
        global_grid.switchToGlobalView();
        distributed_grid.switchToDistributedView();
        distributed_props = std::make_shared<BlackoilPropsAdFromDeck>(*new_props, grid->numCells());
        distributed_state.init(grid->numCells(), grid->numFaces(), state.numPhases());
        // init does not resize surfacevol. Do it manually.
        distributed_state.surfacevol().resize(grid->numCells()*state.numPhases(),
                                             std::numeric_limits<double>::max());
        Opm::BlackoilStateDataHandle state_handle(global_grid, distributed_grid,
                                                  state, distributed_state);
        Opm::BlackoilPropsDataHandle props_handle(*new_props,
                                                  *distributed_props);
        grid->scatterData(state_handle);
        grid->scatterData(props_handle);
    }

    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    const double *grav = use_gravity ? &gravity[0] : 0;

    // Solver for Newton iterations.
    std::unique_ptr<NewtonIterationBlackoilInterface> fis_solver;
    
    boost::any parallel_information;
    Opm::extractParallelGridInformationToISTL(*grid, parallel_information);
    if (param.getDefault("use_cpr", true)) {
        fis_solver.reset(new NewtonIterationBlackoilCPR(param, parallel_information));
    } else {
        fis_solver.reset(new NewtonIterationBlackoilSimple(param, parallel_information));
    }

    Opm::TimeMapConstPtr timeMap(eclipseState->getSchedule()->getTimeMap());
    SimulatorTimer simtimer;

    // initialize variables
    simtimer.init(timeMap);

    Opm::DerivedGeology geology(distributed_grid, *distributed_props, eclipseState, false, grav);

    std::vector<double> threshold_pressures = thresholdPressures(eclipseState, distributed_grid);

    SimulatorFullyImplicitBlackoil<Dune::CpGrid> simulator(param,
                                                           distributed_grid,
                                                           geology,
                                                           *distributed_props,
                                                           rock_comp->isActive() ? rock_comp.get() : 0,
                                                           *fis_solver,
                                                           grav,
                                                           deck->hasKeyword("DISGAS"),
                                                           deck->hasKeyword("VAPOIL"),
                                                           eclipseState,
                                                           outputWriter,
                                                           threshold_pressures);

    if( grid->comm().rank()==0 )
    {
        std::cout << "\n\n================ Starting main simulation loop ===============\n"
                  << std::flush;
    }

    SimulatorReport fullReport = simulator.run(simtimer, must_distribute ? distributed_state : state);

    if( grid->comm().rank()==0 )
    {
        std::cout << "\n\n================    End of simulation     ===============\n\n";
        fullReport.report(std::cout);
    }

    if (output) {
        std::string filename = output_dir + "/walltime.txt";
        std::fstream tot_os(filename.c_str(),std::fstream::trunc | std::fstream::out);
        fullReport.reportParam(tot_os);
        warnIfUnusedParams(param);
    }
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

