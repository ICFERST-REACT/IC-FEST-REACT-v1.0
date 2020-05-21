
#include <stdio.h>
#include <stdlib.h>
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"

int worker_tasks_cc(int *task_number, void * cookie);
int do_something(void *cookie);

double my_basic_callback(double x1, double x2, const char *str, void *cookie);
void register_basic_callback(void *cookie);

void AdvectCpp(std::vector<double> &c, std::vector<double> bc_conc, int ncomps,
		int nxyz, int dim);
double my_basic_callback(double x1, double x2, const char *str, void *cookie);

class my_data {
public:
	PhreeqcRM *PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif
	std::vector<double> *hydraulic_K;
};

int main(void) {

	// --------------------------------------------------------------------------
	// Create PhreeqcRM
	// --------------------------------------------------------------------------

	int nxyz = 40;
	std::vector<double> hydraulic_K;
	for (int i = 0; i < nxyz; i++) {
		hydraulic_K.push_back(i * 2.0);
	}
	my_data some_data;
	some_data.hydraulic_K = &hydraulic_K;

	int nthreads = 3;
	PhreeqcRM phreeqc_rm(nxyz, nthreads);
	some_data.PhreeqcRM_ptr = &phreeqc_rm;
	int status = phreeqc_rm.SetUnitsSolution(2);

	// Open files
	status = phreeqc_rm.SetFilePrefix("Advect_cpp");
	phreeqc_rm.OpenFiles();

	// Set concentration units
	status = phreeqc_rm.SetUnitsSolution(2);     // 1, mg/L; 2, mol/L; 3, kg/kgs
	status = phreeqc_rm.SetUnitsPPassemblage(1); // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = phreeqc_rm.SetUnitsExchange(1); // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = phreeqc_rm.SetUnitsSurface(1); // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = phreeqc_rm.SetUnitsGasPhase(1); // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = phreeqc_rm.SetUnitsSSassemblage(1); // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	status = phreeqc_rm.SetUnitsKinetics(1); // 0, mol/L cell; 1, mol/L water; 2 mol/L rock

	// Set conversion from seconds to user units (days)
	double time_conversion = 1.0 / 86400;
	status = phreeqc_rm.SetTimeConversion(time_conversion);
	// Set representative volume
	std::vector<double> rv;
	rv.resize(nxyz, 1.0);
	status = phreeqc_rm.SetRepresentativeVolume(rv);
	// Set initial porosity
	std::vector<double> por;
	por.resize(nxyz, 0.2);
	status = phreeqc_rm.SetPorosity(por);
	// Set initial saturation
	std::vector<double> sat;
	sat.resize(nxyz, 1.0);
	status = phreeqc_rm.SetSaturation(sat);
	// Set cells to print chemistry when print chemistry is turned on
	std::vector<int> print_chemistry_mask;
	print_chemistry_mask.resize(nxyz, 0);
	for (int i = 0; i < nxyz / 2; i++) {
		print_chemistry_mask[i] = 1;
	}
	status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask);
	// test getters
	const std::vector<int> & print_chemistry_mask1 =
			phreeqc_rm.GetPrintChemistryMask();
	const std::vector<bool> & print_on = phreeqc_rm.GetPrintChemistryOn();
	bool rebalance = phreeqc_rm.GetRebalanceByCell();
	double f_rebalance = phreeqc_rm.GetRebalanceFraction();
	bool so_on = phreeqc_rm.GetSelectedOutputOn();
	int units_exchange = phreeqc_rm.GetUnitsExchange();
	int units_gas_phase = phreeqc_rm.GetUnitsGasPhase();
	int units_kinetics = phreeqc_rm.GetUnitsKinetics();
	int units_pp_assemblage = phreeqc_rm.GetUnitsPPassemblage();
	int units_solution = phreeqc_rm.GetUnitsSolution();
	int units_ss_exchange = phreeqc_rm.GetUnitsSSassemblage();
	int units_surface = phreeqc_rm.GetUnitsSurface();

	// Demonstation of mapping, two equivalent rows by symmetry
	std::vector<int> grid2chem;
	grid2chem.resize(nxyz, -1);
	for (int i = 0; i < nxyz / 2; i++) {
		grid2chem[i] = i;
		grid2chem[i + nxyz / 2] = i;
	}
	status = phreeqc_rm.CreateMapping(grid2chem);
	if (status < 0)
		phreeqc_rm.DecodeError(status);
	int nchem = phreeqc_rm.GetChemistryCellCount();

	// --------------------------------------------------------------------------
	// Set initial conditions
	// --------------------------------------------------------------------------
	// Set printing of chemistry file
	status = phreeqc_rm.SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
	// Load database
	status = phreeqc_rm.LoadDatabase("llnl.dat");

	// Demonstrate add to Basic: Set a function for Basic CALLBACK after LoadDatabase
//	register_basic_callback(&some_data);

	// Demonstration of error handling if ErrorHandlerMode is 0
	if (status != IRM_OK) {
		std::cerr << phreeqc_rm.GetErrorString(); // retrieve error messages if needed
		throw PhreeqcRMStop();
	}

	// Run file to define solutions and reactants for initial conditions, selected output
	bool workers = true; // Worker instances do the reaction calculations for transport
	bool initial_phreeqc = true; // InitialPhreeqc instance accumulates initial and boundary conditions
	bool utility = true;         // Utility instance is available for processing
	status = phreeqc_rm.RunFile(workers, initial_phreeqc, utility,
			"advect.pqi");
	// Clear contents of workers and utility
	initial_phreeqc = false;
	std::string input = "DELETE; -all";
	status = phreeqc_rm.RunString(workers, initial_phreeqc, utility,
			input.c_str());
	// Determine number of components to transport
	int ncomps = phreeqc_rm.FindComponents();
	// Print some of the reaction module information
	{
		std::ostringstream oss;
		oss << "Database:                                         "
				<< phreeqc_rm.GetDatabaseFileName().c_str() << "\n";
		oss << "Number of threads:                                "
				<< phreeqc_rm.GetThreadCount() << "\n";
		oss << "Number of MPI processes:                          "
				<< phreeqc_rm.GetMpiTasks() << "\n";
		oss << "MPI task number:                                  "
				<< phreeqc_rm.GetMpiMyself() << "\n";
		oss << "File prefix:                                      "
				<< phreeqc_rm.GetFilePrefix() << "\n";
		oss << "Number of grid cells in the user's model:         "
				<< phreeqc_rm.GetGridCellCount() << "\n";
		oss << "Number of chemistry cells in the reaction module: "
				<< phreeqc_rm.GetChemistryCellCount() << "\n";
		oss << "Number of components for transport:               "
				<< phreeqc_rm.GetComponentCount() << "\n";
		oss << "Partioning of UZ solids:                          "
				<< phreeqc_rm.GetPartitionUZSolids() << "\n";
		oss << "Error handler mode:                               "
				<< phreeqc_rm.GetErrorHandlerMode() << "\n";
		phreeqc_rm.OutputMessage(oss.str());
	}

	const std::vector<int> &f_map = phreeqc_rm.GetForwardMapping();
	// Get component information
	const std::vector<std::string> &components = phreeqc_rm.GetComponents();
	const std::vector<double> & gfw = phreeqc_rm.GetGfw();
	for (int i = 0; i < ncomps; i++) {
		std::ostringstream strm;
		strm.width(10);
		strm << components[i] << "    " << gfw[i] << "\n";
		phreeqc_rm.OutputMessage(strm.str());
	}
	phreeqc_rm.OutputMessage("\n");

	// Set array of initial conditions
	std::vector<int> ic1, ic2;
	ic1.resize(nxyz * 7, -1);
	ic2.resize(nxyz * 7, -1);
	std::vector<double> f1;
	f1.resize(nxyz * 7, 1.0);
	for (int i = 0; i < nxyz; i++) {
		ic1[i] = 1;              // Solution 1
		ic1[nxyz + i] = 1;      // Equilibrium phases none
		ic1[2 * nxyz + i] = -1;     // Exchange 1
		ic1[3 * nxyz + i] = -1;    // Surface none
		ic1[4 * nxyz + i] = -1;    // Gas phase none
		ic1[5 * nxyz + i] = -1;    // Solid solutions none
		ic1[6 * nxyz + i] = -1;    // Kinetics none
	}
	status = phreeqc_rm.InitialPhreeqc2Module(ic1, ic2, f1);
	// No mixing is defined, so the following is equivalent
	// status = phreeqc_rm.InitialPhreeqc2Module(ic1.data());

	// alternative for setting initial conditions
	// cell number in first argument (-1 indicates last solution, 40 in this case)
	// in advect.pqi and any reactants with the same number--
	// Equilibrium phases, exchange, surface, gas phase, solid solution, and (or) kinetics--
	// will be written to cells 18 and 19 (0 based)
	std::vector<int> module_cells;
	module_cells.push_back(18);
	module_cells.push_back(19);
	status = phreeqc_rm.InitialPhreeqcCell2Module(-1, module_cells);
	// Get temperatures
	const std::vector<double> & tempc = phreeqc_rm.GetTemperature();
	// get current saturation
	std::vector<double> current_sat;
	status = phreeqc_rm.GetSaturation(current_sat);
	// Initial equilibration of cells
	double time = 0.0;
	double time_step = 0.0;
	std::vector<double> c;
	c.resize(nxyz * components.size());
	status = phreeqc_rm.SetTime(time);
	status = phreeqc_rm.SetTimeStep(time_step);
	status = phreeqc_rm.RunCells();
	status = phreeqc_rm.GetConcentrations(c);

	// --------------------------------------------------------------------------
	// Set boundary condition
	// --------------------------------------------------------------------------

	std::vector<double> bc_conc, bc_f1;
	std::vector<int> bc1, bc2;
	int nbound = 1;
	bc1.resize(nbound, 1);          // solution 0 from Initial IPhreeqc instance
	bc2.resize(nbound, -1);                     // no bc2 solution for mixing
	bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
	status = phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, bc1, bc2, bc_f1);



	// --------------------------------------------------------------------------
	// Additional features and finalize
	// --------------------------------------------------------------------------

	// Use utility instance of PhreeqcRM to calculate pH of a mixture
	std::vector<double> c_well;
	c_well.resize(1 * ncomps, 0.0);
	for (int i = 0; i < ncomps; i++) {
		c_well[i] = 0.5 * c[0 + nxyz * i] + 0.5 * c[9 + nxyz * i];
	}
	std::vector<double> tc, p_atm;
	tc.resize(1, 15.0);
	p_atm.resize(1, 3.0);
	IPhreeqc * util_ptr = phreeqc_rm.Concentrations2Utility(c_well, tc, p_atm);
	input = "SELECTED_OUTPUT; -pH;RUN_CELLS; -cells 1";
	int iphreeqc_result;
//	util_ptr->SetOutputFileName("utility_cpp.txt");
	util_ptr->SetOutputFileOn(true);
	iphreeqc_result = util_ptr->RunString(input.c_str());


return EXIT_SUCCESS;
}
