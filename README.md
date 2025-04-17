# PSA Process Optimization Using NSGA-II

This project was part of an exploratory study focused on optimising a Pressure Swing Adsorption (PSA) process for hydrogen separation. The MATLAB script uses the NSGA-II multi-objective evolutionary algorithm to optimise key process variables with respect to purity and recovery objectives.

---

## Background

Pressure Swing Adsorption (PSA) is a cyclic process used to separate gas mixtures based on molecular characteristics and affinity for an adsorbent. This simulation investigates how changes in adsorption time, flow rates, and pressure affect the hydrogen purity and recovery of the system.

---

## Methodology

- **Optimisation Engine**: NSGA-II (Non-dominated Sorting Genetic Algorithm II)
- **Simulation Framework**: 1D dynamic column model for PSA
- **Objective Functions**: Maximise hydrogen purity and recovery
- **Decision Variables**:
  - Adsorption time
  - Feed and purge flow velocities
  - Adsorption pressure

---

## Key Features

- Fully parallelised optimisation using MATLAB’s `parpool` for faster evaluation
- Integration of NSGA-II with a detailed physical PSA model
- Post-processing of Pareto-optimal solutions and visualisation of trade-offs
- Plotting of the Pareto front for purity vs. recovery
- Visualisation of feature influence across rank-1 individuals

---

## Structure

- `PSA_Process.m`: Core simulation model used as the objective function
- `NSGA-II/`: Contains the optimisation engine and dependencies
- `main.m`: Runs NSGA-II with bounds, objectives, and output logic
- Output: Pareto front showing trade-off between hydrogen purity and recovery

---

## Output Example

The model outputs a Pareto front of optimal trade-offs between hydrogen purity and hydrogen recovery. Feature values for each rank-1 individual are stored for deeper analysis and plotting.

---

## Tools Used

- MATLAB (with Parallel Computing Toolbox)
- NSGA-II Toolbox (MATLAB-based)
- Custom PSA solvers for mass, energy, and momentum balance

---

## Author

**Thomas Stone-Wigg**  
Focused on energy systems modeling, sustainability, and applied numerical optimization.
