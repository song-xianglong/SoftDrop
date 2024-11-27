# SoftDrop isolation on exploring QED splitting function.

## Project Overview

This project explores photon isolation in high-energy particle collisions using **Pythia** for event generation and **FastJet** for jet clustering and photon isolation. The study focuses on understanding the QED splitting function for the process $q \to q\gamma$, employing grooming techniques like the **Soft Drop Algorithm**. It also validates theoretical expectations by comparing observed momentum distributions with the QED splitting function.

---

## Directory Structure

- **`code/`**: Contains essential scripts for simulation and analysis.
  - `lundplane_all.cc`: Core file for generating events using Pythia and analyzing their parton shower properties.
  - `macro.C`: ROOT macro for visualizing and plotting results.
  - `compile_pythia.sh`: Shell script to compile and run the Pythia simulations.
- **`manuscript/`**: Includes the report describing the methodology, results, and conclusions of the study.

---

## Instructions for Running the Project

### Prerequisites

1. **Environment**:
   - Access to a Linux environment (e.g., CERN's lxplus).
   - ROOT and Pythia installed in the system.
2. **Dependencies**:
   - FastJet: Ensure version 3.4.1 or later is installed.
   - GCC/G++ for compiling the code.

### Compilation and Execution

1. Compile the simulation code using:
   ```bash
   ./compile_pythia.sh
This script compiles lundplane_all.cc and links it with the necessary libraries.

2. Run the simulation:
   ```bash
   ./lundplane_all
This will generate the event data required for the study.

3. Use `macro.C` for plotting and analyzing the results:
   ```bash
   root -l macro.C
## Additional Information

For a detailed explanation of the physics concepts, simulation parameters, and analysis techniques used in this project, please refer to the `manuscript/` directory.