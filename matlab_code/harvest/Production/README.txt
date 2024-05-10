The scripts in this directory were developed by:
Eric Frizzell (efrizz@umd.edu) (January 2020 - May 2024)
University of Maryland, Aerospace Engineering, 3178 Glenn L. Martin Hall
College Park, MD, 20740, United States
https://github.com/efrizz-umd/

The MATLAB scripts in this folder (the code) process LIGGGHTS output files (*.post). Running the code
results in a *.mat file which contains time history data of a granular channel that has
been subjected to a piston impact. Particle states are averaged into yz bins. We measure two
quantities related to packing fraction (bulk density) changes: 1) channel height and 2) local
packing fraction. The resultant *.mat file is then used by our data analysis scripts, found at our public
repository. This code is provided 'as-is', without any guarantees or warranties of any kind, express or
implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, or
non-infringement. The author shall not be held responsible or liable for any direct, indirect, incidental,
special, exemplary, or consequential damages arising out of the use or inability to use the code. By
using this code, you agree that the author assumes no responsibility for providing technical support
or troubleshooting assistance.

LIGGGHTS:
https://www.cfdem.com/liggghts-open-source-discrete-element-method-particle-simulation-code


This code was used in the following works:
1) Frizzell and Hartzell 2023, Simulation of lateral impulse induced inertial
 dilation at the surface of a vacuum‑exposed granular assembly
https://doi.org/10.1007/s10035-023-01363-6
2) Frizzell and Hartzell 2024, Material parameter influence on the expression of
 Solitary-Wave-Induced Surface Dilation (submitted, Granular Matter)
3) Frizzell and Hartzell 2024, Wave transmission through the megaregolith as a
 mechanism for lunar cold spot formation (submitted, Icarus)

% ----- Workflow overview ----- %
1) relocate the LIGGGHTS post file directories (both dump*.post and dump_computes*.post)
to the same directory as these scripts
2) ensure that the 'mat_files' directory exists
3) update harvest_script_decay.m to match output rate and scaling factor details
4) launch the job script with appropriate fields altered (the job script launches the wrapper script,
harvest_script_decay.m, which manages the process and generates the *.mat)


% ----- Function overview ----- %
# --------------------------------------------- #
harvest_script_decay.m
- wrapper script that manages the data binning and accumulation process
- generates a *.mat file, to be used in our analysis code
- variables to be updated before use: file name, output time step, sim  timestep, scaling factor


# --------------------------------------------- #
HarvestNoPlot_Production.m
- initializes storage for the various bins
- loops through the data in time order and manages updating and storage
- generates the final *.mat file


# --------------------------------------------- #
volumeFraction_updateFinderChunk_overlaps.m
- opens and reads the data
- passes the data to the appropriate binning functions

# --------------------------------------------- #
overlap_Finder_PID.m
- takes the file filled with all particle overlaps and computes the total and averaged
overlap per particle
- uses the particle ID as a lookup for data storage

# --------------------------------------------- #
force_updateChunkFinder_overlaps.m
- establishes grid bounds if not already established
- selects sensor particles in as the closest particle to the center of a bins or uses
the passed in ID to look up the sensor particle
- updates sensor particles state information
- manages the binning scripts (chunkParticles_overlaps.m for all state information,
chunkParticles_voro for the local packing fraction)
- computes the current height of the channel


# --------------------------------------------- #
chunkParticles_overlaps.m
- averages particle data into yz bins

# --------------------------------------------- #
chunkParticles_voro.m
- averages voronoi volume into yz bins

# --------------------------------------------- #
voronize.m
- performs a voronoi tessellation of each particle in the system
- output contains each particle's voronoi volume (used to compute its local packing fraction)
- output also contains corresponding particle positions to facilitate chunking

# --------------------------------------------- #
volsphere.m
-computes volume of a sphere


% ----- RECOMMENDED UPDATES ----- %
1. Presently, this code is designed to be run in serial and can take several hours depending
on the file size (around 250 MB for our tests) and number of files (typically 200-300). However, the code
can (and should) be parallelized which could be accomplished with MATLAB's par4
parallelization of HarvestNoPlot_Production.m.

2. The use of vector arrays no longer makes sense, structures or classes should be implemented
