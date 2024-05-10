Wave transmission through the megaregolith as a mechanism for lunar cold spot formation
Eric Frizzell (efrizz@umd.edu) and Christine Hartzell
University of Maryland, Aerospace Engineering, 3178 Glenn L. Martin Hall
College Park, MD, 20740, United States
Manuscript submitted for publication May, 2024
https://github.com/efrizz-umd/SID_sensitivity

The text files in these folders are LIGGGHTS input scripts used in Frizzell & Hartzell, 2024
Wave transmission through the megaregolith as a mechanism for lunar cold spot formation (submitted).

The main folder, 'input_scripts', contains scripts for initiating impacts within the prepared channels ('decay scripts'). The restart files are then
used as in input to the decay script, allowing multiple parameters to be studied without preforming channel preparation for each test.
The restart files we produced in this work are located at the link below. Many of the restart files used were generated in our prior work and can be found at:
https://zenodo.org/doi/10.5281/zenodo.11003039. We used the channel fill works from that same prior work to generate the boulder scale particles
in the present work, and those restart files can be found at the link below.

### Restart files ###
https://doi.org/10.5281/zenodo.11176402

Some cases were run on our stand-alone lab server (44 core). After compiling LIGGGHTS, the files are run as:
mpirun -np NUM_PROC_HERE ./lmp_auto -echo both -i in.decay

### Piston Impacts ###
Although we have run tests for each parameter at several different values, we typically just include one 'in.decay' file per test type
(ex: in.decay). The user would need to make sure that the input file parameter of interest (in our example, the LIGGGHTS variable
'coefficientFriction' would be updated) is set to the value desired. There are two methods of using restart files with the input scripts. First,
is to use a restart file that was prepared with the exact same material parameters as in the decay script. The second method involves changing a
parameter from what is already represented in the restart file (ex: restart file was prepared with coefficientFriction = 0.6 but you now want to
impact a channel which has particles of coefficientFriction = 1.2). In this case you will want to make sure that the few lines of code before the
impact occurs (where the piston velocity is set) are uncommented to allow the particles to 'relax'. In some cases, changing particle property
requires quite a bit of settling time (like increasing elastic modulus by orders of magnitude).

These cases were run on the UMD HPC cluster Zaratan (128 cores/node) and the included 'job_decay.sub' provides an example of running the piston impact scripts.
The job script is written to be run from the same folder as in the input script.


### References ###
1 - Frizzell, E.S., Hartzell, C.M. Simulation of lateral impulse induced inertial dilation at the surface of a vacuum-exposed granular assembly.
Granular Matter 25, 75 (2023). https://doi.org/10.1007/s10035-023-01363-6

LIGGGHTS: https://www.cfdem.com/liggghtsr-open-source-discrete-element-method-particle-simulation-code
