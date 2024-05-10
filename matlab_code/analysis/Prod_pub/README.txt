The scripts in this directory (the code) were developed by:
Eric Frizzell (efrizz@umd.edu) (January 2020 - May 2024)
University of Maryland, Aerospace Engineering, 3178 Glenn L. Martin Hall
College Park, MD, 20740, United States
https://github.com/efrizz-umd/

% ----- Associated works ----- %
Th code was used in the following works:
1) Frizzell and Hartzell 2023, Simulation of lateral impulse induced inertial
 dilation at the surface of a vacuum‑exposed granular assembly
https://doi.org/10.1007/s10035-023-01363-6
2) Frizzell and Hartzell 2024, Material parameter influence on the expression of
 Solitary-Wave-Induced Surface Dilation (submitted, Granular Matter)
3) Frizzell and Hartzell 2024, Wave transmission through the megaregolith as a
 mechanism for lunar cold spot formation (submitted, Icarus)


% ----- Code description ----- %
The code performs analysis (plotting, wave tracking) on a *.mat file which is the output
of our harvesting function (found at https://github.com/efrizz-umd/runout). The harvested
*.mat file represents the window-averaged (binned) time history of particles within a granular
channel that experienced the lateral impulse of a piston. This code allows for the comparison
across several different test runs, which are generally grouped by varying analysis type (i.e.,
'boulder_velocity' represents a velocity sweep of the piston in boulder scale particles). Running
the code results in various plots that appeared in works 2 and 3 above (described below). Along with
this work we have included the *.mat files that resulted from three of the different comparison cases
considered in work 3 (runout sensitivity to particle size, piston velocity and piston mass). See the
Supplementary files section below for a link to the files. The *.mat files will need to be downloaded
and placed into a sub-directory within the 'mat_files' directory in order to run this code.

% ----- Code output ----- %
Presently, this code will output plots below. The functionality for many more diagnostic
type plots is included, but currently commented out. Assuming you have cleared your workspace
and closed all plots before running the code, then, for the varied cases you select, the plots
show each cases:

Figure 1 - h_init_VFvsdepth
normalized (by particle radius) depth vs packing fraction in the channel (used
to compute average and std. dev. packing fraction in all three works above)

Figure 2 - h_max_force
average particle force in the wave front vs normalized (by radius) radial position
(Fig. 11 in work 2)

Figure 3 - h_bulkdilationbinvsdepth
normalized (by particle radius) depth vs average percent bulk volume change
(Fig. 1, 3 and 9 in work 2)

Figure 4 - h_max_overlap
normalized (by initial value) average maximum overlap vs normalized (by particle
radius) distance traveled
(Fig. 3 in work 3)

Figure 5 - h_max_overlap_vstime
normalized (by initial value) average maximum overlap vs time
(Fig. 4 in work 3)

Figure 6 - h_max_velocity_new
average maximum particle velocity vs normalized (by particle radius) distance traveled
(Fig. 5 in work 3)

Figure 7 - h_max_overlap_vsdepth
normalized depth (by particle radius) vs average initial and maximum h_max_overlap in
a driven solitary wave
(Fig. 6 in work 2)

Figure 8 - h_max_velocity_new_vsdepth
depth vs average maximum particle velocity at a selected location (with a driven Solitary Wave)
(average over all depths in this Figure produces the depth at a single radial distance
in Figure 4 or a single time in Figure 3)

Figure 9 - h_max_force_new_vsdepth
normalized depth (by particle radius) vs average simulated force and predicted particle force
(Fig. 5 in work 2)


% ----- Supplementary files ----- %
The *.mat files can be found at:
https://doi.org/10.5281/zenodo.11176402

% ----- Workflow overview ----- %
1) Place the *.mat files to compare in the mat_files/TEST_TYPE_FOLDER_NAME
2) Ensure that the 'ProductionPlots/TEST_TYPE_FOLDER_NAME' directory exists
4) Update the header information in comp_script_load_only_production_nd_decay.m to
point to the desired files and describe their quality (set the piston velocity used,
create meaningful legend entries, etc.)
4) Run comp_script_load_only_production_nd_decay.m


% ----- Function overview ----- %
# --------------------------------------------- #
comp_script_load_only_production_nd_decay.m
- Wrapper script that manages the process and performs plotting


# --------------------------------------------- #
HarvestAndPlotFun_production_loadonly.m
- Takes an individual *.mat file and finds the averaged quantities


# --------------------------------------------- #
findFinalConditions_vf.m
findInitialConditions_vf.m
findFullAverages_xy.m
findFullaverages.m
findFullMax.m
findInitialConditions.m
findWindowAverages.m
- Performs averaging on an input storage array. There is some redundancy here and
some no longer used functionality that has not yet been removed.

# --------------------------------------------- #
FSW_impulse_fun_v4.m
- computes predicted force in the wave front

# --------------------------------------------- #
finite_diff_d1_fun.m
- Calculates derivative using finite difference scheme

# --------------------------------------------- #
init_overlap_calc_fun.m
- predicted initial overlap vs depth in the channel

# --------------------------------------------- #
least_fit_fun_nearestneighbor.m
- performs fit to data using least squares, then gives the initial slope (used
for finding wave speeds)

# --------------------------------------------- #
printer_cell_fun.m
printerfun.m
- prints data from arrays or cells

# ---------------------------------------------
star_calc.m
- computes equivalent material parameters during a contact

# ---------------------------------------------
terminationfun.m
- assess a termination point, given some criteria (not currently used)


% ----- Disclaimer ----- %
The code is presented as is and is intended only to demonstrate our procedure for generating
those plots listed above. There are many lines of experimental and work in progress code that
has been commented in order to present a minimal version, producing only the relevant plots. However,
an interested user could easily uncomment some of these features in order to add functionality.

This code is provided 'as-is', without any guarantees or warranties of any kind, express or
implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, or
non-infringement. The author shall not be held responsible or liable for any direct, indirect, incidental,
special, exemplary, or consequential damages arising out of the use or inability to use the code. By
using this code, you agree that the author assumes no responsibility for providing technical support
or troubleshooting assistance.
