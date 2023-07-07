#include <algorithm>  // For std::is_sorted
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"


static void declare_parameters(po::options_description &cfg,
                               Param &p)
{
    /* To have a new input parameter declared as such in parameters.hpp,
     *
     *     struct SectionType {
     *         type name;
     *     }
     *     struct Param {
     *         SectionType section;
     *     }
     *
     * add this line in this function:
     *
     *     ("section.name", po::value<type>(&p.section.name), "help string")
     *
     */

    cfg.add_options()
        ("sim.modelname", po::value<std::string>(&p.sim.modelname),
         "Prefix for the output files")

        ("sim.max_steps", po::value<int>(&p.sim.max_steps),
         "Max. number of time steps")
        ("sim.max_time_in_yr", po::value<double>(&p.sim.max_time_in_yr),
         "Max. time (in years)")
        ("sim.output_step_interval", po::value<int>(&p.sim.output_step_interval),
         "Output step interval")
        ("sim.output_time_interval_in_yr", po::value<double>(&p.sim.output_time_interval_in_yr),
         "Output time interval (in years)")

        ("sim.checkpoint_frame_interval", po::value<int>(&p.sim.checkpoint_frame_interval)->default_value(10),
         "How frequent to write checkpoint file (used for restarting simulation)?")
        ("sim.restarting_from_modelname", po::value<std::string>(&p.sim.restarting_from_modelname),
         "Prefix for the checkpoint and save files for restarting.")
        ("sim.restarting_from_frame", po::value<int>(&p.sim.restarting_from_frame),
         "Which output frame to read for restarting?")
        ("sim.is_restarting", po::value<bool>(&p.sim.is_restarting)->default_value(false),
         "Restarting from previous checkpoint file?")

        ("sim.has_initial_checkpoint", po::value<bool>(&p.sim.has_initial_checkpoint)->default_value(true),
         "Output checkpoint file at 0th step?")
        ("sim.has_marker_output", po::value<bool>(&p.sim.has_marker_output)->default_value(true),
         "Output marker coordinate and material?")
        ("sim.has_output_during_remeshing", po::value<bool>(&p.sim.has_output_during_remeshing)->default_value(false),
         "Output immediately before and after remeshing?")
        ("sim.is_outputting_averaged_fields", po::value<bool>(&p.sim.is_outputting_averaged_fields)->default_value(true),
         "Output time-averaged (smoothed) field variables or not. These fields are: velocity, strain rate, and stress.\n"
         "no: output instantaneous fields. The velocity and strain-rate might oscillate temporally.\n"
         "yes: output field variables averaged over mesh.quality_check_step_interval time steps.\n")
        ;

    cfg.add_options()
        ("mesh.meshing_option", po::value<int>(&p.mesh.meshing_option)->default_value(1),
         "How to create the new mesh?\n"
         "1: rectangular box with roughly uniform resolution\n"
         "2: rectangular box with rectangular zone of refined resolution\n"
         "90: read the bounding polygon from a .poly file\n"
         "91: same as 90, but the max. element size is normalized with resolution.\n"
         "95: read the mesh from from a .exo file\n"
         )
        ("mesh.meshing_verbosity", po::value<int>(&p.mesh.meshing_verbosity)->default_value(-1),
         "Output verbose during mesh/remeshing. -1 for no output.")
        ("mesh.meshing_sediment", po::value<bool>(&p.mesh.meshing_sediment)->default_value(false),
         "Refine the elements with sediment marker during remeshing.")
        ("mesh.tetgen_optlevel", po::value<int>(&p.mesh.tetgen_optlevel)->default_value(3),
         "Optimization level for tetgen. 0: no optimization; 1: multiple edge filps; 2: 1 & free vertex deletion; 3: 2 & new vertex insertion. High optimization level could slow down the speed of mesh generation. For 3D only.")

        ("mesh.xlength", po::value<double>(&p.mesh.xlength)->required(),
         "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(&p.mesh.ylength)->required(),
         "Length of y (in meters), for 3D only")
        ("mesh.zlength", po::value<double>(&p.mesh.zlength)->required(),
         "Length of z (in meters)")

        ("mesh.resolution", po::value<double>(&p.mesh.resolution)->required(),
         "Spatial resolution (in meters)")
        ("mesh.smallest_size", po::value<double>(&p.mesh.smallest_size)->default_value(0.01),
         "The size of smallest element relative to an element with typical resolution")
        ("mesh.largest_size", po::value<double>(&p.mesh.largest_size)->default_value(30),
         "The size of largest element relative to an element with typical resolution")
        ("mesh.sediment_size", po::value<double>(&p.mesh.sediment_size)->default_value(1.),
         "The size threshold of the refined sediment elements")
        // for 2D only
        ("mesh.min_angle", po::value<double>(&p.mesh.min_angle)->default_value(32.),
         "Min. angle of all triangles (in degrees), for 2D only")
        // for 3D only
        ("mesh.min_tet_angle", po::value<double>(&p.mesh.min_tet_angle)->default_value(22.),
         "Min. dihedral angle of all tetrahedra (in degrees), for 3D only")
        ("mesh.max_ratio", po::value<double>(&p.mesh.max_ratio)->default_value(2.),
         "Max. radius / length ratio of all tetrahedra, for 3D only")

        // for meshing_option = 2 only
        /* read the value as string, then parse as two numbers later */
        ("mesh.refined_zonex", po::value<std::string>()->default_value("[0.4, 0,6]"),
         "Refining portion of xlength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only")
        ("mesh.refined_zoney", po::value<std::string>()->default_value("[0.4, 0.6]"),
         "Refining portion of ylength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only, for 3D only")
        ("mesh.refined_zonez", po::value<std::string>()->default_value("[0.8, 1]"),
         "Refining portion of zlength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only")

        // for meshing_option = 90 only
        ("mesh.poly_filename", po::value<std::string>(&p.mesh.poly_filename)->default_value("mesh.poly"),
         "Filename of the input polygon, for meshing_option=90 or 91 only.\n"
         "The file format is described in\n"
         "http://www.cs.cmu.edu/~quake/triangle.poly.html (2D) and\n"
         "http://wias-berlin.de/software/tetgen/fformats.poly.html (3D).\n"
         "Limitation: no point attributes, no boundary marker for points, no holes, "
         "and no regional attributes.\n"
         "Users still need to provide these parameters: mesh.xlength, mesh.ylength, and mesh.zlength.")

        // for meshing_option = 95 only
        ("mesh.exo_filename", po::value<std::string>(&p.mesh.exo_filename)->default_value("mesh.exo"),
         "Filename of the input polygon, for meshing_option=95 only.\n"
         "Meshing software like Trelis is supposed to be used for creating a mesh in the ExodusII format.\n"
         "Limitation: Only box-like domains with six sides are supported.\n"
         "Users still need to provide these parameters: mesh.xlength, mesh.ylength, and mesh.zlength.")

        ("mesh.quality_check_step_interval", po::value<int>(&p.mesh.quality_check_step_interval)->default_value(100),
         "How often to check mesh quality?")
        ("mesh.min_quality", po::value<double>(&p.mesh.min_quality)->default_value(0.4),
         "Min. mesh quality before remeshing (between 0 and 1)")
        ("mesh.max_boundary_distortion", po::value<double>(&p.mesh.max_boundary_distortion)->default_value(0.25),
         "Max. distance of boundary distortion before remeshing (in unit of mesh.resolution).")

        ("mesh.remeshing_option", po::value<int>(&p.mesh.remeshing_option)->default_value(0),
         "How to deal with the boundaries during remeshing?\n"
         " 0: no modification on any boundary.\n"
         " 1: move all bottom nodes to initial depth, other boundaries are intact.\n"
         " 2: create a new bottom boundary at the initial depth, other boundaries are intact (2D only).\n"
         "10: no modification on any boundary, except small boundary segments might get merged.\n"
         "11: move all bottom nodes to initial depth, other boundaries are intact, small boundary segments might get merged.\n")

        ("mesh.is_discarding_internal_segments", po::value<bool>(&p.mesh.is_discarding_internal_segments)->default_value(true),
         "Discarding internal segments after initial mesh is created? "
         "Using it when remeshing process can modify segments (e.g. remeshing_option=11).")

        // for mesh optimization with MMG2D/3D
        ("mesh.mmg_debug", po::value<int>(&p.mesh.mmg_debug)->default_value(0),
         "Run MMG remesher in debug mode? No:0; Yes:1\n")
        ("mesh.mmg_verbose", po::value<int>(&p.mesh.mmg_verbose)->default_value(0),
         "Verbosity level of MMG remesher. For debugging, set a value greater than 4.\n")
        ("mesh.mmg_hmax_factor", po::value<double>(&p.mesh.mmg_hmax_factor)->default_value(2.0),
         "Factor multiplied to param.mesh.resolution to set the maximum element size\n")
         ("mesh.mmg_hmin_factor", po::value<double>(&p.mesh.mmg_hmin_factor)->default_value(0.2),
         "Factor multiplied to param.mesh.resolution to set the minimum element size\n")
         ("mesh.mmg_hausd_factor", po::value<double>(&p.mesh.mmg_hausd_factor)->default_value(0.01),
         "Factor multiplied to param.mesh.resolution to set the Hausdorff distance between original and remeshed surfaces.\n")
        ;

    cfg.add_options()
        ("markers.init_marker_option", po::value<int>(&p.markers.init_marker_option)->default_value(1),
         "How to generate markers?\n"
         "1: randomly distributed markers.\n"
         "2: regularly spaced markers.\n")
        ("markers.markers_per_element", po::value<int>(&p.markers.markers_per_element)->default_value(4),
         "Number of markers per element. Used when init_marker_option=1.")
        ("markers.init_marker_spacing", po::value<double>(&p.markers.init_marker_spacing)->default_value(0.3),
         "Spacing of markers (in unit of mesh.resolution). Used when init_marker_option=2.")
        ("markers.min_num_markers_in_element", po::value<int>(&p.markers.min_num_markers_in_element)->default_value(3),
         "When the number of markers in an element is less than this number, a new marker will be replenished in the element.")
        ("markers.replenishment_option", po::value<int>(&p.markers.replenishment_option)->default_value(2),
         "How to determine the mattype of replenished markers?\n"
         "0: always set to 0 (fastest option).\n"
         "1: by the probability of marker mattype of the element or surrounding elements.\n"
         "2: same as the mattype of the nearest marker (slowest option).")
        ("markers.random_seed", po::value<uint>(&p.markers.random_seed)->default_value(1),
         "Random seed of marker position. If 0, the current time is used as the seed.")
        ;

    cfg.add_options()
        ("control.gravity", po::value<double>(&p.control.gravity)->default_value(10),
         "Magnitude of the gravity (in m/s^2)")

        ("control.characteristic_speed",
         po::value<double>(&p.control.characteristic_speed)->default_value(0),
         "Characteristic tectonic speed (in m/s). "
         "It is used to estimate the size of stable time step. "
         "If it is 0, its value is inferred from the imposed boundary velocity. "
         "Set it to other value if the imposed boundary velocity is 0 everywhere.")

        ("control.is_quasi_static", po::value<bool>(&p.control.is_quasi_static)->default_value(true),
         "Is the simulation quasi-static or dynamic? If quasi-static, inertial scaling and strong damping is applied.\n")
        ("control.dt_fraction", po::value<double>(&p.control.dt_fraction)->default_value(1.0),
         "Take dt as a fraction of max. stable time step size (0-1).\n")
        ("control.fixed_dt", po::value<double>(&p.control.fixed_dt)->default_value(0),
         "Fixed dt size (in seconds). If 0, dt sized will be determined dynamically.\n")
        ("control.inertial_scaling", po::value<double>(&p.control.inertial_scaling)->default_value(1e5),
         "Scaling factor for inertial (a large number)")

        ("control.damping_option", po::value<int>(&p.control.damping_option)->default_value(1),
         "How to damp the elastic wave?\n"
         "0: no damping.\n"
         "1: damping/acceleration depends on force and direction of velocity.\n"
         "2: damping depends on force.\n"
         "3: damping/(weaker) acceleration depends on force and direction of velocity.\n"
         "4: Rayleigh damping.\n")
        ("control.damping_factor", po::value<double>(&p.control.damping_factor)->default_value(0.8),
         "A factor for force damping (0-1)")

        ("control.ref_pressure_option", po::value<int>(&p.control.ref_pressure_option)->default_value(0),
         "How to define reference pressure?\n"
         "0: using density of the 0-th element to compute lithostatic pressure.\n"
         "1: computing reference pressure from the PREM model.\n"
         "2: computing reference pressure from the PREM model, modified for continent.\n")
//        ("control.surface_pressure_correction", po::value<bool>(&p.control.surface_pressure_correction)->default_value(false),
//         "Correct the pressure of surface elements"
//         "which has positive stress 1st invariant"
//         "and force the 1st invariant to zero.")
        ("control.is_using_mixed_stress", po::value<bool>(&p.control.is_using_mixed_stress)->default_value(true),
         "If use Nodal Mixed Discretization For Stress")
        ("control.mixed_stress_reference_viscosity", po::value<double>(&p.control.mixed_stress_reference_viscosity)->default_value(1.e19),
         "The reference viscosity for appling mixed stress.")

        ("control.surface_process_option", po::value<int>(&p.control.surface_process_option)->default_value(0),
         "What kind of surface processes? 0: no surface processes. "
         "1: using simple diffusion to modify surface topography. "
         "101: custom function."
         "2: using simple deposition."
         "102: using simple deposition and diffusion.")
        ("control.surface_diffusivity", po::value<double>(&p.control.surface_diffusivity)->default_value(1e-6),
         "Standard diffusition coefficient of surface topography (m^2/s)")
        ("control.surf_diff_ratio_terrig", po::value<double>(&p.control.surf_diff_ratio_terrig)->default_value(1.),
         "Ratio between terrigenous and standard Diffusition coefficients.")
        ("control.surf_diff_ratio_marine", po::value<double>(&p.control.surf_diff_ratio_marine)->default_value(1.),
         "Ratio between marine and standard Diffusition coefficients.")
        ("control.surf_depo_universal", po::value<double>(&p.control.surf_depo_universal)->default_value(0.),
         "Sedimentation rate of universal suspended deposite (in m/s).")
        ("control.surf_base_level",po::value<double>(&p.control.surf_base_level)->default_value(0.e0),
         "Base level of surface processes.")
        ("control.terrig_sediment_diffusivity",po::value<double>(&p.control.terrig_sediment_diffusivity)->default_value(3.17e-6),
         "Submarine diffusion coefficient (=100 m^2/yr, Kaufman et al., 1991) (in m^2/s).")
        ("control.terrig_depth_coefficient",po::value<double>(&p.control.terrig_depth_coefficient)->default_value(5e-4),
         "Submarine diffusion coefficient decay (Kaufman et al., 1991) (in m^-1).")
        ("control.is_reporting_terrigenous_info",po::value<bool>(&p.control.is_reporting_terrigenous_info)->default_value(false),
         "If report surface terrigenous information.")
        ("control.hemipelagic_sedimentation_rate",po::value<double>(&p.control.hemipelagic_sedimentation_rate)->default_value(3.17e-15),
         "At 1km depth, 1e-1 mm/yr ~= 3.17e-12 m/s / 1e3 m = 3.17e-15 m/s")

        ("control.terrig_sediment_volume",po::value<double>(&p.control.terrig_sediment_volume)->default_value(1.),
         "The volume of the sediment bring from continent source (in m^3/s). "
         "5e7 ton/yr ~= 1 m^3/s")
        ("control.terrig_sediment_area",po::value<double>(&p.control.terrig_sediment_area)->default_value(3.17e-7),
         "The are of the sediment bring from continent source (=10 m^2/yr) (in m^2/s). "
         "Assuming the width of the channel is 50 km, 5e7 ton/yr ~= 2.e-5 m^2/s")

        ("control.has_thermal_diffusion", po::value<bool>(&p.control.has_thermal_diffusion)->default_value(true),
         "Does the model have thermal diffusion? If not, temperature is advected, but not diffused.\n")

        ("control.has_hydration_processes", po::value<bool>(&p.control.has_hydration_processes)->default_value(false),
         "Does the model have hydration processes? It is required to model some types of phase changes.")
        ("control.hydration_migration_speed", po::value<double>(&p.control.hydration_migration_speed)->default_value(3e-9),
         "The upward migration speed of hydrous fluid (in m/s).\n")
        ;

    cfg.add_options()
        ("bc.surface_temperature", po::value<double>(&p.bc.surface_temperature)->default_value(273),
         "Surface temperature (in Kelvin)")
        ("bc.mantle_temperature", po::value<double>(&p.bc.mantle_temperature)->default_value(1600),
         "Mantle temperature (in Kelvin)")
        ("bc.has_winkler_foundation", po::value<bool>(&p.bc.has_winkler_foundation)->default_value(true),
         "Using Winkler foundation for the bottom boundary?")
        ("bc.winkler_delta_rho", po::value<double>(&p.bc.winkler_delta_rho)->default_value(0),
         "Excess density of the bottom Winkler foundation (in kg/m^3)")

        ("bc.has_elastic_foundation", po::value<bool>(&p.bc.has_elastic_foundation)->default_value(false),
         "Using elastic foundation for the bottom boundary?")
        ("bc.elastic_foundation_constant", po::value<double>(&p.bc.elastic_foundation_constant)->default_value(1e11),
         "Elastic constant for elastic foundation.")

        ("bc.has_water_loading", po::value<bool>(&p.bc.has_water_loading)->default_value(true),
         "Applying water loading for top boundary that is below sea level?")

        ("bc.vbc_x0", po::value<int>(&p.bc.vbc_x0)->default_value(1),
         "Type of velocity boundary condition for the left/western side. "
         "Odd number indicates the normal component of the velocity is fixed. "
         "Possible type is \n"
         "0: all velocity components free;\n"
         "1: normal component fixed, shear components free;\n"
         "2: normal component free, shear components fixed at 0;\n"
         "3: normal component fixed, shear components fixed at 0;\n"
         "4: normal component free, shear component (not z) fixed, z component fixed at 0, only in 3D;\n"
         "5: normal component fixed at 0, shear component (not z) fixed, z component fixed at 0, only in 3D;\n"
         "7: normal component fixed, shear component (not z) fixed at 0, z component free, only in 3D;\n"
         "11: horizontal normal component fixed, z component and horizontal shear component free, only in 3D;\n"
         "13: horizontal normal component fixed, z component and horizontal shear component fixed at 0, only in 3D;\n")
        ("bc.vbc_x1", po::value<int>(&p.bc.vbc_x1)->default_value(1),
         "Type of boundary condition for the right/eastern side")
        ("bc.vbc_val_x0", po::value<double>(&p.bc.vbc_val_x0)->default_value(-1e-9),
         "Value of boundary condition for left/western side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_division_x0_min", po::value<double>(&p.bc.vbc_val_division_x0_min)->default_value(1.),
         "Ratio of first rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_division_x0_max", po::value<double>(&p.bc.vbc_val_division_x0_max)->default_value(1.),
         "Ratio of second rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_x0_ratio0", po::value<double>(&p.bc.vbc_val_x0_ratio0)->default_value(1.),
         "Ratio of the velocity at the surface of the left boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x0_ratio1", po::value<double>(&p.bc.vbc_val_x0_ratio1)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x0_min of the left boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x0_ratio2", po::value<double>(&p.bc.vbc_val_x0_ratio2)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x1_max of the left boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x0_ratio3", po::value<double>(&p.bc.vbc_val_x0_ratio3)->default_value(1.),
         "Ratio of the velocity at the bottom of the left boundary width (for option = 1 and 3)")

        ("bc.num_vbc_period_x0", po::value<int>(&p.bc.num_vbc_period_x0)->default_value(1),
         "Number of velocity period at x0")
        ("bc.num_vbc_period_x1", po::value<int>(&p.bc.num_vbc_period_x1)->default_value(1),
         "Number of velocity period at x1")
        ("bc.vbc_period_x0_time_in_yr", po::value<std::string>()->default_value("[1.e6]"),
         "End time of BC velocity periods at x0 side '[t0, t1, t2, ...]' (in yr)")
        ("bc.vbc_period_x1_time_in_yr", po::value<std::string>()->default_value("[1.e6]"),
         "End time of BC velocity periods at x1 side '[t0, t1, t2, ...]' (in yr)")
        ("bc.vbc_period_x0_ratio", po::value<std::string>()->default_value("[1.]"),
         "rate ratio of BC velocity periods at x0 side '[r0, r1, r2, ...]'")
        ("bc.vbc_period_x1_ratio", po::value<std::string>()->default_value("[1.]"),
         "rate ratio of BC velocity periods at x1 side '[r0, r1, r2, ...]'")

        ("bc.vbc_val_x1", po::value<double>(&p.bc.vbc_val_x1)->default_value(1e-9),
         "Value of boundary condition for the right/eastern side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_division_x1_min", po::value<double>(&p.bc.vbc_val_division_x1_min)->default_value(1.),
         "Ratio of first rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_division_x1_max", po::value<double>(&p.bc.vbc_val_division_x1_max)->default_value(1.),
         "Ratio of second rate change division from the top to the right boundary width (for option = 8)")
        ("bc.vbc_val_x1_ratio0", po::value<double>(&p.bc.vbc_val_x1_ratio0)->default_value(1.),
         "Ratio of the velocity at the surface of the right boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x1_ratio1", po::value<double>(&p.bc.vbc_val_x1_ratio1)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x0_min of the right boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x1_ratio2", po::value<double>(&p.bc.vbc_val_x1_ratio2)->default_value(1.),
         "Ratio of the velocity at the vbc_val_division_x1_max of the right boundary width (for option = 1 and 3)")
        ("bc.vbc_val_x1_ratio3", po::value<double>(&p.bc.vbc_val_x1_ratio3)->default_value(1.),
         "Ratio of the velocity at the bottom of the right boundary width (for option = 1 and 3)")

        ("bc.vbc_y0", po::value<int>(&p.bc.vbc_y0)->default_value(0),
         "Type of boundary condition for the southern side")
        ("bc.vbc_y1", po::value<int>(&p.bc.vbc_y1)->default_value(0),
         "Type of boundary condition for the northern side")
        ("bc.vbc_val_y0", po::value<double>(&p.bc.vbc_val_y0)->default_value(0),
         "Value of boundary condition for the southern side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_y1", po::value<double>(&p.bc.vbc_val_y1)->default_value(0),
         "Value of boundary condition for the northern side (if velocity, unit is m/s; if stress, unit is Pa)")

        ("bc.vbc_z0", po::value<int>(&p.bc.vbc_z0)->default_value(0),
         "Type of boundary condition for the bottom side")
        ("bc.vbc_z1", po::value<int>(&p.bc.vbc_z1)->default_value(0),
         "Type of boundary condition for the top side")
        ("bc.vbc_val_z0", po::value<double>(&p.bc.vbc_val_z0)->default_value(0),
         "Value of boundary condition for the bottom side (if velocity, unit is m/s; if stress, unit is Pa)")
        ("bc.vbc_val_z1", po::value<double>(&p.bc.vbc_val_z1)->default_value(0),
         "Value of boundary condition for the top side (if velocity, unit is m/s; if stress, unit is Pa)")

        ("bc.vbc_n0", po::value<int>(&p.bc.vbc_n0)->default_value(1),
         "Type of boundary condition for slant boundary #0 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n0", po::value<double>(&p.bc.vbc_val_n0)->default_value(0),
         "Value of boundary condition for slant boundary #0 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ("bc.vbc_n1", po::value<int>(&p.bc.vbc_n1)->default_value(1),
         "Type of boundary condition for slant boundary #1 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n1", po::value<double>(&p.bc.vbc_val_n1)->default_value(0),
         "Value of boundary condition for slant boundary #1 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ("bc.vbc_n2", po::value<int>(&p.bc.vbc_n2)->default_value(1),
         "Type of boundary condition for slant boundary #2 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n2", po::value<double>(&p.bc.vbc_val_n2)->default_value(0),
         "Value of boundary condition for slant boundary #2 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ("bc.vbc_n3", po::value<int>(&p.bc.vbc_n3)->default_value(1),
         "Type of boundary condition for slant boundary #3 (only type 1, 3, 11, 13 are supported).")
        ("bc.vbc_val_n3", po::value<double>(&p.bc.vbc_val_n3)->default_value(0),
         "Value of boundary condition for slant boundary #3 (if velocity, unit is m/s, "
         "outward normal direction is positive; if stress, unit is Pa)")
        ;

    cfg.add_options()
        ("ic.mattype_option", po::value<int>(&p.ic.mattype_option)->default_value(0),
         "How to set the initial material type of markers?\n"
         "0: marker's mattype is determined by regional attribute (currently can be modified only in the code or in a poly file).\n"
         "1: marker's mattype is layered.\n"
         "101: custom mattype.")
        ("ic.num_mattype_layers", po::value<int>(&p.ic.num_mattype_layers)->default_value(2),
         "Number of material layers")
        ("ic.layer_mattypes", po::value<std::string>()->default_value("[0,1]"),
         "Material type of each material layer '[d0, d1, d2, ...]', d0<=d1<=d2...")
        ("ic.mattype_layer_depths", po::value<std::string>()->default_value("[0.5]"),
         "Depths of the interfaces of each material layer '[d0, d1, d2, ...]', d0<=d1<=d2..., (in unit of zlength)")

        ("ic.weakzone_option", po::value<int>(&p.ic.weakzone_option)->default_value(1),
         "How to set the initial weak zone?\n"
         "0: no weak zone.\n"
         "1: planar weak zone with specified azimuth, inclination, halfwidth, min/max depth range, and center location.\n"
         "2: ellipsoidal weak zone with specified center location and semi-axes.\n"
         "3: Gaussian distribution point weak zone with specified center location and standard deviation.\n")
        ("ic.weakzone_plstrain", po::value<double>(&p.ic.weakzone_plstrain)->default_value(0.1),
         "Initial plastic strain in the weak zone")
        ("ic.weakzone_azimuth", po::value<double>(&p.ic.weakzone_azimuth)->default_value(0),
         "Azimuth angle (relative to +y axis) (in degree)")
        ("ic.weakzone_inclination", po::value<double>(&p.ic.weakzone_inclination)->default_value(90),
         "Inclination angle (relative to horizontal plane) (in degree)")
        ("ic.weakzone_halfwidth", po::value<double>(&p.ic.weakzone_halfwidth)->default_value(1.5),
         "Half-width of the weak zone (in unit of mesh.resolution)")
        ("ic.weakzone_y_min", po::value<double>(&p.ic.weakzone_y_min)->default_value(0),
         "Y-direction range (between 0 and 1, in unit of mesh.ylength)")
        ("ic.weakzone_y_max", po::value<double>(&p.ic.weakzone_y_max)->default_value(1),
         "Y-direction range (between 0 and 1, in unit of mesh.ylength)")
        ("ic.weakzone_depth_min", po::value<double>(&p.ic.weakzone_depth_min)->default_value(0),
         "Depth upper (shallower) range (between 0 and 1, in unit of mesh.zlength)")
        ("ic.weakzone_depth_max", po::value<double>(&p.ic.weakzone_depth_max)->default_value(1),
         "Depth lower (deeper) range (between 0 and 1, in unit of mesh.zlength)")
        ("ic.weakzone_xcenter", po::value<double>(&p.ic.weakzone_xcenter)->default_value(0.5),
         "Location of weak zone center in x direction (between 0 and 1, in unit of mesh.xlength)")
        ("ic.weakzone_ycenter", po::value<double>(&p.ic.weakzone_ycenter)->default_value(0.5),
         "Location of weak zone center in y direction (between 0 and 1, in unit of mesh.ylength)")
        ("ic.weakzone_zcenter", po::value<double>(&p.ic.weakzone_zcenter)->default_value(0.5),
         "Location of weak zone center in z direction (between 0 and 1, in unit of mesh.zlength)")
        ("ic.weakzone_xsemi_axis", po::value<double>(&p.ic.weakzone_xsemi_axis)->default_value(1e3),
         "Length of weak zone semi-axis in x direction (in meters)")
        ("ic.weakzone_ysemi_axis", po::value<double>(&p.ic.weakzone_ysemi_axis)->default_value(1e3),
         "Length of weak zone semi-axis in y direction (in meters)")
        ("ic.weakzone_zsemi_axis", po::value<double>(&p.ic.weakzone_zsemi_axis)->default_value(1e3),
         "Length of weak zone semi-axis in z direction (in meters)\n")
        ("ic.weakzone_standard_deviation", po::value<double>(&p.ic.weakzone_standard_deviation)->default_value(1e3),
         "Standard deviation of Gussian distribution point weak zone (in meters)\n")

        ("ic.temperature_option", po::value<int>(&p.ic.temperature_option)->default_value(0),
         "How to set the initial temperature?\n"
         "0: uniform half-space cooling.\n"
         "90: temperature read from an external grid\n")

        // for temperature_option = 0
	("ic.oceanic_plate_age_in_yr", po::value<double>(&p.ic.oceanic_plate_age_in_yr)->default_value(60e6),
         "Age of the oceanic plate (in years), used for the temperature profile on the plate.\n")

        // for temperature_option = 1
        ("ic.continental_plate_age_in_yr", po::value<double>(&p.ic.continental_plate_age_in_yr)->default_value(100.e6),
         "Age of the continental plate (in years), used for the temperature profile on the plate.\n")
        ("ic.radiogenic_crustal_thickness", po::value<double>(&p.ic.radiogenic_crustal_thickness)->default_value(30.e3),
         "Radiogenic crustal thickness (in meters), used for the temperature profile on the plate.\n")
        ("ic.radiogenic_folding_depth", po::value<double>(&p.ic.radiogenic_folding_depth)->default_value(10.e3),
         "Radiogenic folding depth (in meters), used for the temperature profile on the plate.\n")
        ("ic.radiogenic_heating_of_crust", po::value<double>(&p.ic.radiogenic_heating_of_crust)->default_value(1.e-9),
         "Radiogenic heating of crust (in W/kg), used for the temperature profile on the plate.\n")
        ("ic.lithospheric_thickness", po::value<double>(&p.ic.lithospheric_thickness)->default_value(100.e3),
         "Lithospheric thickness (in meters), used for the temperature profile on the plate.\n")

        // for temperature_option = 90
        ("ic.Temp_filename", po::value<std::string>(&p.ic.Temp_filename)->default_value("Thermal.dat"),
         "Filename of the input thermal field, for temperature_option=90.\n"
         "The file format is: <x-coord y-coord (if 3D) z-coord temperature>.\n"
	 "Can be generated by COMSOL.")
	("ic.Nodes_filename", po::value<std::string>(&p.ic.Nodes_filename)->default_value("Coord.dat"),
         "Filename of the mesh coordinates, for temperature_option=90.\n"
         "The nodes may be ordered differently with the thermal file. The file format is: "
         "<x-coord y-coord (if 3D) z-coord>.\n"
         "Can be generated by COMSOL.")
	("ic.Connectivity_filename", po::value<std::string>(&p.ic.Connectivity_filename)->default_value("Connectivity.dat"),
         "Filename of the mesh connectivity, for temperature_option = 90 only.\n"
         "The file format is: <#ofnode0 #ofnode1 #ofnode2 #ofnode3 (if 3D)>.\n"
         "Can be generated by COMSOL.\n")

        ("ic.isostasy_adjustment_time_in_yr", po::value<double>(&p.ic.isostasy_adjustment_time_in_yr)->default_value(0),
         "Time for spinning up isostasy adjustment.\n")
        ;

    cfg.add_options()
        ("mat.rheology_type", po::value<std::string>()->required(),
         "Type of rheology, either 'elastic', 'viscous' (experimental), 'maxwell', "
         "'elasto-plastic', or 'elasto-visco-plastic'.")
        ("mat.is_plane_strain", po::value<bool>(&p.mat.is_plane_strain)->default_value(false),
         "Is the rheology formulation in plane strain (2D elasto-plastic case only)?\n")

        ("mat.phase_change_option", po::value<int>(&p.mat.phase_change_option)->default_value(0),
         "What kind of phase changes?\n"
         "0: no phase changes.\n"
         "1: simple rules of subduction-related phase changes. See SimpleSubduction class in phasechanges.cxx for more details.\n"
         "2: simple rules of middle ocean ridge formation. See SimpleRifting class in phasechanges.cxx for more details.\n"
         "101: custom phase changes.")
        ("mat.num_materials", po::value<int>(&p.mat.nmat)->default_value(1),
         "Number of material types")
        ("mat.mattype_mantle", po::value<int>(&p.mat.mattype_mantle)->default_value(0),
         "Index of mantle material. For continental thermal gradient")
        ("mat.mattype_depleted_mantle", po::value<int>(&p.mat.mattype_depleted_mantle)->default_value(0),
         "Index of depleted mantle material. For phase change of middle ocean ridge")
        ("mat.mattype_partial_melting_mantle", po::value<int>(&p.mat.mattype_partial_melting_mantle)->default_value(0),
         "Index of parital melting mantle material. For phase change of middle ocean ridge")
        ("mat.mattype_crust", po::value<int>(&p.mat.mattype_crust)->default_value(0),
         "Index of crust material.  For continental thermal gradient")
        ("mat.mattype_oceanic_crust", po::value<int>(&p.mat.mattype_oceanic_crust)->default_value(0),
         "Index of oceanic crust. For formation of middle ocean ridge in phasechanges.cxx")
        ("mat.mattype_sed", po::value<int>(&p.mat.mattype_sed)->default_value(0),
         "Index of sediment material. Should be the same as mattype['sediment'] in 2vtk.py for sediment time convertion")
        ("mat.max_viscosity", po::value<double>(&p.mat.visc_max)->default_value(1e24),
         "Max. value of viscosity (in Pa.s)")
        ("mat.min_viscosity", po::value<double>(&p.mat.visc_min)->default_value(1e18),
         "Min. value of viscosity (in Pa.s)")
        ("mat.max_tension", po::value<double>(&p.mat.tension_max)->default_value(1e9),
         "Max. value of tensile stress (in Pa)")
        ("mat.max_thermal_diffusivity", po::value<double>(&p.mat.therm_diff_max)->default_value(5e-6),
         "Max. value of thermal diffusivity (in m^2/s)")
        ("mat.convert_rate_oceanic_crust", po::value<double>(&p.mat.convert_rate_oceanic_crust)->default_value(1.e-13),
         "Converting rate from melting to oceanic crust. For formation of middle ocean ridge in bc.cxx")

        // these parameters need to parsed later
        ("mat.rho0", po::value<std::string>()->default_value("[3210]"),
         "Density of the materials at 0 Pa and 273 K '[d0, d1, d2, ...]' (in kg/m^3)")
        ("mat.alpha", po::value<std::string>()->default_value("[3e-5]"),
         "Volumetic thermal expansion of the materials '[d0, d1, d2, ...]' (in 1/Kelvin)")

        ("mat.bulk_modulus", po::value<std::string>()->default_value("[128.2e9]"),
         "Bulk modulus of the materials '[d0, d1, d2, ...]' (in Pa)")
        ("mat.shear_modulus", po::value<std::string>()->default_value("[80.5e9]"),
         "Shear modulus of the materials '[d0, d1, d2, ...]' (in Pa)")

        ("mat.visc_exponent", po::value<std::string>()->default_value("[3.05]"),
         "Exponents of non-linear viscosity of the materials'[d0, d1, d2, ...]'")
        ("mat.visc_coefficient", po::value<std::string>()->default_value("[1.25e-1]"),
         "Pre-exponent coefficient of non-linear viscosity of the materials '[d0, d1, d2, ...]'")
        ("mat.visc_activation_energy", po::value<std::string>()->default_value("[3.76e5]"),
         "Activation energy of non-linear viscosity of the materials '[d0, d1, d2, ...]' (in J/mol)")

        ("mat.heat_capacity", po::value<std::string>()->default_value("[1e3]"),
         "Heat capacity (isobaric) of the materials '[d0, d1, d2, ...]' (in J/kg/Kelvin)")
        ("mat.therm_cond", po::value<std::string>()->default_value("[3]"),
         "Thermal conductivity of the materials '[d0, d1, d2, ...]' (in W/m/Kelvin)")

        ("mat.pls0", po::value<std::string>()->default_value("[0]"),
         "Plastic strain of the materials where weakening starts '[d0, d1, d2, ...]' (no unit)")
        ("mat.pls1", po::value<std::string>()->default_value("[0.1]"),
         "Plastic strain of the materials where weakening saturates '[d0, d1, d2, ...]' (no unit)")
        ("mat.cohesion0", po::value<std::string>()->default_value("[4e7]"),
         "Cohesion of the materials when weakening starts '[d0, d1, d2, ...]' (in Pa)")
        ("mat.cohesion1", po::value<std::string>()->default_value("[4e6]"),
         "Cohesion of the materials when weakening saturates '[d0, d1, d2, ...]' (in Pa)")
        ("mat.friction_angle0", po::value<std::string>()->default_value("[30]"),
         "Friction angle of the materials when weakening starts '[d0, d1, d2, ...]' (in degree)")
        ("mat.friction_angle1", po::value<std::string>()->default_value("[5]"),
         "Friction angle of the materials when weakening saturates '[d0, d1, d2, ...]' (in degree)")
        ("mat.dilation_angle0", po::value<std::string>()->default_value("[0]"),
         "Dilation angle of the materials when weakening starts '[d0, d1, d2, ...]' (in degree)")
        ("mat.dilation_angle1", po::value<std::string>()->default_value("[0]"),
         "Dilation angle of the materials when weakening saturates '[d0, d1, d2, ...]' (in degree)")
        ;

    /* These parameters will enable additional debugging output. DO NOT document these parameters
       in defaults.cfg. */
    cfg.add_options()
        ("debug.dt", po::value<bool>(&p.debug.dt)->default_value(false),
         "Print all dt criteria")
//        ("debug.has_two_layers_for", po::value<bool>(&p.debug.has_two_layers_for)->default_value(true),
//         "Use two layers of for elem to avoid race condition.")
        ;
}


static void read_parameters_from_file
(const char* filename,
 const po::options_description cfg,
 po::variables_map &vm)
{
    try {
        po::store(po::parse_config_file<char>(filename, cfg), vm);
        po::notify(vm);
    }
    catch (const boost::program_options::multiple_occurrences& e) {
        std::cerr << e.what() << " from option: " << e.get_option_name() << '\n';
        std::exit(1);
    }
    catch (std::exception& e) {
        std::cerr << "Error reading config_file '" << filename << "'\n";
        std::cerr << e.what() << "\n";
        std::exit(1);
    }
}


template<class T>
static int read_numbers(const std::string &input, std::vector<T> &vec, int len)
{
    /* Read 'len' numbers from input.
     * The format of input must be '[n0, n1, n2]' or '[n0, n1, n2,]' (with a trailing ,),
     * for len=3.
     */

    std::istringstream stream(input);
    vec.resize(len);

    char sentinel;

    stream >> sentinel;
    if (sentinel != '[') return 1;

    for (int i=0; i<len; ++i) {
        stream >> vec[i];

        if (i == len-1) break;

        // consume ','
        char sep;
        stream >> sep;
        if (sep != ',') return 1;
    }

    stream >> sentinel;
    if (sentinel == ',') stream >> sentinel;
    if (sentinel != ']') return 1;

    if (! stream.good()) return 1;

    // success
    return 0;
}


template<class T>
static void get_numbers(const po::variables_map &vm, const char *name,
                        std::vector<T> &values, int len, int optional_size=0)
{
    if ( ! vm.count(name) ) {
        std::cerr << "Error: " << name << " is not provided.\n";
        std::exit(1);
    }

    std::string str = vm[name].as<std::string>();
    int err = read_numbers(str, values, len);
    if (err && optional_size) {
        err = read_numbers(str, values, optional_size);
    }

    if (err) {
        std::cerr << "Error: incorrect format for " << name << ",\n"
                  << "       must be '[d0, d1, d2, ...]'\n";
        std::exit(1);
    }
}


static void validate_parameters(const po::variables_map &vm, Param &p)
{
    std::cout << "Checking consistency of input parameters...\n";

    //
    // stopping condition and output interval are based on either model time or step
    //
    if ( ! (vm.count("sim.max_steps") || vm.count("sim.max_time_in_yr")) ) {
        std::cerr << "Must provide either sim.max_steps or sim.max_time_in_yr\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.max_steps") )
        p.sim.max_steps = std::numeric_limits<int>::max();
    if ( ! vm.count("sim.max_time_in_yr") )
        p.sim.max_time_in_yr = std::numeric_limits<double>::max();

    if ( ! (vm.count("sim.output_step_interval") || vm.count("sim.output_time_interval_in_yr")) ) {
        std::cerr << "Must provide either sim.output_step_interval or sim.output_time_interval_in_yr\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.output_step_interval") )
        p.sim.output_step_interval = std::numeric_limits<int>::max();
    if ( ! vm.count("sim.output_time_interval_in_yr") )
        p.sim.output_time_interval_in_yr = std::numeric_limits<double>::max();

    //
    // These parameters are required when restarting
    //
    if (p.sim.is_restarting) {
        if ( ! vm.count("sim.restarting_from_modelname") ) {
            std::cerr << "Must provide sim.restarting_from_modelname when restarting.\n";
            std::exit(1);
        }
        if ( ! vm.count("sim.restarting_from_frame") ) {
            std::cerr << "Must provide sim.restarting_from_frame when restarting.\n";
            std::exit(1);
        }
    }

    if (p.sim.is_outputting_averaged_fields == true)
        if (p.sim.output_step_interval%p.mesh.quality_check_step_interval !=0) {
            std::cerr << "sim.output_step_interval must be a multiple of mesh.quality_check_step_interval!.\n";
            std::exit(1);
    }

    //
    // these parameters are required in mesh.meshing_option == 2
    //
    if (p.mesh.meshing_option == 2) {
        if ( ! vm.count("mesh.refined_zonex") ||
#ifdef THREED
             ! vm.count("mesh.refined_zoney") ||
#endif
             ! vm.count("mesh.refined_zonez") ) {
        std::cerr << "Must provide mesh.refined_zonex, "
#ifdef THREED
                  << "mesh.refined_zoney, "
#endif
                  << "mesh.refined_zonez.\n";
        std::exit(1);
        }

        /* get 2 numbers from the string */
        double_vec tmp;
        int err;
        std::string str;
        str = vm["mesh.refined_zonex"].as<std::string>();
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refine_zonex,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonex.first = tmp[0];
        p.mesh.refined_zonex.second = tmp[1];
#ifdef THREED
        str = vm["mesh.refined_zoney"].as<std::string>();
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refine_zoney,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zoney.first = tmp[0];
        p.mesh.refined_zoney.second = tmp[1];
#endif
        str = vm["mesh.refined_zonez"].as<std::string>();
        err = read_numbers(str, tmp, 2);
        if (err || tmp[0] < 0 || tmp[1] > 1 || tmp[0] > tmp[1]) {
            std::cerr << "Error: incorrect value for mesh.refine_zonez,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonez.first = tmp[0];
        p.mesh.refined_zonez.second = tmp[1];
    }

    if (p.mesh.smallest_size > p.mesh.largest_size) {
        std::cerr << "Error: mesh.smallest_size is greater than mesh.largest_size.\n";
        std::exit(1);
    }

#ifdef THREED
    if (p.mesh.remeshing_option == 2) {
        std::cerr << "Error: mesh.remeshing_option=2 is not available in 3D.\n";
        std::exit(1);
    }
#endif

    //
    // bc
    //
    {
        if ( p.bc.has_winkler_foundation && p.control.gravity == 0 ) {
            p.bc.has_winkler_foundation = 0;
            std::cerr << "Warning: no gravity, Winkler foundation is turned off.\n";
        }
        if ( p.bc.has_winkler_foundation && p.bc.vbc_z0 != 0 ) {
            p.bc.vbc_z0 = 0;
            std::cerr << "Warning: Winkler foundation is turned on, setting bc.vbc_z0 to 0.\n";
        }
        if ( p.bc.has_water_loading && p.control.gravity == 0 ) {
            p.bc.has_water_loading = 0;
            std::cerr << "Warning: no gravity, water loading is turned off.\n";
        }
        if ( p.bc.has_water_loading && p.bc.vbc_z1 != 0 ) {
            p.bc.vbc_z1 = 0;
            std::cerr << "Warning: water loading is turned on, setting bc.vbc_z1 to 0.\n";
        }

        if ( p.bc.vbc_z0 > 3) {
            std::cerr << "Error: bc.vbc_z0 is not 0, 1, 2, or 3.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_z1 > 3) {
            std::cerr << "Error: bc.vbc_z0 is not 0, 1, 2, or 3.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n0 != 1 && p.bc.vbc_n0 != 3 && p.bc.vbc_n0 != 11 && p.bc.vbc_n0 != 13 ) {
            std::cerr << "Error: bc.vbc_n0 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n1 != 1 && p.bc.vbc_n1 != 3 && p.bc.vbc_n1 != 11 && p.bc.vbc_n1 != 13 ) {
            std::cerr << "Error: bc.vbc_n1 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n2 != 1 && p.bc.vbc_n2 != 3 && p.bc.vbc_n2 != 11 && p.bc.vbc_n2 != 13 ) {
            std::cerr << "Error: bc.vbc_n2 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
        if ( p.bc.vbc_n3 != 1 && p.bc.vbc_n3 != 3 && p.bc.vbc_n3 != 11 && p.bc.vbc_n3 != 13 ) {
            std::cerr << "Error: bc.vbc_n3 is not 1, 3, 11, or 13.\n";
            std::exit(1);
        }
    }

    //
    // control
    //
    {
        if ( p.control.dt_fraction < 0 || p.control.dt_fraction > 1 ) {
            std::cerr << "Error: control.dt_fraction must be between 0 and 1.\n";
            std::exit(1);
        }
        if ( p.control.damping_factor < 0 || p.control.damping_factor > 1 ) {
            std::cerr << "Error: control.damping_factor must be between 0 and 1.\n";
            std::exit(1);
        }

    }

    //
    // ic
    //
    {
        if ( p.ic.mattype_option == 1) {
            get_numbers(vm, "ic.layer_mattypes", p.ic.layer_mattypes, p.ic.num_mattype_layers);
            get_numbers(vm, "ic.mattype_layer_depths", p.ic.mattype_layer_depths, p.ic.num_mattype_layers-1);
            // mattype_layer_depths must be already sorted
            if (! std::is_sorted(p.ic.mattype_layer_depths.begin(), p.ic.mattype_layer_depths.end())) {
                std::cerr << "Error: the content of ic.mattype_layer_depths is not ordered from"
                    " small to big values.\n";
                std::exit(1);
            }
        }
    }

    //
    // marker
    //
    {

    }

    //
    // material properties
    //
    {
        std::string str = vm["mat.rheology_type"].as<std::string>();
        if (str == std::string("elastic"))
            p.mat.rheol_type = MatProps::rh_elastic;
        else if (str == std::string("viscous"))
            p.mat.rheol_type = MatProps::rh_viscous;
        else if (str == std::string("maxwell"))
            p.mat.rheol_type = MatProps::rh_maxwell;
        else if (str == std::string("elasto-plastic"))
            p.mat.rheol_type = MatProps::rh_ep;
        else if (str == std::string("elasto-visco-plastic"))
            p.mat.rheol_type = MatProps::rh_evp;
        else {
            std::cerr << "Error: unknown rheology: '" << str << "'\n";
            std::exit(1);
        }

#ifdef THREED
        if ( p.mat.is_plane_strain ) {
            p.mat.is_plane_strain = false;
            std::cerr << "Warning: mat.is_plane_strain is not avaiable in 3D.\n";
        }
#endif

        if (p.mat.phase_change_option != 0 && p.mat.nmat == 1) {
            std::cerr << "Error: mat.phase_change_option is chosen, but mat.num_materials is 1.\n";
            std::exit(1);
        }
        if (p.mat.phase_change_option == 1 && p.mat.nmat < 8) {
            std::cerr << "Error: mat.phase_change_option is 1, but mat.num_materials is less than 8.\n";
            std::exit(1);
        }

        if (p.mat.nmat < 1) {
            std::cerr << "Error: mat.num_materials must be greater than 0.\n";
            std::exit(1);
        }

        if (p.mat.nmat == 1 && p.control.ref_pressure_option != 0) {
            p.control.ref_pressure_option = 0;
            std::cerr << "Warning: mat.num_materials is 1, using simplest control.ref_pressure_option.\n";
        }
        if (p.mat.nmat == 1 && p.markers.replenishment_option != 1) {
            p.markers.replenishment_option = 1;
            std::cerr << "Warning: mat.num_materials is 1, using simplest markers.replenishment_option.\n";
        }

        get_numbers(vm, "bc.vbc_period_x0_time_in_yr", p.bc.vbc_period_x0_time_in_yr, p.bc.num_vbc_period_x0, 1);
        get_numbers(vm, "bc.vbc_period_x1_time_in_yr", p.bc.vbc_period_x1_time_in_yr, p.bc.num_vbc_period_x1, 1);
        get_numbers(vm, "bc.vbc_period_x0_ratio", p.bc.vbc_period_x0_ratio, p.bc.num_vbc_period_x0, 1);
        get_numbers(vm, "bc.vbc_period_x1_ratio", p.bc.vbc_period_x1_ratio, p.bc.num_vbc_period_x1, 1);

        get_numbers(vm, "mat.rho0", p.mat.rho0, p.mat.nmat, 1);
        get_numbers(vm, "mat.alpha", p.mat.alpha, p.mat.nmat, 1);

        get_numbers(vm, "mat.bulk_modulus", p.mat.bulk_modulus, p.mat.nmat, 1);
        get_numbers(vm, "mat.shear_modulus", p.mat.shear_modulus, p.mat.nmat, 1);

        get_numbers(vm, "mat.visc_exponent", p.mat.visc_exponent, p.mat.nmat, 1);
        get_numbers(vm, "mat.visc_coefficient", p.mat.visc_coefficient, p.mat.nmat, 1);
        get_numbers(vm, "mat.visc_activation_energy", p.mat.visc_activation_energy, p.mat.nmat, 1);

        get_numbers(vm, "mat.heat_capacity", p.mat.heat_capacity, p.mat.nmat, 1);
        get_numbers(vm, "mat.therm_cond", p.mat.therm_cond, p.mat.nmat, 1);

        get_numbers(vm, "mat.pls0", p.mat.pls0, p.mat.nmat, 1);
        get_numbers(vm, "mat.pls1", p.mat.pls1, p.mat.nmat, 1);
        get_numbers(vm, "mat.cohesion0", p.mat.cohesion0, p.mat.nmat, 1);
        get_numbers(vm, "mat.cohesion1", p.mat.cohesion1, p.mat.nmat, 1);
        get_numbers(vm, "mat.friction_angle0", p.mat.friction_angle0, p.mat.nmat, 1);
        get_numbers(vm, "mat.friction_angle1", p.mat.friction_angle1, p.mat.nmat, 1);
        get_numbers(vm, "mat.dilation_angle0", p.mat.dilation_angle0, p.mat.nmat, 1);
        get_numbers(vm, "mat.dilation_angle1", p.mat.dilation_angle1, p.mat.nmat, 1);
    }

}


void get_input_parameters(const char* filename, Param& p)
{
    po::options_description cfg("Config file options");
    po::variables_map vm;

    declare_parameters(cfg, p);
    // print help message
    if (std::strncmp(filename, "-h", 3) == 0 ||
        std::strncmp(filename, "--help", 7) == 0) {
        std::cout << cfg;
        std::exit(0);
    }
    read_parameters_from_file(filename, cfg, vm);
    validate_parameters(vm, p);
}
