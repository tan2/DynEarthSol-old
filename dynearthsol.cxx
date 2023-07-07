#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include <limits>


#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "binaryio.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "ic.hpp"
#include "input.hpp"
#include "matprops.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "phasechanges.hpp"
#include "remeshing.hpp"
#include "rheology.hpp"
#include "utils.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

void init_var(const Param& param, Variables& var)
{
    var.time = 0;
    var.steps = 0;
    var.func_time.output_time = 0;
    var.func_time.remesh_time = 0;
    var.func_time.start_time = get_nanoseconds();

    for (int i=0;i<nbdrytypes;++i)
        var.bfacets[i] = new std::vector< std::pair<int,int> >;
    for (int i=0;i<nbdrytypes;++i)
        var.bnodes[i] = new int_vec;
    var.bnormals = new array_t(nbdrytypes);

    if (param.control.characteristic_speed == 0)
        var.max_vbc_val = find_max_vbc(param.bc);
        // todo max_vbc_val change with boundary period.
    else
        var.max_vbc_val = param.control.characteristic_speed;

    // XXX: Hard coded boundary flag. If the order of ibound?? is changed
    //      in the future, the following lines have to be updated as well.
    var.vbc_types[0] = param.bc.vbc_x0;
    var.vbc_types[1] = param.bc.vbc_x1;
    var.vbc_types[2] = param.bc.vbc_y0;
    var.vbc_types[3] = param.bc.vbc_y1;
    var.vbc_types[4] = param.bc.vbc_z0;
    var.vbc_types[5] = param.bc.vbc_z1;
    var.vbc_types[6] = param.bc.vbc_n0;
    var.vbc_types[7] = param.bc.vbc_n1;
    var.vbc_types[8] = param.bc.vbc_n2;
    var.vbc_types[9] = param.bc.vbc_n3;

    var.vbc_values[0] = param.bc.vbc_val_x0;
    var.vbc_values[1] = param.bc.vbc_val_x1;
    var.vbc_values[2] = param.bc.vbc_val_y0;
    var.vbc_values[3] = param.bc.vbc_val_y1;
    var.vbc_values[4] = param.bc.vbc_val_z0;
    var.vbc_values[5] = param.bc.vbc_val_z1;
    var.vbc_values[6] = param.bc.vbc_val_n0;
    var.vbc_values[7] = param.bc.vbc_val_n1;
    var.vbc_values[8] = param.bc.vbc_val_n2;
    var.vbc_values[9] = param.bc.vbc_val_n3;

    var.vbc_vertical_div_x0[0] = 0.;
    var.vbc_vertical_div_x0[1] = param.bc.vbc_val_division_x0_min;
    var.vbc_vertical_div_x0[2] = param.bc.vbc_val_division_x0_max;
    var.vbc_vertical_div_x0[3] = 1.;
    var.vbc_vertical_div_x1[0] = 0.;
    var.vbc_vertical_div_x1[1] = param.bc.vbc_val_division_x1_min;
    var.vbc_vertical_div_x1[2] = param.bc.vbc_val_division_x1_max;
    var.vbc_vertical_div_x1[3] = 1.;

    var.vbc_vertical_ratio_x0[0] = param.bc.vbc_val_x0_ratio0;
    var.vbc_vertical_ratio_x0[1] = param.bc.vbc_val_x0_ratio1;
    var.vbc_vertical_ratio_x0[2] = param.bc.vbc_val_x0_ratio2;
    var.vbc_vertical_ratio_x0[3] = param.bc.vbc_val_x0_ratio3;
    var.vbc_vertical_ratio_x1[0] = param.bc.vbc_val_x1_ratio0;
    var.vbc_vertical_ratio_x1[1] = param.bc.vbc_val_x1_ratio1;
    var.vbc_vertical_ratio_x1[2] = param.bc.vbc_val_x1_ratio2;
    var.vbc_vertical_ratio_x1[3] = param.bc.vbc_val_x1_ratio3;
}


void init(const Param& param, Variables& var)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    std::cout << "Initializing mesh and field data...\n";

    create_new_mesh(param, var);
    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);
    create_elemmarkers(param, var);
    create_markers(param, var);

    allocate_variables(param, var);
//    var.markersets[0]->create_marker_in_elem(var);
//    var.markersets[0]->create_melt_markers(param.mat.mattype_partial_melting_mantle,var.melt_markers);

    create_top_elems(var);
    create_surface_info(param,var,var.surfinfo);

    for(int i=0; i<var.nnode; i++)
        for(int d=0; d<NDIMS; d++)
            (*var.coord0)[i][d] = (*var.coord)[i][d];

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    *var.volume_old = *var.volume;
    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.tmp_result);
    compute_shape_fn(var, *var.shpdx, *var.shpdy, *var.shpdz);

    create_boundary_normals(var, *var.bnormals, var.edge_vectors);
    apply_vbcs(param, var, *var.vel);

    // temperature should be init'd before stress and strain
    initial_temperature(param, var, *var.temperature);
    initial_stress_state(param, var, *var.stress, *var.stressyy, *var.strain, var.compensation_pressure);
    initial_weak_zone(param, var, *var.plstrain);

    phase_changes_init(param, var);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void restart(const Param& param, Variables& var)
{
    std::cout << "Initializing mesh and field data from checkpoints...\n";

    /* Reading info file */
    {
        char filename[256];
        std::snprintf(filename, 255, "%s.info", param.sim.restarting_from_modelname.c_str());
        std::FILE *f = std::fopen(filename, "r");
        int frame, steps, nnode, nelem, nseg;
        while (1) {
            int n = std::fscanf(f, "%d %d %*f %*f %*f %d %d %d\n",
                                &frame, &steps, &nnode, &nelem, &nseg);
            if (n != 5) {
                std::cerr << "Error: reading info file: " << filename << '\n';
                std::exit(2);
            }
            if (frame == param.sim.restarting_from_frame)
                break;
        }

        var.steps = steps;
        var.nnode = nnode;
        var.nelem = nelem;
        var.nseg = nseg;

        std::fclose(f);
    }

    char filename_save[256];
    std::snprintf(filename_save, 255, "%s.save.%06d",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    BinaryInput bin_save(filename_save);
    std::cout << "  Reading " << filename_save << "...\n";

    char filename_chkpt[256];
    std::snprintf(filename_chkpt, 255, "%s.chkpt.%06d",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    BinaryInput bin_chkpt(filename_chkpt);
    std::cout << "  Reading " << filename_chkpt << "...\n";

    //
    // Following the same procedure in init()
    //

    // Reading mesh, replacing create_new_mesh()
    {
        var.coord = new array_t(var.nnode);
        bin_save.read_array(*var.coord, "coordinate");
        var.connectivity = new conn_t(var.nelem);
        bin_save.read_array(*var.connectivity, "connectivity");

        var.segment = new segment_t(var.nseg);
        bin_chkpt.read_array(*var.segment, "segment");
        var.segflag = new segflag_t(var.nseg);
        bin_chkpt.read_array(*var.segflag, "segflag");
        // Note: regattr is not needed for restarting
        // var.regattr = new regattr_t(var.nelem);
        // bin_chkpt.read_array(*var.regattr, "regattr", var.nelem);
    }

    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);
    create_elemmarkers(param, var);

    // Replacing create_markers()
    var.markersets.push_back(new MarkerSet(param, var, bin_chkpt, std::string("markerset")));
    if (param.control.has_hydration_processes) {
        var.hydrous_marker_index = var.markersets.size();
        var.markersets.push_back(new MarkerSet(param, var, bin_chkpt, std::string("hydrous-markerset")));
    }

    allocate_variables(param, var);

    create_top_elems(var);
//    var.markersets[0]->create_marker_in_elem(var);
//    var.markersets[0]->create_melt_markers(param.mat.mattype_partial_melting_mantle,var.melt_markers);
    create_surface_info(param,var,var.surfinfo);

    bin_save.read_array(*var.coord0, "coord0");

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    bin_chkpt.read_array(*var.volume_old, "volume_old");
    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.tmp_result);

    compute_shape_fn(var, *var.shpdx, *var.shpdy, *var.shpdz);

    create_boundary_normals(var, *var.bnormals, var.edge_vectors);

    // Initializing field variables
    {
        bin_save.read_array(*var.vel, "velocity");
        bin_save.read_array(*var.temperature, "temperature");
        bin_save.read_array(*var.strain_rate, "strain-rate");
        bin_save.read_array(*var.strain, "strain");
        bin_save.read_array(*var.stress, "stress");
        bin_save.read_array(*var.plstrain, "plastic strain");

        if (param.mat.is_plane_strain)
            bin_chkpt.read_array(*var.stressyy, "stressyy");
    }
    apply_vbcs(param, var, *var.vel);

    // Misc. items
    {
        double_vec tmp(2);
        bin_chkpt.read_array(tmp, "time compensation_pressure");
        var.time = tmp[0];
        var.compensation_pressure = tmp[1];

        // the following fields are not required for restarting
        bin_save.read_array(*var.force, "force");
    }

    phase_changes_init(param, var);
}


void update_mesh(const Param& param, Variables& var)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    update_coordinate(var, *var.coord);
    surface_processes(param, var, *var.coord, *var.stress, *var.strain, *var.strain_rate, \
                      *var.plstrain, var.surfinfo, var.markersets, *var.elemmarkers);
//    var.markersets[0]->update_marker_in_elem(var);
//    var.markersets[0]->create_melt_markers(param.mat.mattype_partial_melting_mantle,var.melt_markers);

#ifdef USE_NPROF
    nvtxRangePushA("swap vectors");
#endif
    #pragma serial
    {
        double_vec *tmp = var.volume;
        var.volume = var.volume_old;
        var.volume_old = tmp;
    }
//    var.volume->swap(*var.volume_old);
#ifdef USE_NPROF
    nvtxRangePop();
#endif

    compute_volume(var, *var.volume);
    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.tmp_result);
    compute_shape_fn(var, *var.shpdx, *var.shpdy, *var.shpdz);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void isostasy_adjustment(const Param &param, Variables &var)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    std::cout << "Adjusting isostasy for " << param.ic.isostasy_adjustment_time_in_yr << " yrs...\n";

    var.dt = compute_dt(param, var);
    int iso_steps = param.ic.isostasy_adjustment_time_in_yr*YEAR2SEC / var.dt;

    for (int n=0; n<iso_steps; n++) {
        update_strain_rate(var, *var.strain_rate);
        compute_dvoldt(var, *var.ntmp, *var.tmp_result_sg);
        compute_edvoldt(var, *var.ntmp, *var.edvoldt);
        update_stress(param ,var, *var.stress, *var.stressyy, *var.dpressure,
            *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
            *var.strain_rate);
        update_force(param, var, *var.force, *var.tmp_result);
        update_velocity(var, *var.vel);

        // do not apply vbc to allow free boundary

        // displacment is vertical only
        #pragma omp parallel for default(none)          \
            shared(var, param)
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0; j<NDIMS-1; ++j) {
                (*var.vel)[i][j] = 0;
            }
            if (param.bc.has_winkler_foundation == false &&
                (*var.bcflag)[i] & BOUNDZ0) {
                // holding bottom surface fixed
                (*var.vel)[i][NDIMS-1] = 0;
            }
        }

        update_mesh(param, var);

    }
    std::cout << "Adjusted isostasy for " << iso_steps << " steps.\n";
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


int main(int argc, const char* argv[])
{
    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        std::cout << "       " << argv[0] << " -h or --help\n";
        return -1;
    }

    Param param;
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    init_var(param, var);

    Output output(param, var.func_time.start_time,
                  (param.sim.is_restarting) ? param.sim.restarting_from_frame : 0);

    if (! param.sim.is_restarting) {
        init(param, var);

        if (param.ic.isostasy_adjustment_time_in_yr > 0) {
            // output.write_exact(var);
            isostasy_adjustment(param, var);
        }
        if (param.sim.has_initial_checkpoint)
            output.write_checkpoint(param, var);
    }
    else {
        restart(param, var);
    }

    var.dt = compute_dt(param, var);
    output.write_exact(var);

    double starting_time = var.time; // var.time & var.steps might be set in restart()
    double starting_step = var.steps;
    int next_regular_frame = 1;  // excluding frames due to output_during_remeshing

    std::cout << "Starting simulation...\n";
    do {
#ifdef USE_NPROF
        nvtxRangePush("dynearthsol");
#endif
        var.steps ++;
        var.time += var.dt;

        if (param.control.has_thermal_diffusion)
            update_temperature(param, var, *var.temperature, *var.ntmp, *var.tmp_result);

        update_strain_rate(var, *var.strain_rate);
        compute_dvoldt(var, *var.ntmp, *var.tmp_result_sg);
        compute_edvoldt(var, *var.ntmp, *var.edvoldt);

        update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
            *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
            *var.strain_rate);

	// Nodal Mixed Discretization For Stress
        if (param.control.is_using_mixed_stress)
            NMD_stress(param, var, *var.ntmp, *var.stress, *var.tmp_result_sg);

        update_force(param, var, *var.force, *var.tmp_result);
        update_velocity(var, *var.vel);
        apply_vbcs(param, var, *var.vel);
        update_mesh(param, var);

        // elastic stress/strain are objective (frame-indifferent)
        if (var.mat->rheol_type & MatProps::rh_elastic)
            rotate_stress(var, *var.stress, *var.strain);

        const int slow_updates_interval = 10;
        if (var.steps % slow_updates_interval == 0) {
            // The functions inside this if-block are expensive in computation is expensive,
            // and only changes slowly. Don't have to do it every time step
            phase_changes(param, var);

            if (param.control.has_hydration_processes)
                advect_hydrous_markers(param, var, 10*var.dt,
                                       *var.markersets[var.hydrous_marker_index],
                                       *var.hydrous_elemmarkers);

            var.dt = compute_dt(param, var);
        }

        if (param.sim.is_outputting_averaged_fields)
            output.average_fields(var);

        if (( (param.sim.output_step_interval != std::numeric_limits<int>::max() &&
               (var.steps - starting_step) == next_regular_frame * param.sim.output_step_interval)
              ||
              (param.sim.output_time_interval_in_yr != std::numeric_limits<double>::max() &&
               (var.time - starting_time) > next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC)
             )
            // time or step output requirements are met
            &&
            ((! param.sim.is_outputting_averaged_fields) ||
                (param.sim.is_outputting_averaged_fields &&
                    (var.steps % param.mesh.quality_check_step_interval == 0)))
            // When is_outputting_averaged_fields is turned on, the output cannot be
            // done at arbitrary time steps.
            ) {

            if (next_regular_frame % param.sim.checkpoint_frame_interval == 0)
                output.write_checkpoint(param, var);

            int64_t time_tmp = get_nanoseconds();
            output.write(var);
            var.func_time.output_time += get_nanoseconds() - time_tmp;

            next_regular_frame ++;
        }

        if (var.steps % param.mesh.quality_check_step_interval == 0) {
            int quality_is_bad, bad_quality_index;
            quality_is_bad = bad_mesh_quality(param, var, bad_quality_index);
            if (quality_is_bad) {

                if (param.sim.has_output_during_remeshing) {
                    int64_t time_tmp = get_nanoseconds();
                    output.write_exact(var);
                    var.func_time.output_time += get_nanoseconds() - time_tmp;
                }

                int64_t time_tmp = get_nanoseconds();
                remesh(param, var, quality_is_bad);
                var.func_time.remesh_time += get_nanoseconds() - time_tmp;

                if (param.sim.has_output_during_remeshing) {
                    int64_t time_tmp = get_nanoseconds();
                    output.write_exact(var);
                    var.func_time.output_time += get_nanoseconds() - time_tmp;
                }
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    std::cout << "Ending simulation.\n";
    int64_t duration_ns = get_nanoseconds() - var.func_time.start_time;
    std::cout << "Time summary...\n  Execute: ";
    print_time_ns(duration_ns);
    std::cout << "\n  Remesh : ";
    print_time_ns(var.func_time.remesh_time);
    std::cout << " (" <<  std::setw(5) << std::fixed << std::setprecision(2) << std::setfill(' ')
        << 100.*var.func_time.remesh_time/duration_ns << "%)\n";
    std::cout << "  Output : ";
    print_time_ns(var.func_time.output_time);
    std::cout << " (" <<  std::setw(5) <<  std::fixed << std::setprecision(2) << std::setfill(' ')
        << 100./var.func_time.output_time/duration_ns << "%)\n";
    return 0;
}
