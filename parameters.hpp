#ifndef DYNEARTHSOL3D_PARAMETERS_HPP
#define DYNEARTHSOL3D_PARAMETERS_HPP

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "constants.hpp"
#include "array2d.hpp"

typedef std::pair<double,double> double_pair;

typedef std::vector<double> double_vec;
typedef std::vector<int> int_vec;
typedef std::vector<int_vec> int_vec2D;
typedef std::vector<uint> uint_vec;

typedef Array2D<double,NDIMS> array_t;
typedef Array2D<double,NSTR> tensor_t;
typedef Array2D<double,NODES_PER_ELEM> shapefn;
typedef Array2D<double,1> regattr_t;
typedef Array2D<double,NODES_PER_ELEM*3> elem_cache;

typedef Array2D<int,NODES_PER_ELEM> conn_t;
typedef Array2D<int,NDIMS> segment_t;
typedef Array2D<int,1> segflag_t;

// forward declaration
class PhaseChange;


//
// Structures for input parameters
//
struct Sim {
    double max_time_in_yr;
    double output_time_interval_in_yr;
    int max_steps;
    int output_step_interval;
    int checkpoint_frame_interval;
    int restarting_from_frame;
    bool is_outputting_averaged_fields;
    bool is_restarting;
    bool has_initial_checkpoint;
    bool has_output_during_remeshing;
    bool has_marker_output;

    std::string modelname;
    std::string restarting_from_modelname;
};

struct Mesh {
    int meshing_option;
    int meshing_verbosity;
    int tetgen_optlevel;
    int quality_check_step_interval;

    double xlength, ylength, zlength;
    double resolution;
    double smallest_size;
    double largest_size;
    double min_angle;  // for 2D only
    double min_tet_angle, max_ratio; // for 3D only
    double min_quality;
    double max_boundary_distortion;

    double_pair refined_zonex, refined_zoney, refined_zonez;
    std::string poly_filename;
    std::string exo_filename;

    bool is_discarding_internal_segments;
    int remeshing_option;

    // Parameters for mesh optimizer MMG
    int mmg_debug;
    int mmg_verbose;
    double mmg_hmax_factor;
    double mmg_hmin_factor;
    double mmg_hausd_factor;
};

struct Control {
    double gravity;
    double characteristic_speed;
    double inertial_scaling;
    double dt_fraction;
    double fixed_dt;
    double damping_factor;
    int damping_option;
    int ref_pressure_option;

    int surface_process_option;
    double surface_diffusivity;

    bool is_quasi_static;
    bool has_thermal_diffusion;

    bool has_hydration_processes;
    double hydration_migration_speed;
};

struct BC {
    double surface_temperature;
    double mantle_temperature;

    double winkler_delta_rho;
    bool has_winkler_foundation;

    double elastic_foundation_constant;
    bool has_elastic_foundation;

    bool has_water_loading;

    int vbc_x0;
    int vbc_x1;
    int vbc_y0;
    int vbc_y1;
    int vbc_z0;
    int vbc_z1;
    int vbc_n0;
    int vbc_n1;
    int vbc_n2;
    int vbc_n3;

    double vbc_val_x0;
    double vbc_val_x1;
    double vbc_val_y0;
    double vbc_val_y1;
    double vbc_val_z0;
    double vbc_val_z1;
    double vbc_val_n0;
    double vbc_val_n1;
    double vbc_val_n2;
    double vbc_val_n3;
};

struct IC {
    int mattype_option;
    int num_mattype_layers;
    int_vec layer_mattypes;
    double_vec mattype_layer_depths;

    int weakzone_option;
    double weakzone_plstrain;
    double weakzone_azimuth;
    double weakzone_inclination;
    double weakzone_halfwidth;
    double weakzone_y_min;
    double weakzone_y_max;
    double weakzone_depth_min;
    double weakzone_depth_max;
    double weakzone_xcenter;
    double weakzone_ycenter;
    double weakzone_zcenter;
    double weakzone_xsemi_axis;
    double weakzone_ysemi_axis;
    double weakzone_zsemi_axis;

    int temperature_option;
    std::string Temp_filename;
    std::string Nodes_filename;
    std::string Connectivity_filename;
    double oceanic_plate_age_in_yr;

    double isostasy_adjustment_time_in_yr;
};

struct Mat {
    int rheol_type;
    int phase_change_option;
    int nmat;
    bool is_plane_strain;
    double visc_min;
    double visc_max;
    double tension_max;
    double therm_diff_max;

    double_vec rho0;
    double_vec alpha;

    double_vec bulk_modulus;
    double_vec shear_modulus;

    double_vec visc_exponent;
    double_vec visc_coefficient;
    double_vec visc_activation_energy;

    double_vec heat_capacity;
    double_vec therm_cond;

    // plastic parameters
    double_vec pls0, pls1;
    double_vec cohesion0, cohesion1;
    double_vec friction_angle0, friction_angle1;
    double_vec dilation_angle0, dilation_angle1;
};

struct Time {
    int64_t remesh_time;
    int64_t output_time;
    int64_t start_time;
};

struct Markers {
    int init_marker_option;
    int markers_per_element;
    int min_num_markers_in_element;
    int replenishment_option;
    uint random_seed;
    double init_marker_spacing;
};

struct Debug {
    bool dt;
};

struct Param {
    Sim sim;
    Mesh mesh;
    Control control;
    BC bc;
    IC ic;
    Mat mat;
    Markers markers;
    Debug debug;
};


//
// Structures for model variables
//
class MatProps;
class MarkerSet;
struct Variables {
    double time;
    double dt;
    int steps;
    Time func_time;

    int nnode;
    int nelem;
    int nseg;

    double max_vbc_val;
    double compensation_pressure;

    // These 5 arrays are allocated by external library
    array_t *coord;
    conn_t *connectivity;
    segment_t *segment;
    segflag_t *segflag;
    regattr_t *regattr;

    uint_vec *bcflag;
    int_vec bnodes[nbdrytypes];
    std::vector< std::pair<int,int> > bfacets[nbdrytypes];
    double bnormals[nbdrytypes][NDIMS];
    int vbc_types[nbdrytypes];
    double vbc_values[nbdrytypes];
    std::map<std::pair<int,int>, double*> edge_vectors;

    int_vec2D *support;

    double_vec *volume, *volume_old, *volume_n;
    double_vec *mass, *tmass;
    double_vec *edvoldt;
    double_vec *temperature, *plstrain, *delta_plstrain;
    double_vec *stressyy, *dpressure;
    double_vec *ntmp;

    array_t *vel, *force, *coord0;
    tensor_t *strain_rate, *strain, *stress;
    shapefn *shpdx, *shpdy, *shpdz;
    elem_cache *tmp_result;
    double_vec *tmp_result_sg;

    MatProps *mat;

    std::vector<MarkerSet*> markersets;
    int hydrous_marker_index;

    int_vec2D *elemmarkers; // for marksersets[0] (mattype markers)
    Array2D<int,1> *hydrous_elemmarkers; // for markersets[hydrous_marker_index] (hydrous markers)

    PhaseChange *phch;
};

#endif
