#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static void allocate_variables(Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    var.jacobian = new double_vec(n);
    var.ejacobian = new double_vec(e);

    var.temperature = new double_vec(n);
    var.plstrain = new double_vec(e);
    var.tmp0 = new double_vec(std::max(n,e));

    var.vel = new double2d(boost::extents[n][NDIMS]);
    var.force = new double2d(boost::extents[n][NDIMS]);

    var.strain_rate = new double2d(boost::extents[e][NSTR]);
    var.strain = new double2d(boost::extents[e][NSTR]);
    var.stress = new double2d(boost::extents[e][NSTR]);

    var.shpdx = new double2d(boost::extents[e][NODES_PER_ELEM]);
    if (NDIMS == 3) var.shpdy = new double2d(boost::extents[e][NODES_PER_ELEM]);
    var.shpdz = new double2d(boost::extents[e][NODES_PER_ELEM]);
}


static void create_matprops(const Param &par, Variables &var)
{
    // TODO: get material properties from cfg file
    var.mat = new MatProps(1, MatProps::rh_maxwell);
}


void initial_stress_state(const Param &param, const Variables &var,
                          double2d &stress, double2d &strain,
                          double &compensation_pressure)
{
    if (param.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    // XXX: compute reference pressure correctly
    const double rho = var.mat->density(0);
    const double ks = var.mat->bulkm(0);
    compensation_pressure = rho * param.gravity * param.mesh.zlength;
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = - rho * param.gravity * zcenter;
            strain[e][i] = - rho * param.gravity * zcenter / ks / NDIMS;
        }
    }
}


void initial_temperature(const Param &param, const Variables &var, double_vec &temperature)
{
    const double oceanic_plate_age = 1e6 * YEAR2SEC;
    const double diffusivity = 1e-6;

    for (int i=0; i<var.nnode; ++i) {
        double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * oceanic_plate_age);
        temperature[i] = param.bc.surface_temperature +
            (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
    }
}


void apply_vbcs(const Param &param, const Variables &var, double2d &vel)
{
    // TODO: adding different types of vbcs later

    // diverging x-boundary
    for (int i=0; i<var.nnode; ++i) {
        int flag = (*var.bcflag)[i];

        // X
        if (flag & BOUNDX0) {
            vel[i][0] = -param.bc.max_vbc_val;
        }
        else if (flag & BOUNDX1) {
            vel[i][0] = param.bc.max_vbc_val;
        }
#ifdef THREED
        // Y
        if (flag & BOUNDY0) {
            vel[i][1] = 0;
        }
        else if (flag & BOUNDY1) {
            vel[i][1] = 0;
        }
#endif
        // Z
        if (flag & BOUNDZ0) {
            //vel[i][NDIMS-1] = 0;
        }
        else if (flag & BOUNDZ1) {
            vel[i][NDIMS-1] = 0;
        }
    }
}


void init(const Param& param, Variables& var)
{
    void create_matprops(const Param&, Variables&);

    create_new_mesh(param, var);
    allocate_variables(var);
    create_matprops(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.volume, *var.volume_n);
    *var.volume_old = *var.volume;
    compute_mass(param, *var.coord, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume,
                     *var.shpdx, *var.shpdy, *var.shpdz);
    // XXX
    //create_jacobian();

    initial_stress_state(param, var, *var.stress, *var.strain, var.compensation_pressure);
    initial_temperature(param, var, *var.temperature);
    apply_vbcs(param, var, *var.vel);
}


void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{
    // diffusion matrix
    double D[NODES_PER_ELEM][NODES_PER_ELEM];

    tdot.assign(var.nnode, 0);
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volumn
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int j=0; j<NODES_PER_ELEM; ++j) {
                if (NDIMS == 3) {
                    D[i][j] = ((*var.shpdx)[e][i] * (*var.shpdx)[e][j] +
                               (*var.shpdy)[e][i] * (*var.shpdy)[e][j] +
                               (*var.shpdz)[e][i] * (*var.shpdz)[e][j]);
                }
                else {
                    D[i][j] = ((*var.shpdx)[e][i] * (*var.shpdx)[e][j] +
                               (*var.shpdz)[e][i] * (*var.shpdz)[e][j]);
                }
            }
        }
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double diffusion = 0;
            for (int j=0; j<NODES_PER_ELEM; ++j)
                diffusion += D[i][j] * temperature[conn[j]];

            tdot[conn[i]] += diffusion * kv;
        }
    }

    for (int n=0; n<var.nnode; ++n) {
        if ((*var.bcflag)[n] & BOUNDZ1)
            temperature[n] = param.bc.surface_temperature;
        else
            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
    }
}


void update_strain_rate(const Variables& var, double2d& strain_rate)
{
    double *v[NODES_PER_ELEM];

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        const double *shpdx = &(*var.shpdx)[e][0];
        const double *shpdy = &(*var.shpdy)[e][0];
        const double *shpdz = &(*var.shpdz)[e][0];
        double *s = &(*var.strain_rate)[e][0];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = &(*var.vel)[conn[i]][0];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        // YY component
        n = 1;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][1] * shpdy[i];
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][NDIMS-1] * shpdz[i];

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] + v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] + v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] + v[i][2] * shpdy[i]);
#endif
    }
}


void update_force() {};
void rotate_stress() {};


void update_velocity(const Variables& var, double2d& vel)
{
    const double* m = &(*var.volume)[0];
    // flatten 2d arrays to simplify indexing
    const double* f = var.force->data();
    double* v = vel.data();
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        int n = i / NDIMS;
        v[i] += var.dt * f[i] / m[n];
    }
}


static void update_coordinate(const Variables& var, double2d_ref& coord)
{
    double* x = var.coord->data();
    const double* v = var.vel->data();
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        x[i] += v[i] * var.dt;
    }

    // surface_processes()
}


void update_mesh(const Param& param, Variables& var)
{
    update_coordinate(var, *var.coord);

    var.volume->swap(*var.volume_old);
    compute_volume(*var.coord, *var.connectivity, *var.volume, *var.volume_n);
    compute_mass(param, *var.coord, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume,
                     *var.shpdx, *var.shpdy, *var.shpdz);
}


int main(int argc, const char* argv[])
{
    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        return -1;
    }

    Param param;
    void get_input_parameters(const char*, Param&);
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    var.time = 0;
    var.steps = 0;
    var.frame = 0;

    if (! param.sim.is_restarting) {
        init(param, var);
        output(param, var);
        var.frame ++;
    }
    else {
        restart();
        var.frame ++;
    }

    var.dt = compute_dt(param, var);

    do {
        var.steps ++;
        var.time += var.dt;

        update_temperature(param, var, *var.temperature, *var.tmp0);
        update_strain_rate(var, *var.strain_rate);
        update_stress(var, *var.stress, *var.strain, *var.plstrain);
        update_force();
        update_velocity(var, *var.vel);
        apply_vbcs(param, var, *var.vel);
        update_mesh(param, var);
        // dt computation is expensive, and dt only changes slowly
        // don't have to do it every time step
        if (var.steps % 10 == 0) var.dt = compute_dt(param, var);
        rotate_stress();

        if ( (var.steps == var.frame * param.sim.output_step_interval) ||
             (var.time > var.frame * param.sim.output_time_interval_in_yr * YEAR2SEC) ) {
            output(param, var);
            var.frame ++;
        }

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    return 0;
}
