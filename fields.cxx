#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "matprops.hpp"
#include "fields.hpp"


void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    var.edvoldt = new double_vec(e);

//    var.marker_in_elem = new int_vec2D(e);

    {
        // these fields are reallocated during remeshing interpolation
        var.temperature = new double_vec(n);
        var.coord0 = new array_t(n);
        var.plstrain = new double_vec(e);
        var.delta_plstrain = new double_vec(e);
        var.vel = new array_t(n, 0);
        var.strain = new tensor_t(e, 0);
        var.stress = new tensor_t(e, 0);
        var.stressyy = new double_vec(e, 0);
    }

    var.ntmp = new double_vec(n);
    if (param.control.is_using_mixed_stress)
        var.dpressure = new double_vec(e);

    var.viscosity = new double_vec(e,param.mat.visc_max);

    var.force = new array_t(n, 0);

    var.strain_rate = new tensor_t(e, 0);

    var.shpdx = new shapefn(e);
    if (NDIMS == 3) var.shpdy = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);

    var.tmp_result = new std::vector<double_vec>(NODES_PER_ELEM*3);
    for (int i =0;i<NODES_PER_ELEM*3;i++)
        (*var.tmp_result)[i].resize(e);
}


void reallocate_variables(const Param& param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    delete var.volume;
    delete var.volume_old;
    delete var.volume_n;
    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    delete var.mass;
    delete var.tmass;
    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    delete var.edvoldt;
    var.edvoldt = new double_vec(e);

//    delete var.marker_in_elem;
//    var.marker_in_elem = new int_vec2D(e);

    delete var.ntmp;
    var.ntmp = new double_vec(n);
    if (param.control.is_using_mixed_stress) {
        delete var.dpressure;
        var.dpressure = new double_vec(e);
    }
    delete var.viscosity;
    var.viscosity = new double_vec(e,param.mat.visc_max);
    delete var.force;
    var.force = new array_t(n, 0);

    delete var.strain_rate;
    var.strain_rate = new tensor_t(e, 0);

    delete var.shpdx;
    delete var.shpdz;
    var.shpdx = new shapefn(e);
    if (NDIMS == 3) {
        delete var.shpdy;
        var.shpdy = new shapefn(e);
    }
    var.shpdz = new shapefn(e);

    delete var.mat;
    var.mat = new MatProps(param, var);

    for (int i =0;i<NODES_PER_ELEM*3;i++)
        (*var.tmp_result)[i].resize(e);

}


void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot, double_vec2D &tmp_result)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    tdot.assign(var.nnode, 0);

    const int var_nelem=var.nelem;
    const conn_t *var_connectivity=var.connectivity;
    const shapefn *var_shpdx=var.shpdx;
    const shapefn *var_shpdz=var.shpdz; 
    const shapefn *var_shpdy=var.shpdy;
    double_vec *var_volume = var.volume;
    MatProps *var_mat = var.mat;

#ifdef USE_NPROF
    nvtxRangePushA("ACC 1");
#endif

    #pragma omp parallel for default(none)      \
        shared(var,temperature,tmp_result)
    #pragma acc kernels
    for (int e=0;e<var_nelem;e++) {
        // diffusion matrix
        double D[NODES_PER_ELEM][NODES_PER_ELEM];

        const int *conn = (*var_connectivity)[e];
        double kv = var_mat->k(e) *  (*var_volume)[e]; // thermal conductivity * volume
        const double *shpdx = (*var_shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var_shpdy)[e];
#endif
        const double *shpdz = (*var_shpdz)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int j=0; j<NODES_PER_ELEM; ++j) {
#ifdef THREED
                D[i][j] = (shpdx[i] * shpdx[j] +
                            shpdy[i] * shpdy[j] +
                            shpdz[i] * shpdz[j]);
#else
                D[i][j] = (shpdx[i] * shpdx[j] +
                            shpdz[i] * shpdz[j]);
#endif
            }
        }

        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double diffusion = 0;
            for (int j=0; j<NODES_PER_ELEM; ++j)
                diffusion += D[i][j] * temperature[conn[j]];
            tmp_result[i][e] = diffusion * kv;
        }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif

    const int var_nnode = var.nnode;
    const int_vec2D *var_support = var.support;
    const uint_vec *var_bcflag = var.bcflag;
    const double surface_temperature = param.bc.surface_temperature;
    const double_vec *var_tmass = var.tmass;
    const double var_dt = var.dt;

#ifdef USE_NPROF
    nvtxRangePushA("ACC 2");
#endif

    #pragma omp parallel for default(none)      \
        shared(param,var,tdot,temperature,tmp_result)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        for( auto e = (*var_support)[n].begin(); e < (*var_support)[n].end(); ++e) {
            const int *conn = (*var_connectivity)[*e];
            for (int i=0;i<NODES_PER_ELEM;i++) {
                if (n == conn[i]) {
                    tdot[n] += tmp_result[i][*e];
                    break;
                }
            }
        }
    // Combining temperature update and bc in the same loop for efficiency,
    // since only the top boundary has Dirichlet bc, and all the other boundaries
    // have no heat flux bc.
        if ((*var_bcflag)[n] & BOUNDZ1)
            temperature[n] = surface_temperature;
        else
            temperature[n] -= tdot[n] * var_dt / (*var_tmass)[n];
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif

    // Combining temperature update and bc in the same loop for efficiency,
    // since only the top boundary has Dirichlet bc, and all the other boundaries
    // have no heat flux bc.
//    #pragma omp parallel for default(none)      \
        shared(var, param, tdot, temperature)
//    for (int n=0; n<var.nnode; ++n) {
//        if ((*var.bcflag)[n] & BOUNDZ1)
//            temperature[n] = param.bc.surface_temperature;
//        else
//            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
//    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void update_strain_rate(const Variables& var, tensor_t& strain_rate)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double *v[NODES_PER_ELEM];

    // extract all variables
    const int var_nelem=var.nelem;
    const conn_t *var_connectivity=var.connectivity;
    const shapefn *var_shpdx=var.shpdx;
    const shapefn *var_shpdz=var.shpdz; 
    const shapefn *var_shpdy=var.shpdy;
    array_t *var_vel=var.vel;

    #pragma omp parallel for default(none)      \
        shared(var, strain_rate) private(v)
    #pragma acc parallel loop private(v)
    for (int e=0; e<var_nelem; ++e) {
        const int *conn = (*var_connectivity)[e];
        const double *shpdx = (*var_shpdx)[e];
        const double *shpdz = (*var_shpdz)[e];
        double *s = strain_rate[e];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = (*var_vel)[conn[i]];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        const double *shpdy = (*var_shpdy)[e];
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
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


static void apply_damping(const Param& param, const Variables& var, array_t& force)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // flatten 2d arrays to simplify indexing
    double* ff = force.data();
    const double* v = var.vel->data();
    const double small_vel = 1e-13;

    const int bound = var.nnode*NDIMS;
    const double damping_factor = param.control.damping_factor;

    switch (param.control.damping_option) {
    case 0:
        // no damping, stress field can become very noisy
        break;
    case 1:
        // damping when force and velocity are parallel
        // acclerating when force and velocity are anti-parallel
#ifdef LLVM
        #pragma omp parallel for default(none)          \
            shared(var, param, ff, v, small_vel)
#else
        #pragma omp parallel for default(none)          \
            shared(var, param, ff, v)
#endif
        #pragma acc kernels
        for (int i=0; i<bound; ++i) {
            if (std::fabs(v[i]) > small_vel) {
                ff[i] -= damping_factor * std::copysign(ff[i], v[i]);
            }
        }
        break;
    case 2:
        // damping prop. to force
        #pragma omp parallel for default(none)          \
            shared(var, param, ff, v)
        #pragma acc kernels
        for (int i=0; i<bound; ++i) {
            ff[i] -= damping_factor * ff[i];
        }
        break;
    case 3:
        // damping when force and velocity are parallel
        // weakly acclerating when force and velocity are anti-parallel
        #pragma omp parallel for default(none)          \
            shared(var, param, ff, v)
        #pragma acc kernels
        for (int i=0; i<bound; ++i) {
            if ((ff[i]<0) == (v[i]<0)) {
                // strong damping
                ff[i] -= damping_factor * ff[i], v[i];
            }
            else {
                // weak acceleration
                ff[i] += (1 - damping_factor) * ff[i];
            }
        }
        break;
    case 4:
        // rayleigh damping
        break;
    default:
        std::cerr << "Error: unknown damping_option: " << param.control.damping_option << '\n';
        std::exit(1);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

/*
#pragma acc routine seq
static double rho(const conn_t &var_connectivity, \
    const double_vec &var_temperature, const int_vec2D &var_elemmarkers, \
    const double_vec &rho0, const double_vec &alpha, int nmat,  int e)
{
    const double celsius0 = 273;

    // average temperature of this element
    double T = 0;
    const int *conn = var_connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += var_temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    double TinCelsius = T - celsius0;
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        // TODO: compressibility
        result += rho0[m] * (1 - alpha[m] * TinCelsius) * var_elemmarkers[e][m];
        n += var_elemmarkers[e][m];
    }
    return result / n;
}
*/

void update_force(const Param& param, const Variables& var, array_t& force, double_vec2D& tmp_result)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    std::fill_n(force.data(), var.nnode*NDIMS, 0);


    const int var_nelem = var.nelem;
    const conn_t *var_connectivity = var.connectivity;
    const double_vec *var_temperature = var.temperature;
    const int_vec2D *var_elemmarkers = var.elemmarkers;
    const shapefn *var_shpdx = var.shpdx;
#ifdef THREED
    const shapefn *var_shpdy = var.shpdy;
#endif
    const shapefn *var_shpdz = var.shpdz;
    tensor_t *var_stress = var.stress;
    double_vec *var_volume = var.volume;
    const MatProps *var_mat = var.mat;
    const double gravity = param.control.gravity;
    const int nmat = param.mat.nmat;
    const double_vec *rho0 = &param.mat.rho0;
    const double_vec *alpha = &param.mat.alpha;

#ifdef USE_NPROF
    nvtxRangePushA("ACC 1");
#endif

    #pragma omp parallel for default(none)      \
        shared(var,param,tmp_result)
    #pragma acc kernels
    for (int e=0;e<var_nelem;e++) {

        const int *conn = (*var_connectivity)[e];
        const double *shpdx = (*var_shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var_shpdy)[e];
#endif
        const double *shpdz = (*var_shpdz)[e];
        double *s = (*var_stress)[e];
        double vol = (*var_volume)[e];

        double buoy = 0;
        if (gravity != 0)
            buoy = var_mat->rho(e) * gravity / NODES_PER_ELEM;

        for (int i=0; i<NODES_PER_ELEM; ++i) {
#ifdef THREED
            tmp_result[i][e] = (s[0]*shpdx[i] + s[3]*shpdy[i] + s[4]*shpdz[i]) * vol;
            tmp_result[i+NODES_PER_ELEM][e] = (s[3]*shpdx[i] + s[1]*shpdy[i] + s[5]*shpdz[i]) * vol;
            tmp_result[i+NODES_PER_ELEM*2][e] = (s[4]*shpdx[i] + s[5]*shpdy[i] + s[2]*shpdz[i] + buoy) * vol;
#else
            tmp_result[i][e] = (s[0]*shpdx[i] + s[2]*shpdz[i]) * vol;
            tmp_result[i+NODES_PER_ELEM][e] = (s[2]*shpdx[i] + s[1]*shpdz[i] + buoy) * vol;
#endif
        }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif

    const int var_nnode = var.nnode;
    const int_vec2D *var_support = var.support;

#ifdef USE_NPROF
    nvtxRangePushA("ACC 2");
#endif

    #pragma omp parallel for default(none)      \
        shared(var,force,tmp_result)
    #pragma acc kernels
    for (int n=0;n<var_nnode;n++) {
        double *f = force[n];
        for( auto e = (*var_support)[n].begin(); e < (*var_support)[n].end(); ++e) {
            const int *conn = (*var_connectivity)[*e];
            for (int i=0;i<NODES_PER_ELEM;i++) {
                if (n == conn[i]) {
                    for (int j=0;j<NDIMS;j++)
                        f[j] -= tmp_result[i+NODES_PER_ELEM*j][*e];
                    break;
                }
            }
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif

    apply_stress_bcs(param, var, force);

    if (param.control.is_quasi_static) {
        apply_damping(param, var, force);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void update_velocity(const Variables& var, array_t& vel)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    const double* m = &(*var.mass)[0];
    // flatten 2d arrays to simplify indexing
    const double* f = var.force->data();
    double* v = vel.data();

    const int var_nnode = var.nnode;
    const double var_dt = var.dt;

    #pragma omp parallel for default(none) \
        shared(var, m, f, v)
    #pragma acc kernels
    for (int i=0; i<var_nnode*NDIMS; ++i) {
        int n = i / NDIMS;
        v[i] += var_dt * f[i] / m[n];
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void update_coordinate(const Variables& var, array_t& coord)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double* x = var.coord->data();
    const double* v = var.vel->data();

    // for gpu parallelization dt and bound need to be sent to 
    const int bound=var.nnode*NDIMS;
    const double var_dt=var.dt;

    #pragma omp parallel for default(none) \
        shared(var, x, v)
    #pragma acc parallel loop
    for (int i=0; i<bound; ++i) {
        x[i] += v[i] * var_dt;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

#ifdef THREED
    #pragma acc routine seq 
    void jaumann_rate_3d(double *s, double dt, double w3, double w4, double w5)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[3] * w3 - 2.0 * s[4] * w4;
        s_inc[1] = 2.0 * s[3] * w3 - 2.0 * s[5] * w5;
        s_inc[2] = 2.0 * s[4] * w4 + 2.0 * s[5] * w5;
        s_inc[3] = s[0] * w3 - s[1] * w3 - s[4] * w5 - s[5] * w4;
        s_inc[4] = s[0] * w4 - s[2] * w4 + s[3] * w5 - s[5] * w3;
        s_inc[5] = s[1] * w5 - s[2] * w5 + s[3] * w4 + s[4] * w3;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

#else
    #pragma acc routine seq
    void jaumann_rate_2d(double *s, double dt, double w2)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[2] * w2;
        s_inc[1] = 2.0 * s[2] * w2;
        s_inc[2] = s[0] * w2 - s[1] * w2;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

#endif

}


void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // The spin rate tensor, W, and the Cauchy stress tensor, S, are
    //     [  0  w3  w4]     [s0 s3 s4]
    // W = [-w3   0  w5],  S=[s3 s1 s5].
    //     [-w4 -w5   0]     [s4 s5 s2]
    //
    // Stress (and strain) increment based on the Jaumann rate is
    // dt*(S*W-W*S). S*W-W*S is also symmetric.
    // So, following the indexing of stress tensor,
    // we get
    // sj[0] = dt * ( -2 * s3 * w3 - 2 * s4 * w4)
    // sj[1] = dt * (  2 * s3 * w3 - 2 * s5 * w5)
    // sj[2] = dt * (  2 * s4 * w4 + 2 * s5 * w5)
    // sj[3] = dt * ( s0 * w3 - s1 * w3 - s4 * w5 - s5 * w4)
    // sj[4] = dt * ( s0 * w4 - s2 * w4 + s3 * w5 - s5 * w3)
    // sj[5] = dt * ( s1 * w5 - s2 * w5 + s3 * w4 + s4 * w3)

    const int var_nelem=var.nelem;
    const double var_dt=var.dt;

    const shapefn *var_shpdx=var.shpdx;
    const shapefn *var_shpdz=var.shpdz;
#ifdef THREED
    const shapefn *var_shpdy=var.shpdy;
#endif
    const conn_t *var_connectivity=var.connectivity;
    array_t *var_vel=var.vel;

    #pragma omp parallel for default(none) \
        shared(var, stress, strain)
    #pragma acc parallel loop
    for (int e=0; e<var_nelem; ++e) {
        const int *conn = (*var_connectivity)[e];

#ifdef THREED

        double w3, w4, w5;
        {
            const double *shpdx = (*var_shpdx)[e];
            const double *shpdy = (*var_shpdy)[e];
            const double *shpdz = (*var_shpdz)[e];

            double *v[NODES_PER_ELEM];
            for (int i=0; i<NODES_PER_ELEM; ++i)
                v[i] = (*var_vel)[conn[i]];

            w3 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w3 += 0.5 * (v[i][0] * shpdy[i] - v[i][1] * shpdx[i]);

            w4 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w4 += 0.5 * (v[i][0] * shpdz[i] - v[i][NDIMS-1] * shpdx[i]);

            w5 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w5 += 0.5 * (v[i][1] * shpdz[i] - v[i][NDIMS-1] * shpdy[i]);
        }

        jaumann_rate_3d(stress[e], var_dt, w3, w4, w5);
        jaumann_rate_3d(strain[e], var_dt, w3, w4, w5);

#else

        double w2;
        {
            const double *shpdx = (*var_shpdx)[e];
            const double *shpdz = (*var_shpdz)[e];

            double *v[NODES_PER_ELEM];
            for (int i=0; i<NODES_PER_ELEM; ++i)
                v[i] = (*var_vel)[conn[i]];

            w2 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w2 += 0.5 * (v[i][NDIMS-1] * shpdx[i] - v[i][0] * shpdz[i]);
        }

        jaumann_rate_2d(stress[e], var_dt, w2);
        jaumann_rate_2d(strain[e], var_dt, w2);

#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

