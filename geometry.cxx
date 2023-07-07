#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include <cmath>
#include <limits>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "geometry.hpp"
#include "bc.hpp"
#include "mesh.hpp"


/* Given two points, returns the distance^2 */
double dist2(const double* a, const double* b)
{
    double sum = 0;
    for (int i=0; i<NDIMS; ++i) {
        double d = b[i] - a[i];
        sum += d * d;
    }
    return sum;
}


/* Given four 3D points, returns the (signed) volume of the enclosed
   tetrahedron */
#pragma acc routine seq
static double tetrahedron_volume(const double *d0,
                                 const double *d1,
                                 const double *d2,
                                 const double *d3)
{
    double x01 = d0[0] - d1[0];
    double x12 = d1[0] - d2[0];
    double x23 = d2[0] - d3[0];

    double y01 = d0[1] - d1[1];
    double y12 = d1[1] - d2[1];
    double y23 = d2[1] - d3[1];

    double z01 = d0[2] - d1[2];
    double z12 = d1[2] - d2[2];
    double z23 = d2[2] - d3[2];

    return (x01*(y23*z12 - y12*z23) +
            x12*(y01*z23 - y23*z01) +
            x23*(y12*z01 - y01*z12)) / 6;
}


/* Given two points, returns the area of the enclosed triangle */
#pragma acc routine seq
static double triangle_area(const double *a,
                            const double *b,
                            const double *c)
{
    double ab0, ab1, ac0, ac1;

    // ab: vector from a to b
    ab0 = b[0] - a[0];
    ab1 = b[1] - a[1];
    // ac: vector from a to c
    ac0 = c[0] - a[0];
    ac1 = c[1] - a[1];

#ifndef THREED
    // area = norm(cross product of ab and ac) / 2
    return std::fabs(ab0*ac1 - ab1*ac0) / 2;
#else
    double ab2, ac2;
    ab2 = b[2] - a[2];
    ac2 = c[2] - a[2];

    // vector components of ab x ac
    double d0, d1, d2;
    d0 = ab1*ac2 - ab2*ac1;
    d1 = ab2*ac0 - ab0*ac2;
    d2 = ab0*ac1 - ab1*ac0;

    // area = norm(cross product of ab and ac) / 2
    return std::sqrt(d0*d0 + d1*d1 + d2*d2) / 2;
#endif
}

void compute_volume(const double **coord, double &volume)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    const double *a = coord[0];
    const double *b = coord[1];
    const double *c = coord[2];
#ifdef THREED
    const double *d = coord[3];
    volume = tetrahedron_volume(a, b, c, d);
#else
    volume = triangle_area(a, b, c);
#endif
#ifdef USE_NPROF
    nvtxRangePop();
#endif

}

void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    const int bound = volume.size();

    #pragma omp parallel for default(none)      \
        shared(coord, connectivity, volume)
    #pragma acc parallel loop 
    for (int e=0; e<bound; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        const double *d = coord[n3];
        volume[e] = tetrahedron_volume(a, b, c, d);
#else
        volume[e] = triangle_area(a, b, c);
#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void compute_volume(const Variables &var,
                    double_vec &volume)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    const int var_nelem = var.nelem;
    const conn_t &connectivity = *var.connectivity;
    const array_t &coord = *var.coord;


    #pragma omp parallel for default(none)      \
        shared(coord, connectivity, volume)
    #pragma acc parallel loop 
    for (int e=0; e<var_nelem; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        const double *d = coord[n3];
        volume[e] = tetrahedron_volume(a, b, c, d);
#else
        volume[e] = triangle_area(a, b, c);
#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void compute_dvoldt(const Variables &var, double_vec &dvoldt, double_vec &tmp_result_sg)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    /* dvoldt is the volumetric strain rate, weighted by the element volume,
     * lumped onto the nodes.
     */
//    std::fill_n(dvoldt.begin(), var.nnode, 0);

    const int var_nelem = var.nelem;
    const conn_t *var_connectivity = var.connectivity;
    const tensor_t *var_strain_rate= var.strain_rate;
    const double_vec *var_volume = var.volume;
    const double_vec *var_volume_n = var.volume_n;
    const int var_nnode = var.nnode;
    const int_vec2D *var_support = var.support;

    #pragma omp parallel for default(none)      \
        shared(var, tmp_result_sg,var_connectivity,var_strain_rate,var_volume)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {
        const int *conn = (*var_connectivity)[e];
        const double *strain_rate= (*var_strain_rate)[e];
        // TODO: try another definition:
        // dj = (volume[e] - volume_old[e]) / volume_old[e] / dt
        double dj = trace(strain_rate);
        tmp_result_sg[e] = dj * (*var_volume)[e];
    }

    #pragma omp parallel for default(none)      \
        shared(var,dvoldt,tmp_result_sg,var_support,var_volume_n)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        dvoldt[n] = 0.;
        for( auto e = (*var_support)[n].begin(); e < (*var_support)[n].end(); ++e)
	        dvoldt[n] += tmp_result_sg[*e];
        dvoldt[n] /= (*var_volume_n)[n];
    }

    // std::cout << "dvoldt:\n";
    // print(std::cout, dvoldt);
    // std::cout << "\n";
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    /* edvoldt is the averaged (i.e. smoothed) dvoldt on the element.
     * It is used in update_stress() to prevent mesh locking.
     */

    const int var_nelem = var.nelem;
    const conn_t *var_connectivity = var.connectivity;

    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, edvoldt,var_connectivity)
    #pragma acc parallel loop
    for (int e=0; e<var_nelem; ++e) {
        const int *conn = (*var_connectivity)[e];
        double dj = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dj += dvoldt[n];
        }
        edvoldt[e] = dj / NODES_PER_ELEM;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
    // std::cout << "edvoldt:\n";
    // print(std::cout, edvoldt);
    // std::cout << "\n";
}


void NMD_stress(const Param& param, const Variables &var,
    double_vec &dp_nd, tensor_t& stress, double_vec &tmp_result_sg)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // dp_nd is the pressure change, weighted by the element volume,
    // lumped onto the nodes.

//    double **centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements
/*
    // weight with inverse distance
    if(false) {
        #pragma omp parallel for default(none) shared(var,centroid,tmp_result)
        for (int e=0;e<var.nelem;e++) {
            const int *conn = (*var.connectivity)[e];
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                const double *d = (*var.coord)[conn[i]];
                tmp_result[i][e] = 1. / sqrt( dist2(d, centroid[e])  );
                tmp_result[i + NODES_PER_ELEM][e] = tmp_result[i][e] * (*var.dpressure)[e];
            }
        }

        #pragma omp parallel for default(none) shared(var,dp_nd,tmp_result)
        for (int n=0;n<var.nnode;n++) {
            double dist_inv_sum = 0.;
            for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
                const int *conn = (*var.connectivity)[*e];
                for (int i=0;i<NODES_PER_ELEM;i++) {
                    if (n == conn[i]) {
                        dist_inv_sum += tmp_result[ i ][*e];
                        dp_nd[n] += tmp_result[i + NODES_PER_ELEM][*e];
                        break;
                    }
                }
            }
            dp_nd[n] /= dist_inv_sum;
        }

    // weight with volumn
    } else {
        */
    const int var_nelem = var.nelem;
    const conn_t *var_connectivity = var.connectivity;
    const double_vec *var_dpressure = var.dpressure;
    const double_vec *var_volume = var.volume;
    const double_vec *var_volume_n = var.volume_n;
    const int var_nnode = var.nnode;
    const int_vec2D *var_support = var.support;

    #pragma omp parallel for default(none)      \
        shared(var,tmp_result_sg,var_connectivity,var_dpressure,var_volume)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {
        const int *conn = (*var_connectivity)[e];
        double dp = (*var_dpressure)[e];
        tmp_result_sg[e] = dp * (*var_volume)[e];
    }

    #pragma omp parallel for default(none)      \
        shared(var,dp_nd,tmp_result_sg,var_support,var_volume_n)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        dp_nd[n] = 0;
        for( auto e = (*var_support)[n].begin(); e < (*var_support)[n].end(); ++e)
            dp_nd[n] += tmp_result_sg[*e];
        dp_nd[n] /= (*var_volume_n)[n];
    }
//    }


    const int rheol_type = param.mat.rheol_type;
    const double_vec& var_viscosity = *var.viscosity;
    const double ref_visc = param.control.mixed_stress_reference_viscosity;

    /* dp_el is the averaged (i.e. smoothed) dp_nd on the element.
     */
    #pragma omp parallel for default(none)      \
        shared(param, var, dp_nd, stress,var_viscosity,var_connectivity,var_dpressure)
    #pragma acc parallel loop
    for (int e=0; e<var_nelem; ++e) {

        double factor;
        switch (rheol_type) {
        case MatProps::rh_viscous:
        case MatProps::rh_maxwell:
        case MatProps::rh_evp:
            if (var_viscosity[e] < ref_visc)
                factor = 0.;
            else
                factor = std::min(var_viscosity[e] / (ref_visc * 10.), 1.);
            break;
        default:
            factor = 1;
        }

        const int *conn = (*var_connectivity)[e];
        double dp = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dp += dp_nd[n];
        }
        double dp_el = dp / NODES_PER_ELEM;

    	double* s = stress[e];
	    double dp_orig = (*var_dpressure)[e];
        double ddp = ( - dp_orig + dp_el ) / NDIMS * factor;
	    for (int i=0; i<NDIMS; ++i)
            s[i] += ddp;
    }

//    delete [] centroid[0];
//    delete [] centroid;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


double compute_dt(const Param& param, const Variables& var)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    const int nelem = var.nelem;
    const conn_t& connectivity = *var.connectivity;
    const array_t& coord = *var.coord;
    const double_vec& volume = *var.volume;    

    const MatProps *var_mat = var.mat;
    const double var_mat_visc_min= var.mat->visc_min;
    const double var_mat_therm_diff_max = var.mat->therm_diff_max; 
    bool has_thermal_diffusion= param.control.has_thermal_diffusion;


    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();

#ifdef LLVM
    #pragma omp parallel for reduction(min:minl,dt_maxwell,dt_diffusion)    \
        default(none) shared(param, var, nelem, connectivity, coord, volume,var_mat,has_thermal_diffusion)
    #pragma acc parallel loop reduction(min:minl,dt_maxwell,dt_diffusion)
#else
    #pragma omp parallel for reduction(min:minl,dt_maxwell,dt_diffusion)    \
        default(none) shared(param,var, connectivity, coord, volume,var_mat,has_thermal_diffusion)
    #pragma acc parallel loop reduction(min:minl, dt_maxwell, dt_diffusion)
#endif
    for (int e=0; e<nelem; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

        // min height of this element
        double minh;
#ifdef THREED
        {
            int n3 = connectivity[e][3];
            const double *d = coord[n3];

            // max facet area of this tet
            double maxa = std::max(std::max(triangle_area(a, b, c),
                                            triangle_area(a, b, d)),
                                   std::max(triangle_area(c, d, a),
                                            triangle_area(c, d, b)));
            minh = 3 * volume[e] / maxa;
        }
#else
        {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * volume[e] / maxl;
        }
#endif
        dt_maxwell = std::min(dt_maxwell,
                              0.5 * var_mat_visc_min / (1e-40 + var_mat->shearm(e)));
        if (has_thermal_diffusion)
            dt_diffusion = std::min(dt_diffusion,
                                    0.5 * minh * minh / var_mat_therm_diff_max);
        minl = std::min(minl, minh);
    }

    double max_vbc_val;
    if (param.control.characteristic_speed == 0) {
        max_vbc_val = var.max_vbc_val;
        if (param.control.surface_process_option > 0)
            max_vbc_val = std::max(max_vbc_val, var.surfinfo.max_surf_vel);
    }
    else
        max_vbc_val = param.control.characteristic_speed;

    double dt_advection = 0.5 * minl / max_vbc_val;
    double dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (max_vbc_val * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[0] / param.mat.rho0[0]);
    double dt = std::min(std::min(dt_elastic, dt_maxwell),
                         std::min(dt_advection, dt_diffusion)) * param.control.dt_fraction;
    if (param.debug.dt) {
        std::cout << "step #" << var.steps << "  dt: " << dt_maxwell << " " << dt_diffusion
                  << " " << dt_advection << " " << dt_elastic << " sec\n";
    }
    if (dt <= 0) {
        std::cerr << "Error: dt <= 0!  " << dt_maxwell << " " << dt_diffusion
                  << " " << dt_advection << " " << dt_elastic << "\n";
        std::exit(11);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
    return dt;
}


void compute_mass(const Param &param, const Variables &var,
                  double max_vbc_val, double_vec &volume_n,
                  double_vec &mass, double_vec &tmass, elem_cache &tmp_result)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // volume_n is (node-averaged volume * NODES_PER_ELEM)
    // volume_n.assign(volume_n.size(), 0);
    // mass.assign(mass.size(), 0);
    // tmass.assign(tmass.size(), 0);

    const double pseudo_speed = max_vbc_val * param.control.inertial_scaling;

    const bool is_quasi_static = param.control.is_quasi_static;
    const bool has_thermal_diffusion = param.control.has_thermal_diffusion;
    const int var_nelem = var.nelem;
    const double_vec *var_volume = var.volume;
    const MatProps *var_mat = var.mat;
    const int var_nnode = var.nnode;
    const int_vec2D *var_support = var.support;

#ifdef LLVM
    #pragma omp parallel for default(none)      \
        shared(var, param, pseudo_speed, tmp_result)
#else
    #pragma omp parallel for default(none)      \
        shared(var, param, tmp_result,var_mat,var_volume)
#endif
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {
        double *tr = tmp_result[e];
        double rho = (is_quasi_static) ?
            (*var_mat).bulkm(e) / (pseudo_speed * pseudo_speed) :  // pseudo density for quasi-static sim
            (*var_mat).rho(e);                                     // true density for dynamic sim
        double m = rho * (*var_volume)[e] / NODES_PER_ELEM;
        double tm = (*var_mat).rho(e) * (*var_mat).cp(e) * (*var_volume)[e] / NODES_PER_ELEM;
        tr[0] = (*var_volume)[e];
        tr[1] = m;
        if (has_thermal_diffusion)
            tr[2] = tm;
    }

    #pragma omp parallel for default(none)      \
        shared(param,var,volume_n,mass,tmass,tmp_result,var_support)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        volume_n[n]=0;
        mass[n]=0;
        tmass[n]=0;
        for( auto e = (*var_support)[n].begin(); e < (*var_support)[n].end(); ++e) {
            double *tr = tmp_result[*e];
            volume_n[n] += tr[0];
            mass[n] += tr[1];
            if (has_thermal_diffusion)
                tmass[n] += tr[2];
        }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void compute_shape_fn(const Variables &var, shapefn &shpdx, shapefn &shpdy, shapefn &shpdz)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    const int var_nelem = var.nelem;
    const conn_t *var_connectivity = var.connectivity;
    const array_t *var_coord = var.coord;
    const double_vec *var_volume = var.volume;

    #pragma omp parallel for default(none)      \
        shared(var, shpdx, shpdy, shpdz, var_connectivity, var_coord, var_volume)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {

        int n0 = (*var_connectivity)[e][0];
        int n1 = (*var_connectivity)[e][1];
        int n2 = (*var_connectivity)[e][2];

        const double *d0 = (*var_coord)[n0];
        const double *d1 = (*var_coord)[n1];
        const double *d2 = (*var_coord)[n2];

#ifdef THREED
        {
            int n3 = (*var_connectivity)[e][3];
            const double *d3 = (*var_coord)[n3];

            double iv = 1 / (6 * (*var_volume)[e]);

            double x01 = d0[0] - d1[0];
            double x02 = d0[0] - d2[0];
            double x03 = d0[0] - d3[0];
            double x12 = d1[0] - d2[0];
            double x13 = d1[0] - d3[0];
            double x23 = d2[0] - d3[0];

            double y01 = d0[1] - d1[1];
            double y02 = d0[1] - d2[1];
            double y03 = d0[1] - d3[1];
            double y12 = d1[1] - d2[1];
            double y13 = d1[1] - d3[1];
            double y23 = d2[1] - d3[1];

            double z01 = d0[2] - d1[2];
            double z02 = d0[2] - d2[2];
            double z03 = d0[2] - d3[2];
            double z12 = d1[2] - d2[2];
            double z13 = d1[2] - d3[2];
            double z23 = d2[2] - d3[2];

            shpdx[e][0] = iv * (y13*z12 - y12*z13);
            shpdx[e][1] = iv * (y02*z23 - y23*z02);
            shpdx[e][2] = iv * (y13*z03 - y03*z13);
            shpdx[e][3] = iv * (y01*z02 - y02*z01);

            shpdy[e][0] = iv * (z13*x12 - z12*x13);
            shpdy[e][1] = iv * (z02*x23 - z23*x02);
            shpdy[e][2] = iv * (z13*x03 - z03*x13);
            shpdy[e][3] = iv * (z01*x02 - z02*x01);

            shpdz[e][0] = iv * (x13*y12 - x12*y13);
            shpdz[e][1] = iv * (x02*y23 - x23*y02);
            shpdz[e][2] = iv * (x13*y03 - x03*y13);
            shpdz[e][3] = iv * (x01*y02 - x02*y01);
        }
#else
        {
            double iv = 1 / (2 * (*var_volume)[e]);

            shpdx[e][0] = iv * (d1[1] - d2[1]);
            shpdx[e][1] = iv * (d2[1] - d0[1]);
            shpdx[e][2] = iv * (d0[1] - d1[1]);

            shpdz[e][0] = iv * (d2[0] - d1[0]);
            shpdz[e][1] = iv * (d0[0] - d2[0]);
            shpdz[e][2] = iv * (d1[0] - d0[0]);
        }
#endif
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


double elem_quality(const array_t &coord, const conn_t &connectivity,
                    const double_vec &volume, int e)
{
    /* This function returns the quality (0~1) of the element.
     * The quality of an equidistant (i.e. best quality) tetrahedron/triangle is 1.
     */
    double quality;
    double vol = volume[e];
    int n0 = connectivity[e][0];
    int n1 = connectivity[e][1];
    int n2 = connectivity[e][2];

    const double *a = coord[n0];
    const double *b = coord[n1];
    const double *c = coord[n2];

#ifdef THREED
    {
        int n3 = connectivity[e][3];
        const double *d = coord[n3];
        double normalization_factor = 216 * std::sqrt(3);

        double area_sum = (triangle_area(a, b, c) +
                           triangle_area(a, b, d) +
                           triangle_area(c, d, a) +
                           triangle_area(c, d, b));
        quality = normalization_factor * vol * vol / (area_sum * area_sum * area_sum);
    }
#else
    {
        double normalization_factor = 4 * std::sqrt(3);

        double dist2_sum = dist2(a, b) + dist2(b, c) + dist2(a, c);
        quality = normalization_factor * vol / dist2_sum;
    }
#endif

    return quality;
}


double worst_elem_quality(const array_t &coord, const conn_t &connectivity,
                          const double_vec &volume, int &worst_elem)
{
    double q = 1;
    worst_elem = 0;
    for (std::size_t e=0; e<volume.size(); e++) {
        double quality = elem_quality(coord, connectivity, volume, e);
        if (quality < q) {
            q = quality;
            worst_elem = e;
        }
    }
    return q;
}


