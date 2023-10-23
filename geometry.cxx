#include <cmath>
#include <limits>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "geometry.hpp"


/* Given two points, returns the distance^2 */
double dist2(const double* a, const double* b)
{
    double sum = 0;;
    for (int i=0; i<NDIMS; ++i) {
        double d = b[i] - a[i];
        sum += d * d;
    }
    return sum;
}


/* Given four 3D points, returns the (signed) volume of the enclosed
   tetrahedron */
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


void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume)
{
    #pragma omp parallel for default(none)      \
        shared(coord, connectivity, volume)
    for (std::size_t e=0; e<volume.size(); ++e) {
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
}


void compute_dvoldt(const Variables &var, double_vec &dvoldt, double_vec &tmp_result_sg)
{
    /* dvoldt is the volumetric strain rate, weighted by the element volume,
     * lumped onto the nodes.
     */
    const double_vec& volume = *var.volume;
    const double_vec& volume_n = *var.volume_n;
    const conn_t& connectivity = *var.connectivity;
    const tensor_t& strain_rate = *var.strain_rate;
    const int_vec2D& support = *var.support;
    int var_nelem = var.nelem;
    int var_nnode = var.nnode;

    #pragma omp parallel for default(none)      \
        shared(tmp_result_sg,strain_rate,volume,var_nelem)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {
        const double *srate= strain_rate[e];
        // TODO: try another definition:
        // dj = (volume[e] - volume_old[e]) / volume_old[e] / dt
        double dj = trace(srate);
        tmp_result_sg[e] = dj * volume[e];
    }

    #pragma omp parallel for default(none)      \
        shared(dvoldt,tmp_result_sg,support,volume_n,var_nnode)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        dvoldt[n] = 0.;
        for( auto e = support[n].begin(); e < support[n].end(); ++e)
	        dvoldt[n] += tmp_result_sg[*e];
        dvoldt[n] /= volume_n[n];
    }

    // std::cout << "dvoldt:\n";
    // print(std::cout, dvoldt);
    // std::cout << "\n";
}


void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt)
{
    /* edvoldt is the averaged (i.e. smoothed) dvoldt on the element.
     * It is used in update_stress() to prevent mesh locking.
     */
    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, edvoldt)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double dj = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dj += dvoldt[n];
        }
        edvoldt[e] = dj / NODES_PER_ELEM;
    }

    // std::cout << "edvoldt:\n";
    // print(std::cout, edvoldt);
    // std::cout << "\n";
}


void NMD_stress(const Variables &var, double_vec &dp_nd, tensor_t& stress,
                double_vec &tmp_result_sg)
{
    /* dp_nd is the pressure change, weighted by the element volume,
     * lumped onto the nodes.
     */
    const double_vec& volume = *var.volume;
    const double_vec& volume_n = *var.volume_n;
    const conn_t& connectivity = *var.connectivity;
    const double_vec& dpressure = *var.dpressure;
    const int_vec2D& support = *var.support;
    // const double_vec& viscosity = *var.viscosity;
    int var_nnode = var.nnode;
    int var_nelem = var.nelem;
    // const int rheol_type = param.mat.rheol_type;
    // const double ref_visc = param.control.mixed_stress_reference_viscosity;

    #pragma omp parallel for default(none)      \
        shared(tmp_result_sg,connectivity,dpressure,volume,var_nelem)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {
        const int *conn = connectivity[e];
        double dp = dpressure[e];
        tmp_result_sg[e] = dp * volume[e];
    }

    #pragma omp parallel for default(none)      \
        shared(dp_nd,tmp_result_sg,support,volume_n,var_nnode)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        dp_nd[n] = 0;
        for( auto e = support[n].begin(); e < support[n].end(); ++e)
            dp_nd[n] += tmp_result_sg[*e];
        dp_nd[n] /= volume_n[n];
    }



    /* dp_el is the averaged (i.e. smoothed) dp_nd on the element.
     */
    #pragma omp parallel for default(none)      \
        shared(connectivity,dp_nd,stress,dpressure,var_nelem)
        // shared(param, dp_nd, stress,viscosity,connectivity,dpressure)
    #pragma acc parallel loop
    for (int e=0; e<var_nelem; ++e) {

        // double factor;
        // switch (rheol_type) {
        // case MatProps::rh_viscous:
        // case MatProps::rh_maxwell:
        // case MatProps::rh_evp:
        //     if (viscosity[e] < ref_visc)
        //         factor = 0.;
        //     else
        //         factor = std::min(viscosity[e] / (ref_visc * 10.), 1.);
        //     break;
        // default:
        //     factor = 1;
        // }

        const int *conn = connectivity[e];
        double dp = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dp += dp_nd[n];
        }
        double dp_el = dp / NODES_PER_ELEM;

    	double* s = stress[e];
	    double dp_orig = dpressure[e];
        double ddp = ( - dp_orig + dp_el ) / NDIMS;// * factor;
	    for (int i=0; i<NDIMS; ++i)
            s[i] += ddp;
    }
}


double compute_dt(const Param& param, const Variables& var)
{
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    const int nelem = var.nelem;
    const conn_t& connectivity = *var.connectivity;
    const array_t& coord = *var.coord;
    const double_vec& volume = *var.volume;

    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();

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
                              0.5 * var.mat->visc_min / (1e-40 + var.mat->shearm(e)));
        if (param.control.has_thermal_diffusion)
            dt_diffusion = std::min(dt_diffusion,
                                    0.5 * minh * minh / var.mat->therm_diff_max);
	minl = std::min(minl, minh);
    }

    double dt_advection = 0.5 * minl / var.max_vbc_val;
    double dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (var.max_vbc_val * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(var.mat->bulkm(0) / var.mat->rho(0));
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
    return dt;
}


void compute_mass(const Param &param, const Variables& var,
                  const conn_t &connectivity,
                  const double_vec &volume, const MatProps &mat,
                  double max_vbc_val, double_vec &volume_n,
                  double_vec &mass, double_vec &tmass, elem_cache &tmp_result,
                  int_vec2D &support)
{
    // volume_n is (node-averaged volume * NODES_PER_ELEM)

    double pseudo_speed = max_vbc_val * param.control.inertial_scaling;

    bool is_quasi_static = param.control.is_quasi_static;
    bool has_thermal_diffusion = param.control.has_thermal_diffusion;
    int var_nelem = var.nelem;
    int var_nnode = var.nnode;

    #pragma omp parallel for default(none)      \
        shared(tmp_result,mat,volume,var_nelem,pseudo_speed,is_quasi_static, \
               has_thermal_diffusion)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {
        double *tr = tmp_result[e];
        double rho = (is_quasi_static) ?
            mat.bulkm(e) / (pseudo_speed * pseudo_speed) :  // pseudo density for quasi-static sim
            mat.rho(e);                                     // true density for dynamic sim
        double m = rho * volume[e] / NODES_PER_ELEM;
        double tm = mat.rho(e) * mat.cp(e) * volume[e] / NODES_PER_ELEM;
        tr[0] = volume[e];
        tr[1] = m;
        if (has_thermal_diffusion)
            tr[2] = tm;
    }

    #pragma omp parallel for default(none)      \
        shared(volume_n,mass,tmass,tmp_result,support,var_nnode,has_thermal_diffusion)
    #pragma acc parallel loop
    for (int n=0;n<var_nnode;n++) {
        volume_n[n]=0;
        mass[n]=0;
        tmass[n]=0;
        for( auto e = support[n].begin(); e < support[n].end(); ++e) {
            double *tr = tmp_result[*e];
            volume_n[n] += tr[0];
            mass[n] += tr[1];
            if (has_thermal_diffusion)
                tmass[n] += tr[2];
        }
    }
}


void compute_shape_fn(const Variables& var, const array_t &coord, const conn_t &connectivity,
                      const double_vec &volume,
                      shapefn &shpdx, shapefn &shpdy, shapefn &shpdz)
{
    int var_nelem = var.nelem;

    #pragma omp parallel for default(none)      \
        shared(shpdx, shpdy, shpdz, connectivity, coord, volume, var_nelem)
    #pragma acc parallel loop
    for (int e=0;e<var_nelem;e++) {

        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *d0 = coord[n0];
        const double *d1 = coord[n1];
        const double *d2 = coord[n2];

#ifdef THREED
        {
            int n3 = connectivity[e][3];
            const double *d3 = coord[n3];

            double iv = 1 / (6 * volume[e]);

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
            double iv = 1 / (2 * volume[e]);

            shpdx[e][0] = iv * (d1[1] - d2[1]);
            shpdx[e][1] = iv * (d2[1] - d0[1]);
            shpdx[e][2] = iv * (d0[1] - d1[1]);

            shpdz[e][0] = iv * (d2[0] - d1[0]);
            shpdz[e][1] = iv * (d0[0] - d2[0]);
            shpdz[e][2] = iv * (d1[0] - d0[0]);
        }
#endif
    }
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


