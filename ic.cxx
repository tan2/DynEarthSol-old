#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "ic-read-temp.hpp"
#include "ic.hpp"

namespace {

    class Zone
    {
    public:
        virtual ~Zone() {};
        virtual bool contains(const double x[NDIMS]) const = 0;
    };

    class Empty_zone : public Zone
    {
    public:
        bool contains(const double x[NDIMS]) const {return false;}
    };


    class Planar_zone : public Zone
    {
    private:
        const double az, incl;
        const double halfwidth; // in meter
#ifdef THREED
        const double ymin, ymax; // in meter
#endif
        const double zmin, zmax; // in meter
        const double *x0;

    public:
        Planar_zone(const double center[NDIMS], double azimuth, double inclination, double halfwidth_,
#ifdef THREED
                    double ymin_, double ymax_,
#endif
                    double zmin_, double zmax_) :
            az(std::tan(azimuth * DEG2RAD)), incl(1/std::tan(inclination * DEG2RAD)), halfwidth(halfwidth_),
#ifdef THREED
            ymin(ymin_), ymax(ymax_),
#endif
            zmin(zmin_), zmax(zmax_),
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {}

        bool contains(const double x[NDIMS]) const
        {
            // Is x within halfwidth distance to a plane cutting through x0?
            return (x[NDIMS-1] > zmin &&
                    x[NDIMS-1] < zmax &&
#ifdef THREED
                    x[1] > ymin &&
                    x[1] < ymax &&
#endif
                    std::fabs( (x[0] - x0[0])
#ifdef THREED
                               - az * (x[1] - x0[1])
#endif
                               + incl * (x[NDIMS-1] - x0[NDIMS-1]) ) < halfwidth );
        }
    };


    class Ellipsoidal_zone : public Zone
    {
    private:
        const double *x0;
        double semi_axis2[NDIMS];

    public:
        Ellipsoidal_zone(const double center[NDIMS], const double semi_axis[NDIMS]) :
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {
            for(int i=0; i<NDIMS; i++)
                semi_axis2[i] =  semi_axis[i] * semi_axis[i];
        }

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])/semi_axis2[0]
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])/semi_axis2[1]
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])/semi_axis2[NDIMS-1] < 1 );
        }
    };


    class Gaussian_distribution_point_zone : public Zone
    {
    private:
        const double *x0;
        const double standard_deviation; // in meter

    public:
        Gaussian_distribution_point_zone(const double center[NDIMS], double standard_deviation_) :
            x0(center), // Copy the pointer only, not the data. The caller needs to keep center alive.
            standard_deviation(standard_deviation_)
        {}

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1]) < (standard_deviation * standard_deviation * 16.) );
        }
    };

    class Value
    {
    public:
        virtual ~Value() {};
        virtual double contains(const double x[NDIMS]) const = 0;
    };

    class Constant_value : public Value
    {
    public:
        double contains(const double x[NDIMS]) const {return 1.;}
    };

    class Gaussian_distribution_point_value : public Value
    {
    private:
        const double *x0;
        const double standard_deviation; // in meter

    public:
        Gaussian_distribution_point_value(const double center[NDIMS], double standard_deviation_) :
            x0(center), // Copy the pointer only, not the data. The caller needs to keep center alive.
            standard_deviation(standard_deviation_)
        {}

        double contains(const double x[NDIMS]) const
        {
            return exp(-( ((x[0] - x0[0])*(x[0] - x0[0])
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])) / (2.*standard_deviation*standard_deviation)));
        }
    };


} // anonymous namespace


void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, tensor_t &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        double p = ref_pressure(param, zcenter);
        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
        }

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = -p;
            strain[e][i] = -p / ks / NDIMS;
        }
        if (param.mat.is_plane_strain)
            stressyy[e] = -p;
    }

    compensation_pressure = ref_pressure(param, -param.mesh.zlength);
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    Zone *weakzone;
    Value *weakvalue;

    // TODO: adding different types of weak zone
    double weakzone_center[NDIMS]; // this variable must outlive weakzone
    switch (param.ic.weakzone_option) {
    case 0:
        weakzone = new Empty_zone();
        weakvalue = new Constant_value();
        break;
    case 1:
        // a planar weak zone, cut through top center
        weakzone_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        weakzone_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        weakzone_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Planar_zone(weakzone_center,
                                   param.ic.weakzone_azimuth,
                                   param.ic.weakzone_inclination,
                                   param.ic.weakzone_halfwidth * param.mesh.resolution,
#ifdef THREED
                                   param.ic.weakzone_y_min * param.mesh.ylength,
                                   param.ic.weakzone_y_max * param.mesh.ylength,
#endif
                                   -param.ic.weakzone_depth_max * param.mesh.zlength,
                                   -param.ic.weakzone_depth_min * param.mesh.zlength);
        weakvalue = new Constant_value();
        break;
    case 2:
        // a ellipsoidal weak zone
        double semi_axis[NDIMS];
        weakzone_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
        semi_axis[0] = param.ic.weakzone_xsemi_axis;
#ifdef THREED
        weakzone_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
        semi_axis[1] = param.ic.weakzone_ysemi_axis;
#endif
        weakzone_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        semi_axis[NDIMS-1] = param.ic.weakzone_zsemi_axis;
        weakzone = new Ellipsoidal_zone(weakzone_center, semi_axis);
        weakvalue = new Constant_value();
        break;
    case 3:
        // a Gaussian distribution point weak zone
        weakzone_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        weakzone_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        weakzone_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Gaussian_distribution_point_zone(weakzone_center, param.ic.weakzone_standard_deviation);
        weakvalue = new Gaussian_distribution_point_value(weakzone_center, param.ic.weakzone_standard_deviation);
        break;
    default:
        std::cerr << "Error: unknown weakzone_option: " << param.ic.weakzone_option << '\n';
        std::exit(1);
    }

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[NDIMS] = {0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }

        if (weakzone->contains(center))
            plstrain[e] = param.ic.weakzone_plstrain * weakvalue->contains(center);
        // Find the most abundant marker mattype in this element
        // int_vec &a = (*var.elemmarkers)[e];
        // int material = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    }

    delete weakzone;
    delete weakvalue;
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    switch(param.ic.temperature_option) {
    case 0:
        {
            const double age = param.ic.oceanic_plate_age_in_yr * YEAR2SEC;
            const MatProps &mat = *var.mat;
            const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element

            for (int i=0; i<var.nnode; ++i) {
                double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * age);
                temperature[i] = param.bc.surface_temperature +
                    (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
            }
            break;
        }
    case 1:
        {
            // Continental geotherm
            const double pi = 3.14159265358979323846;

            const int dens_c = param.mat.rho0[param.mat.mattype_crust];
            const int dens_m = param.mat.rho0[param.mat.mattype_mantle];
            const double cond_c = param.mat.therm_cond[std::min(int(param.mat.therm_cond.size())-1, param.mat.mattype_crust)];
            const double cond_m = param.mat.therm_cond[std::min(int(param.mat.therm_cond.size())-1, param.mat.mattype_mantle)];
            const double diff_m = cond_m/1000./dens_m;

            const double age = param.ic.continental_plate_age_in_yr * YEAR2SEC;
            const double hs = param.ic.radiogenic_heating_of_crust; 
            const double hr = param.ic.radiogenic_folding_depth;
            const double hc = param.ic.radiogenic_crustal_thickness;
            const double hl = param.ic.lithospheric_thickness;

            const double t_top = param.bc.surface_temperature;
            const double t_bot = param.bc.mantle_temperature;

            const double tr = dens_c * hs * hr*hr / cond_c * exp(1.-exp(-hc/hr));
            const double q_m = (t_bot - t_top - tr) / (hc / cond_c+(hl-hc) / cond_m);
            const double tm  = t_top + (q_m/cond_c) * hc + tr;
            const double tau_d = hl*hl / (pi*pi*diff_m);

            for (int i=0; i<var.nnode; ++i) {
                double y = -(*var.coord)[i][NDIMS-1];
                double tss;
                // steady state part
                if (y <= hc)
                    tss = t_top + (q_m/cond_c)*y + (dens_c*hs*hr*hr/cond_c) * exp(1.-exp(-y/hr));
                else 
                    tss = tm + (q_m/cond_m) * (y - hc);
                // time-dependent part
                double tt = 0.;
                double pp = -1.;
                double an;
                for (int k=1;k<101;k++) {
                    an = 1.*k;
                    pp = -pp;
                    tt = tt +pp/(an)*exp(-an*an*age/tau_d)*sin(pi*k*(hl-y)/hl);
                }

                temperature[i] = tss + 2./pi*(t_bot-t_top)*tt;

                if (temperature[i] > t_bot || y >= hl)
                    temperature[i] = t_bot;

                if (y == 0.)
                    temperature[i] = t_top;
            }
            break;
        }
    case 90:
        read_external_temperature_from_comsol(param, var, *var.temperature);
        break;
    default:
        std::cout << "Error: unknown ic.temperature option: " << param.ic.temperature_option << '\n';
        std::exit(1);
    }
}


