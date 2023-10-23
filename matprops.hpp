#ifndef DYNEARTHSOL3D_MATPROPS_HPP
#define DYNEARTHSOL3D_MATPROPS_HPP

#include "parameters.hpp"

#ifdef ACC

#pragma acc routine seq
double ref_pressure(const Param& param, double z);


using VectorBase = double_vec;

class MatProps
{
public:
    MatProps(const Param& param, const Variables& var);
    ~MatProps();

    const int rheol_type;
    const int nmat;

    #pragma acc routine seq
    double bulkm(int e) const;
    #pragma acc routine seq
    double shearm(int e) const;
    #pragma acc routine seq
    double visc(int e) const;

    #pragma acc routine seq
    double rho(int e) const;
    #pragma acc routine seq
    double cp(int e) const;
    #pragma acc routine seq
    double k(int e) const;

    #pragma acc routine seq
    void plastic_props(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max) const;

    const bool is_plane_strain;
    const double visc_min;
    const double visc_max;
    const double tension_max;
    const double therm_diff_max;

    const static int rh_elastic = 1 << 0;
    const static int rh_viscous = 1 << 1;
    const static int rh_plastic = 1 << 2;
    const static int rh_plastic2d = rh_plastic | 1 << 3;
    const static int rh_maxwell = rh_elastic | rh_viscous;
    const static int rh_ep = rh_elastic | rh_plastic;
    const static int rh_evp = rh_elastic | rh_viscous | rh_plastic;

private:

    // alias to field variables in var
    // ie. var.mat.temperature == var.temperature
    const array_t &coord;
    const conn_t &connectivity;
    const double_vec &temperature;
    const tensor_t &stress;
    const tensor_t &strain_rate;
    const int_vec2D &elemmarkers;

    VectorBase rho0, alpha;
    VectorBase bulk_modulus, shear_modulus;
    VectorBase visc_exponent, visc_coefficient, visc_activation_energy;
    VectorBase heat_capacity, therm_cond;
    VectorBase pls0, pls1;
    VectorBase cohesion0, cohesion1;
    VectorBase friction_angle0, friction_angle1;
    VectorBase dilation_angle0, dilation_angle1;

    #pragma acc routine seq
    void plastic_weakening(int e, double pls,
                           double &cohesion, double &friction_angle,
                           double &dilation_angle, double &hardening) const;
};

#else

double ref_pressure(const Param& param, double z);

class VectorBase
{
public:
    // factory method
    static VectorBase* create(const double_vec&, int);
    virtual ~VectorBase() {};
    virtual double operator[] (std::size_t) const = 0;
    virtual std::size_t size() const = 0;
};


class MatProps
{
public:
    MatProps(const Param& param, const Variables& var);
    ~MatProps();

    const int rheol_type;
    const int nmat;

    double bulkm(int e) const;
    double shearm(int e) const;
    double visc(int e) const;

    double rho(int e) const;
    double cp(int e) const;
    double k(int e) const;

    void plastic_props(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max) const;

    const bool is_plane_strain;
    const double visc_min;
    const double visc_max;
    const double tension_max;
    const double therm_diff_max;

    const static int rh_elastic = 1 << 0;
    const static int rh_viscous = 1 << 1;
    const static int rh_plastic = 1 << 2;
    const static int rh_plastic2d = rh_plastic | 1 << 3;
    const static int rh_maxwell = rh_elastic | rh_viscous;
    const static int rh_ep = rh_elastic | rh_plastic;
    const static int rh_evp = rh_elastic | rh_viscous | rh_plastic;

private:

    // alias to field variables in var
    // ie. var.mat.temperature == var.temperature
    const array_t &coord;
    const conn_t &connectivity;
    const double_vec &temperature;
    const tensor_t &stress;
    const tensor_t &strain_rate;
    const int_vec2D &elemmarkers;

    const VectorBase *rho0, *alpha;
    const VectorBase *bulk_modulus, *shear_modulus;
    const VectorBase *visc_exponent, *visc_coefficient, *visc_activation_energy;
    const VectorBase *heat_capacity, *therm_cond;
    const VectorBase *pls0, *pls1;
    const VectorBase *cohesion0, *cohesion1;
    const VectorBase *friction_angle0, *friction_angle1;
    const VectorBase *dilation_angle0, *dilation_angle1;

    void plastic_weakening(int e, double pls,
                           double &cohesion, double &friction_angle,
                           double &dilation_angle, double &hardening) const;
};

#endif

#endif
