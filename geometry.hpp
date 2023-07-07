#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

double dist2(const double* a, const double* b);
void compute_volume(const double **coord, double &volume);
void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume);
void compute_volume(const Variables &var, double_vec &volume);

void compute_dvoldt(const Variables &var, double_vec &dvoldt,
                    double_vec &tmp_result_sg);

void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt);

void NMD_stress(const Param& param, const Variables &var, double_vec &dp_nd,
                tensor_t& stress, double_vec &tmp_result_sg);

double compute_dt(const Param& param, const Variables& var);

void compute_mass(const Param &param, const Variables &var,
                  double max_vbc_val, double_vec &volume_n,
                  double_vec &mass, double_vec &tmass,
                  elem_cache &tmp_result);

void compute_shape_fn(const Variables &var, shapefn &shpdx,
                      shapefn &shpdy, shapefn &shpdz);

double elem_quality(const array_t &coord, const conn_t &connectivity,
                    const double_vec &volume, int e);

double worst_elem_quality(const array_t &coord, const conn_t &connectivity,
                          const double_vec &volume, int &worst_elem);

#endif
