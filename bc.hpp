#ifndef DYNEARTHSOL3D_BC_HPP
#define DYNEARTHSOL3D_BC_HPP

bool is_on_boundary(const Variables &var, int node);
double find_max_vbc(const BC &bc);
void create_boundary_normals(const Variables &var, array_t &bnormals,
                             std::map<std::pair<int,int>, double*>  &edge_vectors);
void apply_vbcs(const Param &param, const Variables &var, array_t &vel);
void apply_stress_bcs(const Param& param, const Variables& var, array_t& force);
void surface_plstrain_diffusion(const Param &param, const Variables& var, double_vec& plstrain);
void correct_surface_element(const Variables& var, const double_vec& dhacc, MarkerSet& ms, tensor_t& stress, \
                              tensor_t& strain, tensor_t& strain_rate, double_vec& plstrain);
void surface_processes(const Param& param, const Variables& var, array_t& coord, tensor_t& stress, tensor_t& strain, \
                       tensor_t& strain_rate, double_vec& plstrain, SurfaceInfo& surfinfo, \
                       std::vector<MarkerSet*> &markersets, int_vec2D& elemmarkers);

#endif
