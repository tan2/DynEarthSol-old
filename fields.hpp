#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP

void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot);
void update_strain_rate(const Variables& var, tensor_t& strain_rate);
void update_force(const Param& param, const Variables& var, array_t& force);
void update_velocity(const Variables& var, array_t& vel);
void update_coordinate(const Variables& var, array_t& coord);
void rotate_stress();


#endif