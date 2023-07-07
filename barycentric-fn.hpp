#ifndef DYNEARTHSOL3D_BARYCENTRIC_FN_HPP
#define DYNEARTHSOL3D_BARYCENTRIC_FN_HPP

#include "array2d.hpp"
#include "constants.hpp"
#include "parameters.hpp"

class Barycentric_transformation {

    /* Performing barycentric transformation to a point.
     *
     * The derivation of the formula can be found in
     * http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
     */
    typedef Array2D<double,NODES_PER_ELEM*NDIMS> coeff_t;
    coeff_t coeff_;

public:

    Barycentric_transformation(const array_t &coord,
                               const conn_t &connectivity,
                               const double_vec &volume);
    Barycentric_transformation(const int_vec &elem,
                               const array_t &coord,
                               const conn_t &connectivity,
                               const double_vec &volume);
    Barycentric_transformation(const double** coord,
                               const double volume);
    ~Barycentric_transformation();

    void transform(const double *point, int e, double *result) const;
    bool is_inside_elem(const double *point, int elem) const;
    bool is_inside(const double *result) const;

private:

    inline int index(int node, int dim) const;

#ifdef THREED
    void compute_coeff3d(const double *a,
                         const double *b,
                         const double *c,
                         const double *d,
                         double volume,
                         double *coeff_e);
#else
    void compute_coeff2d(const double *a,
                         const double *b,
                         const double *c,
                         double area,
                         double *coeff_e);
#endif

};

#endif
