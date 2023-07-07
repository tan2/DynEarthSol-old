
#include "barycentric-fn.hpp"


Barycentric_transformation::Barycentric_transformation(const array_t &coord,
                                                       const conn_t &connectivity,
                                                       const double_vec &volume)
    : coeff_(connectivity.size())
{
    #pragma omp parallel for default(none) \
        shared(coord, connectivity, volume)
    for (std::size_t e=0; e<connectivity.size(); ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        const double *d = coord[n3];

        compute_coeff3d(a, b, c, d, volume[e], coeff_[e]);
#else
        compute_coeff2d(a, b, c, volume[e], coeff_[e]);
#endif
    }
}

Barycentric_transformation::Barycentric_transformation(const int_vec &elem,
                                                       const array_t &coord,
                                                       const conn_t &connectivity,
                                                       const double_vec &volume)
    : coeff_(elem.size())
{
    #pragma omp parallel for default(none) \
        shared(elem, coord, connectivity, volume)
    for (std::size_t i=0; i<elem.size(); ++i) {
        int e = elem[i];
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        const double *d = coord[n3];

        compute_coeff3d(a, b, c, d, volume[i], coeff_[i]);
#else
        compute_coeff2d(a, b, c, volume[i], coeff_[i]);
#endif
    }
}

Barycentric_transformation::Barycentric_transformation(const double** coord,
                                                       const double volume)
    : coeff_(1)
{
    const double *a = coord[0];
    const double *b = coord[1];
    const double *c = coord[2];

#ifdef THREED
    const double *d = coord[3];

    compute_coeff3d(a, b, c, d, volume, coeff_[0]);
#else
    compute_coeff2d(a, b, c, volume, coeff_[0]);
#endif
}


Barycentric_transformation::~Barycentric_transformation() {};


void Barycentric_transformation::transform(const double *point, int e, double *result) const
{
    const double *cf = coeff_[e];
    for (int d=0; d<NDIMS; d++) {
        result[d] = cf[index(0,d)];
        for (int i=0; i<NDIMS; i++)
            result[d] += cf[index(i+1,d)]*point[i];
    }
}


bool Barycentric_transformation::is_inside_elem(const double *point, int elem) const
{
    double r[NDIMS];
    transform(point, elem, r);
    return is_inside(r);
}


bool Barycentric_transformation::is_inside(const double *r) const
{
#ifdef THREED

    // 3D has larger round-off error in coeff_
    // => needs greater tolerance
    const double tolerance = 5e-11;

    if (r[0] >= -tolerance &&
        r[1] >= -tolerance &&
        r[2] >= -tolerance &&
        (r[0] + r[1] + r[2]) <= 1 + tolerance)
        return 1;
#else

    const double tolerance = 1e-12;

    if (r[0] >= -tolerance &&
        r[1] >= -tolerance &&
        (r[0] + r[1]) <= 1 + tolerance)
        return 1;
#endif
    return 0;
}


inline int Barycentric_transformation::index(int node, int dim) const
{
    return node*NDIMS + dim;
}


#ifdef THREED

void Barycentric_transformation::compute_coeff3d(const double *a,
                                                 const double *b,
                                                 const double *c,
                                                 const double *d,
                                                 double volume,
                                                 double *coeff_e)
{
    double det = 6 * volume;

    coeff_e[index(0,0)] = (b[0] * (c[1]*d[2] - d[1]*c[2]) +
                           c[0] * (d[1]*b[2] - b[1]*d[2]) +
                           d[0] * (b[1]*c[2] - c[1]*b[2])) / det;
    coeff_e[index(0,1)] = (a[0] * (d[1]*c[2] - c[1]*d[2]) +
                           c[0] * (a[1]*d[2] - d[1]*a[2]) +
                           d[0] * (c[1]*a[2] - a[1]*c[2])) / det;
    coeff_e[index(0,2)] = (a[0] * (b[1]*d[2] - d[1]*b[2]) +
                           b[0] * (d[1]*a[2] - a[1]*d[2]) +
                           d[0] * (a[1]*b[2] - b[1]*a[2])) / det;

    coeff_e[index(1,0)] = ((d[1] - b[1]) * (c[2] - b[2]) -
                           (c[1] - b[1]) * (d[2] - b[2])) / det;
    coeff_e[index(1,1)] = ((c[1] - a[1]) * (d[2] - c[2]) -
                           (d[1] - c[1]) * (c[2] - a[2])) / det;
    coeff_e[index(1,2)] = ((b[1] - d[1]) * (a[2] - d[2]) -
                           (a[1] - d[1]) * (b[2] - d[2])) / det;

    coeff_e[index(2,0)] = ((d[2] - b[2]) * (c[0] - b[0]) -
                           (c[2] - b[2]) * (d[0] - b[0])) / det;
    coeff_e[index(2,1)] = ((c[2] - a[2]) * (d[0] - c[0]) -
                           (d[2] - c[2]) * (c[0] - a[0])) / det;
    coeff_e[index(2,2)] = ((b[2] - d[2]) * (a[0] - d[0]) -
                           (a[2] - d[2]) * (b[0] - d[0])) / det;

    coeff_e[index(3,0)] = ((d[0] - b[0]) * (c[1] - b[1]) -
                           (c[0] - b[0]) * (d[1] - b[1])) / det;
    coeff_e[index(3,1)] = ((c[0] - a[0]) * (d[1] - c[1]) -
                           (d[0] - c[0]) * (c[1] - a[1])) / det;
    coeff_e[index(3,2)] = ((b[0] - d[0]) * (a[1] - d[1]) -
                           (a[0] - d[0]) * (b[1] - d[1])) / det;
}

#else

void Barycentric_transformation::compute_coeff2d(const double *a,
                                                 const double *b,
                                                 const double *c,
                                                 double area,
                                                 double *coeff_e)
{
    double det = 2 * area;

    coeff_e[index(0,0)] = (b[0]*c[1] - b[1]*c[0]) / det;
    coeff_e[index(0,1)] = (c[0]*a[1] - c[1]*a[0]) / det;
    coeff_e[index(1,0)] = (b[1] - c[1]) / det;
    coeff_e[index(1,1)] = (c[1] - a[1]) / det;
    coeff_e[index(2,0)] = (c[0] - b[0]) / det;
    coeff_e[index(2,1)] = (a[0] - c[0]) / det;
}

#endif


