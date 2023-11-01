#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "bc.hpp"


namespace {

void normal_vector_of_facet(int f, const int *conn, const array_t &coord,
                            double *normal, double &zcenter)
{
    int n0 = conn[NODE_OF_FACET[f][0]];
    int n1 = conn[NODE_OF_FACET[f][1]];
#ifdef THREED
    int n2 = conn[NODE_OF_FACET[f][2]];

    // vectors: n0-n1 and n0-n2
    double v01[NDIMS], v02[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = coord[n1][i] - coord[n0][i];
        v02[i] = coord[n2][i] - coord[n0][i];
    }

    // the outward normal vector is parallel to the cross product
    // of v01 and v02, with magnitude = area of the triangle
    normal[0] = (v01[1] * v02[2] - v01[2] * v02[1]) / 2;
    normal[1] = (v01[2] * v02[0] - v01[0] * v02[2]) / 2;
    normal[2] = (v01[0] * v02[1] - v01[1] * v02[0]) / 2;

    zcenter = (coord[n0][2] + coord[n1][2] + coord[n2][2]) / NODES_PER_FACET;
#else
    // the edge vector
    double v01[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = coord[n1][i] - coord[n0][i];
    }

    // the normal vector to the edge, pointing outward
    normal[0] = v01[1];
    normal[1] = -v01[0];

    zcenter = (coord[n0][1] + coord[n1][1]) / NODES_PER_FACET;
#endif
}

}


bool is_on_boundary(const Variables &var, int node)
{
    uint flag = (*var.bcflag)[node];
    return flag & BOUND_ANY;
}


double find_max_vbc(const BC &bc)
{
    double max_vbc_val = 0;
    if (bc.vbc_x0 % 2 == 1) // odd number indicates fixed velocity component
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x0));
    if (bc.vbc_x1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x1));
    if (bc.vbc_y0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y0));
    if (bc.vbc_y1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y1));
    if (bc.vbc_z0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z0));
    if (bc.vbc_z1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z1));
    if (bc.vbc_n0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n0));
    if (bc.vbc_n1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n1));
    if (bc.vbc_n2 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n2));
    if (bc.vbc_n3 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n3));

    return max_vbc_val;
}


void create_boundary_normals(const Variables &var, double bnormals[nbdrytypes][NDIMS],
                             std::map<std::pair<int,int>, double*>  &edge_vectors)
{
    /* This subroutine finds the outward normal unit vectors of boundaries.
     * There are two types of boundaries: ibound{x,y,z}? and iboundn?.
     * The first type ibound{x,y,z}? has well defined "normal" direction
     * as their name implied. Thus they can be curved. (Ex. the top boundary
     * can have topography.)
     * The second type iboundn?, on the other hand, must be a plane so that
     * their vbc are well defined.
     *
     * If the normal component of the boundary velocity is fixed, the boundary
     * normal will not change with time.
     */

    for (int i=0; i<nbdrytypes; i++) {
        double normal[NDIMS] = {0};
        if (var.bfacets[i].size() == 0) continue;

        for (auto j=var.bfacets[i].begin(); j<var.bfacets[i].end(); ++j) {
            int e = j->first;
            int f = j->second;
            double tmp;
            normal_vector_of_facet(f, (*var.connectivity)[e], *var.coord,
                                   normal, tmp);
            // make an unit vector
            double len = 0;
            for(int d=0; d<NDIMS; d++)
                len += normal[d] * normal[d];

            len = std::sqrt(len);
            for(int d=0; d<NDIMS; d++)
                normal[d] = normal[d] / len;

            if (j == var.bfacets[i].begin()) {
                for(int d=0; d<NDIMS; d++)
                    bnormals[i][d] = normal[d];

                if (i < iboundn0) break; // slant boundaries start at iboundn0, other boundary can be curved
            }
            else {
                // Make sure the boundary is a plane, ie. all facets have the same normal vector.
                const double eps2 = 1e-12;
                double diff2 = 0;
                for(int d=0; d<NDIMS; d++)
                    diff2 += (bnormals[i][d] - normal[d]) * (bnormals[i][d] - normal[d]);
                if (diff2 > eps2) {
                    std::cerr << "Error: boundary " << i << " is curved.\n";
                    std::cerr << "Got  -  Expected\n";
                    std::cerr << bnormals[i][0] << " - " << normal[0] << '\n';
                    std::cerr << bnormals[i][1] << " - " << normal[1] << '\n';
#ifdef THREED
                    std::cerr << bnormals[i][2] << " - " << normal[2] << '\n';
#endif
                    std::exit(1);
                }
            }
        }
    }

    for (int i=0; i<nbdrytypes; i++) {
        if (var.bfacets[i].size() == 0) continue;

        const double eps = 1e-15;
        for (int j=i+1; j<nbdrytypes; j++) {
            if (var.bfacets[j].size() == 0) continue;
            double *s = new double[NDIMS];  // intersection of two boundaries
                                            // whole-application lifetime, no need to delete manually
#ifdef THREED
            // quick path: both walls are vertical
            if (std::abs(var.bnormals[i][NDIMS-1]) < eps &&
                std::abs(var.bnormals[j][NDIMS-1]) < eps) {
                s[0] = s[1] = 0;
                s[NDIMS-1] = 1;
            }
            else {
                // cross product of 2 normal vectors
                s[0] = var.bnormals[i][1]*var.bnormals[j][2] - var.bnormals[i][2]*var.bnormals[j][1];
                s[1] = var.bnormals[i][2]*var.bnormals[j][0] - var.bnormals[i][0]*var.bnormals[j][2];
                s[2] = var.bnormals[i][0]*var.bnormals[j][1] - var.bnormals[i][1]*var.bnormals[j][0];
            }
#else
            // 2D
            s[0] = 0;
            s[1] = 1;
#endif
            edge_vectors[std::make_pair(i, j)] = s;
        }
    }
}


void apply_vbcs(const Param &param, const Variables &var, array_t &vel)
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // meaning of vbc flags (odd: free; even: fixed) --
    // 0: all components free
    // 1: normal component fixed, shear components free
    // 2: normal component free, shear components fixed at 0
    // 3: normal component fixed, shear components fixed at 0
    // 4: normal component free, shear component (not z) fixed, only in 3D
    // 5: normal component fixed at 0, shear component (not z) fixed, only in 3D

    const BC &bc = param.bc;

    // diverging x-boundary
    #pragma omp parallel for default(none) \
        shared(bc, var, vel)
    #pragma acc parallel loop
    for (int i=0; i<var.nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(var, i)) continue;

        uint flag = (*var.bcflag)[i];
        double *v = vel[i];

        //
        // X
        //
        if (flag & BOUNDX0) {
            switch (bc.vbc_x0) {
            case 0:
                break;
            case 1:
                v[0] = bc.vbc_val_x0;
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
                v[0] = bc.vbc_val_x0;
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc.vbc_val_x0;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc.vbc_val_x0;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc.vbc_val_x0;
                v[1] = 0;
                break;
#endif
            }
        }
        if (flag & BOUNDX1) {
            switch (bc.vbc_x1) {
            case 0:
                break;
            case 1:
                v[0] = bc.vbc_val_x1;
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
                v[0] = bc.vbc_val_x1;
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc.vbc_val_x1;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc.vbc_val_x1;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc.vbc_val_x1;
                v[1] = 0;
                break;
#endif
            }
        }
#ifdef THREED
        //
        // Y
        //
        if (flag & BOUNDY0) {
            switch (bc.vbc_y0) {
            case 0:
                break;
            case 1:
                v[1] = bc.vbc_val_y0;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = 0;
                v[1] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc.vbc_val_y0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc.vbc_val_y0;
                break;
            }
        }
        if (flag & BOUNDY1) {
            switch (bc.vbc_y1) {
            case 0:
                break;
            case 1:
                v[1] = bc.vbc_val_y1;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = bc.vbc_val_y1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y1;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc.vbc_val_y1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc.vbc_val_y1;
                break;
            }
        }
#endif

        //
        // N
        //
        for (int ib=iboundn0; ib<=iboundn3; ib++) {
            const double eps = 1e-15;
            const double *n = var.bnormals[ib]; // unit normal vector

            if (flag & (1 << ib)) {
                double fac = 0;
                switch (var.vbc_types[ib]) {
                case 1:
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS; d++)
                            vn += v[d] * n[d];  // normal velocity

                        for (int d=0; d<NDIMS; d++)
                            v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        vn += v[d] * n[d];  // normal velocity

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    auto edge = var.edge_vectors.at(std::make_pair(ic, ib));
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 3:
                    for (int d=0; d<NDIMS; d++)
                        v[d] = var.vbc_values[ib] * n[d];  // v must be normal to n
                    break;
                case 11:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS-1; d++)
                            vn += v[d] * n[d];  // normal velocity

                        for (int d=0; d<NDIMS-1; d++)
                            v[d] += (var.vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS-1; d++)
                                        vn += v[d] * n[d];  // normal velocity

                                    for (int d=0; d<NDIMS-1; d++)
                                        v[d] += (var.vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    auto edge = var.edge_vectors.at(std::make_pair(ic, ib));
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 13:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    for (int d=0; d<NDIMS-1; d++)
                        v[d] = var.vbc_values[ib] * fac * n[d];
                    v[NDIMS-1] = 0;
                    break;
                }
            }
        }

        //
        // Z, must be dealt last
        //

        // fast path: vz is usually free in the models
        if (bc.vbc_z0==0 && bc.vbc_z1==0) continue;

        if (flag & BOUNDZ0) {
            switch (bc.vbc_z0) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc.vbc_val_z0;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc.vbc_val_z0;
                break;
            }
        }
        if (flag & BOUNDZ1) {
            switch (bc.vbc_z1) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc.vbc_val_z1;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc.vbc_val_z1;
                break;
            }
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void apply_stress_bcs(const Param& param, const Variables& var, array_t& force)
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // TODO: add general stress (Neumann) bcs

    if (param.control.gravity == 0) return;

    //
    // Gravity-induced (hydrostatic and lithostatic) stress BCs
    //
    for (int i=0; i<nbdrytypes; i++) {
        if (var.vbc_types[i] != 0 &&
            var.vbc_types[i] != 2 &&
            var.vbc_types[i] != 4) continue;

        if (i==iboundz0 && !param.bc.has_winkler_foundation) continue;
        if (i==iboundz1 && !param.bc.has_water_loading) continue;

        const auto& bdry = var.bfacets[i];
        const auto& coord = *var.coord;
        const int nbound = static_cast<int>(bdry.size());
        array_t tmp_val(nbound);
        segment_t tmp_ind(nbound);

        // loops over all bdry facets
        for (int n=0; n<nbound; ++n) {
            // this facet belongs to element e
            int e = bdry[n].first;
            // this facet is the f-th facet of e
            int f = bdry[n].second;
            const int *conn = (*var.connectivity)[e];

            // the outward-normal vector
            double normal[NDIMS];
            // the z-coordinate of the facet center
            double zcenter;

            normal_vector_of_facet(f, conn, *var.coord, normal, zcenter);

            double p;
            if (i==iboundz0 && param.bc.has_winkler_foundation) {
                // Winkler foundation for the bottom boundary
                p = var.compensation_pressure -
                    (var.mat->rho(e) + param.bc.winkler_delta_rho) *
                    param.control.gravity * (zcenter + param.mesh.zlength);
            }
            else if (i==iboundz1 && param.bc.has_water_loading) {
                // hydrostatic water loading for the surface boundary
                const double sea_level = 0;
                p = 0;
                if (zcenter < sea_level) {
                    // below sea level
                    const double sea_water_density = 1030;
                    p = sea_water_density * param.control.gravity * (sea_level - zcenter);
                }
            }
            else {
                // sidewalls
                p = ref_pressure(param, zcenter);
            }

            // lithostatc support - Archimed force (normal to the surface)
            for (int j=0; j<NODES_PER_FACET; ++j)
                tmp_ind[n][j] = conn[NODE_OF_FACET[f][j]];
            for (int d=0; d<NDIMS; ++d)
                tmp_val[n][d] = p * normal[d] / NODES_PER_FACET;
        }

        for (int jj=0; jj<var.bnodes[i].size(); jj++) {
            int node = var.bnodes[i][jj];
            for (int facet=0; facet<nbound; ++facet)
                for (int j=0; j<NODES_PER_FACET; ++j)
                    if (tmp_ind[facet][j] == node) {
                        for (int d=0; d<NDIMS; ++d)
                            force[node][d] -= tmp_val[facet][d];
                        break;
                    }
        }
    }

    if (param.bc.has_elastic_foundation) {
        /* A restoration force on the bottom nodes proportional to total vertical displacement */
        for (auto i=var.bnodes[iboundz0].begin(); i<var.bnodes[iboundz0].end(); ++i) {
            int n = *i;
            force[n][NDIMS-1] -= param.bc.elastic_foundation_constant * ((*var.coord)[n][NDIMS-1] - (*var.coord0)[n][NDIMS-1]);
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

    void simple_diffusion(const Variables& var, array_t& coord,
                          double surface_diffusivity)
    {
#ifdef USE_NPROF
        nvtxRangePush(__FUNCTION__);
#endif
        /* Diffusing surface topography to simulate the effect of erosion and
         * sedimentation.
         */

        const int top_bdry = iboundz1;
        const auto& top = var.bfacets[top_bdry];

        const int_vec& top_nodes = var.bnodes[top_bdry];
        const std::size_t ntop = top_nodes.size();
        double_vec total_dx(var.nnode, 0);
        double_vec total_slope(var.nnode, 0);

        // loops over all top facets
        for (std::size_t i=0; i<top.size(); ++i) {
            // this facet belongs to element e
            int e = top[i].first;
            // this facet is the f-th facet of e
            int f = top[i].second;

            const int *conn = (*var.connectivity)[e];
            int n0 = (*var.connectivity)[e][NODE_OF_FACET[f][0]];
            int n1 = (*var.connectivity)[e][NODE_OF_FACET[f][1]];

#ifdef THREED
            int n2 = (*var.connectivity)[e][NODE_OF_FACET[f][2]];

            double projected_area;
            {
                // normal vector of this facet
                double normal[NDIMS];

                // two vectors n0-n1 and n0-n2
                // n is the cross product of these two vectors
                // the length of n is 2 * triangle area
                double x01, y01, z01, x02, y02, z02;
                x01 = coord[n1][0] - coord[n0][0];
                y01 = coord[n1][1] - coord[n0][1];
                z01 = coord[n1][2] - coord[n0][2];
                x02 = coord[n2][0] - coord[n0][0];
                y02 = coord[n2][1] - coord[n0][1];
                z02 = coord[n2][2] - coord[n0][2];

                normal[0] = y01*z02 - z01*y02;
                normal[1] = z01*x02 - x01*z02;
                normal[2] = x01*y02 - y01*x02;

                /* the area of this facet is:
                 *   0.5 * std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
                 *
                 * theta is the angle between this facet and horizontal
                 *   tan_theta = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]) / normal[2]
                 *   cos_theta = normal[2] / (2 * area)
                 *
                 * the area projected on the horizontal plane is:
                 *   projected_area = area * cos_theta = 0.5 * normal[2]
                 */
                projected_area = 0.5 * normal[2];
            }

            total_dx[n0] += projected_area;
            total_dx[n1] += projected_area;
            total_dx[n2] += projected_area;

            double shp2dx[NODES_PER_FACET], shp2dy[NODES_PER_FACET];
            double iv = 1 / (2 * projected_area);
            shp2dx[0] = iv * (coord[n1][1] - coord[n2][1]);
            shp2dx[1] = iv * (coord[n2][1] - coord[n0][1]);
            shp2dx[2] = iv * (coord[n0][1] - coord[n1][1]);
            shp2dy[0] = iv * (coord[n2][0] - coord[n1][0]);
            shp2dy[1] = iv * (coord[n0][0] - coord[n2][0]);
            shp2dy[2] = iv * (coord[n1][0] - coord[n0][0]);

            double D[NODES_PER_FACET][NODES_PER_FACET];
            for (int j=0; j<NODES_PER_FACET; j++) {
                for (int k=0; k<NODES_PER_FACET; k++) {
                    D[j][k] = (shp2dx[j] * shp2dx[k] +
                               shp2dy[j] * shp2dy[k]);
                }
            }

            const int n[NODES_PER_FACET] = {n0, n1, n2};
            for (int j=0; j<NODES_PER_FACET; j++) {
                double slope = 0;
                for (int k=0; k<NODES_PER_FACET; k++)
                    slope += D[j][k] * coord[n[k]][2];

                total_slope[n[j]] += slope * projected_area;
            }

            // std::cout << i << ' ' << n0 << ' ' << n1 << ' ' << n2 << "  "
            //           << projected_area << "  " << slope << '\n';
#else
            /* The 1D diffusion operation is implemented ad hoc, not using FEM
             * formulation (e.g. computing shape function derivation on the edges).
             */

            double dx = std::fabs(coord[n1][0] - coord[n0][0]);
            total_dx[n0] += dx;
            total_dx[n1] += dx;

            double slope = (coord[n1][1] - coord[n0][1]) / dx;
            total_slope[n0] -= slope;
            total_slope[n1] += slope;

            // std::cout << i << ' ' << n0 << ' ' << n1 << "  " << dx << "  " << slope << '\n';
#endif
        }

        double max_dh = 0;
        for (std::size_t i=0; i<ntop; ++i) {
            // we don't treat edge nodes specially, i.e. reflecting bc is used for erosion.
            int n = top_nodes[i];
            double dh = surface_diffusivity * var.dt * total_slope[n] / total_dx[n];
            coord[n][NDIMS-1] -= dh;
            max_dh = std::max(max_dh, std::fabs(dh));
            // std::cout << n << "  dh:  " << dh << '\n';
        }

        // std::cout << "max erosion / sedimentation rate (cm/yr):  "
        //           << max_dh / var.dt * 100 * YEAR2SEC << '\n';
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void custom_surface_processes(const Variables& var, array_t& coord)
    {
        const int top_bdry = iboundz1;
        const int_vec& top_nodes = var.bnodes[top_bdry];
        const std::size_t ntop = top_nodes.size();

        // loops over all top nodes
        for (std::size_t i=0; i<ntop; ++i) {
            int n = top_nodes[i];

            // coordinate of this node
            double x = coord[n][0];
            double z = coord[n][NDIMS-1];

            // compute topography changes here ...
            // positive dh means sedimentation, negative dh means erosion
            double dh = 0;
            {
                // dh == ....
            }

            coord[n][NDIMS-1] += dh;
        }
    }

}


void surface_processes(const Param& param, const Variables& var, array_t& coord)
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    switch (param.control.surface_process_option) {
    case 0:
        // no surface process
        break;
    case 1:
        simple_diffusion(var, coord, param.control.surface_diffusivity);
        break;
    case 101:
        custom_surface_processes(var, coord);
        break;
    default:
        std::cout << "Error: unknown surface process option: " << param.control.surface_process_option << '\n';
        std::exit(1);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


