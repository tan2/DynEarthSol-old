#include <algorithm>
#include <cstring>
#include <functional>
#include <iostream>
#include <numeric>
#include <unordered_map>

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
#include "brc-interpolation.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "nn-interpolation.hpp"
#include "utils.hpp"
#include "markerset.hpp"
#include "remeshing.hpp"

#ifdef ADAPT

/* Copyright (C) 2009 Imperial College London.

 Please see the AUTHORS file in the main source directory for a full
 list of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 USA
*/
#include <cmath>
#include <vector>

#include "vtk.h"
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>

#include "ErrorMeasure.h"
#include "Adaptivity.h"
#include "DiscreteGeometryConstraints.h"
#include <assert.h>

#endif // end of ifdef ADAPT

#ifdef USEMMG
#ifdef THREED
#include "mmg/src/mmg3d/libmmg3d.h"
#else
#include "mmg/src/mmg2d/libmmg2d.h"
#endif
#define MAX0(a,b)     (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX0(a,b)) > (MAX0(c,d))) ? (MAX0(a,b)) : (MAX0(c,d)))
#endif

namespace {


// equilateral triangle area = 0.433*s^2, equilateral tetrahedron volume = 0.118*s^3
#ifdef THREED
const double sizefactor = 0.118;
#else
const double sizefactor = 0.433;
#endif

const int DELETED_FACET = -1;
const int DEBUG = 0;

bool is_boundary(uint flag)
{
    return flag & BOUND_ANY;
}


bool is_bottom(uint flag)
{
    return flag & BOUNDZ0;
}


bool is_corner(uint flag)
{
    uint f = flag & BOUND_ANY;
    if (!f) return 0;

    // A corner node will have multiple bits (2 in 2D; 3 or more in 3D) set in its flag.
    int nbits = 0;
    for (int j=0; j<nbdrytypes; j++) {
        // counting how many bits are set
        if (f & (1<<j)) nbits++;
    }

#ifdef THREED
    return (nbits >= NDIMS);
#else
    return (nbits == NDIMS);
#endif
}


bool is_bottom_corner(uint flag)
{
    if ((flag & BOUNDZ0) && is_corner(flag)) return 1;
    return 0;
}


void flatten_bottom(const uint_vec &old_bcflag, double *qcoord,
                    double bottom, int_vec &points_to_delete, double min_dist)
{
    // find old nodes that are on or close to the bottom boundary

    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_bottom(flag)) {
            // restore edge nodes to initial depth
            qcoord[i*NDIMS + NDIMS-1] = bottom;
        }
        else if (! is_boundary(flag) &&
                 std::fabs(qcoord[i*NDIMS + NDIMS-1] - bottom) < min_dist) {
            points_to_delete.push_back(i);
        }
    }
}


void new_bottom(const uint_vec &old_bcflag, double *qcoord,
                double bottom_depth, int_vec &points_to_delete, double min_dist,
                int *segment, int *segflag, int nseg)
{
    /* deleting nodes that are on or close to the bottom boundary,
     * excluding nodes on the side walls
     */

#ifdef THREED
    /* In 3D, if bottom nodes were deleted, the facets on the side boundaries
     * near the bottom are affected as well. The code does not deal with this
     * complexity and can only work in 2D.
     */
    std::cerr << "Error: new_bottom() does not work in 3D.\n";
    std::exit(1);
#endif

    int_vec bottom_corners;
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_bottom(flag)) {
            if(is_bottom_corner(flag))
                bottom_corners.push_back(i);
            else
                points_to_delete.push_back(i);
        }
        else if (! is_boundary(flag) &&
                 std::fabs(qcoord[i*NDIMS + NDIMS-1] - bottom_depth) < min_dist) {
            points_to_delete.push_back(i);
        }
    }

    if (DEBUG) {
        std::cout << "bottom points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
        std::cout << "segment before delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag before delete: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }

    // must have 2 bottom corners in 2D
    if (bottom_corners.size() != 2) {
        std::cerr << "Error: cannot find all bottom corners before remeshing. n_bottom_corners = "
                  << bottom_corners.size() << " (2 expected).\n";
        std::cout << "bottom corners: ";
        print(std::cout, bottom_corners);
        std::cout << '\n';
        std::exit(11);
    }

    // move the corners to the same depth
    for (std::size_t i=0; i<bottom_corners.size(); i++) {
        int n = bottom_corners[i];
        qcoord[n*NDIMS + NDIMS-1] = bottom_depth;
    }

    // mark all bottom nodes as deleted
    for (int i=0; i<nseg; ++i) {
        if (static_cast<uint>(segflag[i]) == BOUNDZ0) {
            for (int j=0; j<NODES_PER_FACET; j++)
                segment[i*NODES_PER_FACET + j] = DELETED_FACET;
        }
    }

    // create new bottom segments from corner nodes
    for (int i=0; i<nseg; ++i) {
        if (static_cast<uint>(segflag[i]) == BOUNDZ0) {
            segment[i*NODES_PER_FACET + 0] = bottom_corners[0];
            segment[i*NODES_PER_FACET + 1] = bottom_corners[1];
            break;
        }
    }

    if (DEBUG) {
        std::cout << "bottom corners: ";
        print(std::cout, bottom_corners);
        std::cout << '\n';
        std::cout << "segment with new bottom: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag with new bottom: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }
}


// local cmp functor
struct cmp {
    const array_t &coord;
    const int d;
    cmp (const array_t &coord_, int dim) : coord(coord_), d(dim) {};
    bool operator()(const int &a, const int &b) {return coord[a][d] < coord[b][d];}
};


typedef std::pair<int,int> edge_t;
struct equal_to1
{
    bool operator()(const edge_t &lhs, const edge_t &rhs) const {
        return (lhs.first == rhs.first && lhs.second == rhs.second);
    }
};


struct hash1
{
    std::size_t operator()(const edge_t &k) const {
        // first as the upper half
        const int halfbytes = 4;
        return (static_cast<std::size_t>(k.first) << sizeof(std::size_t)*halfbytes) | k.second;
    }
};


void assemble_bdry_polygons(const Variables &var, const array_t &old_coord,
                             const conn_t &old_connectivity,
                             int_vec (&bdry_polygons)[nbdrytypes])
{
    /* bdry_polygons[i] contains a polygon, ie. a list of vertex, enclosing the i-th boundary */

#ifdef THREED
    const int nodes_per_edge = 2;  // an edge has 2 nodes
    const int edges_per_facet = 3;  // a facet has 3 edges
    const int edgenodes[edges_per_facet][nodes_per_edge] = { {1, 2}, {2, 0}, {0, 1} };

    for (int ibound=0; ibound<nbdrytypes; ibound++) {
        if (var.bnodes[ibound].size() == 0) continue;  // skip empty boundary

        //
        // Collecting edges on each boundary.
        //

        std::unordered_map<edge_t, int, hash1, equal_to1> edges;
        for (std::size_t i=0; i<var.bfacets[ibound].size(); i++) {
            const auto &facet = var.bfacets[ibound][i];
            int e = facet.first;
            int f = facet.second;
            const int *conn = old_connectivity[e];

            for (int j=0; j<edges_per_facet; j++) {
                int n0 = conn[ NODE_OF_FACET[f][ edgenodes[j][0] ] ];
                int n1 = conn[ NODE_OF_FACET[f][ edgenodes[j][1] ] ];
                if (n0 > n1) {  // ensure n0 < n1
                    int tmp;
                    tmp = n0;
                    n0 = n1;
                    n1 = tmp;
                }
                edge_t g = std::make_pair(n0, n1);
                auto search = edges.find(g);
                if (search == edges.end()) {
                    edges[g] = 1;
                }
                else {
                    search->second ++;
                }
            }
        }
        if (DEBUG > 1) {
            std::cout << ibound << "-th edge:\n";
            for (auto kk=edges.begin(); kk!=edges.end(); ++kk) {
                std::cout << kk->first.first << ",\t" << kk->first.second << "\t: " << kk->second << '\n';
            }
            std::cout << '\n';
        }

        //
        // Collecting edges enclosing the boundary
        //
        std::vector<const edge_t*> enclosing_edges;
        for (auto kk=edges.begin(); kk!=edges.end(); ++kk) {
            int count = kk->second;
            if (count == 2) {
                // this is an "internal" edge
                continue;
            }
            else if (count == 1) {
                // this edge encloses the boundary
                enclosing_edges.push_back(&(kk->first));
            }
            else {
                // not possible
                std::cout << "Error: an edge is belonged to more than 2 facets. The mesh is corrupted.\n";
                std::exit(11);
            }
        }

        //
        // Connecting edges to form a polygon
        //

        int_vec &polygon = bdry_polygons[ibound];
        const auto g0 = enclosing_edges.begin();
        int head = (*g0)->first;
        int tail = (*g0)->second;
        polygon.push_back(head);
        polygon.push_back(tail);
        auto g1 = g0+1;
        while (head != tail) {
            for (auto g=g1; g!=enclosing_edges.end(); ++g) {
                if ((*g)->first == tail) {
                    tail = (*g)->second;
                    polygon.push_back(tail);
                }
                else if ((*g)->second == tail) {
                    tail = (*g)->first;
                    polygon.push_back(tail);
                }
            }
        }
        // the starting point and end point must be the same
        if (polygon.front() != polygon.back()) {
            std::cout << "Error: boundary polygon is not closed. The mesh is corrupted.\n";
            std::exit(11);
        }
        polygon.pop_back();  // removed the duplicating end point

        if (DEBUG > 1) {
            std::cout << "nodes for " << ibound << "-th boundary polygon:\n";
            print(std::cout, polygon);
            std::cout << '\n';
        }
    }

#endif
}


void find_tiny_element(const Param &param, const double_vec &volume,
                       int_vec &tiny_elems)
{
    const double smallest_vol = param.mesh.smallest_size * sizefactor * std::pow(param.mesh.resolution, NDIMS);

    for (std::size_t e=0; e<volume.size(); e++) {
        if (volume[e] < smallest_vol)
            tiny_elems.push_back(e);
    }

    if (DEBUG) {
        std::cout << "tiny elements: ";
        print(std::cout, tiny_elems);
        std::cout << '\n';
    }
}


void find_points_of_tiny_elem(const array_t &coord, const conn_t &connectivity,
                              const double_vec &volume, const int_vec &tiny_elems,
                              int npoints, const double *old_points,
                              const uint_vec &old_bcflag, int_vec &points_to_delete,
                              bool excl_func(uint))
{
    // collecting the nodes of tiny_elems
    int tiny_nelem = tiny_elems.size();
    array_t tiny_coord(tiny_nelem * NODES_PER_ELEM);
    conn_t tiny_conn(tiny_nelem);
    double_vec tiny_vol(tiny_nelem);
    int ii = 0;
    for (int ee=0; ee<tiny_nelem; ++ee) {
        int e = tiny_elems[ee];

        tiny_vol[ee] = volume[e];

        const int *conn = connectivity[e];
        for (int j=0; j<NODES_PER_ELEM; ++j) {
            int n = conn[j];
            tiny_conn[ee][j] = ii;

            for (int d=0; d<NDIMS; ++d) {
                tiny_coord[ii][d] = coord[n][d];
            }
            ii ++;
        }
    }

    Barycentric_transformation bary(tiny_coord, tiny_conn, tiny_vol);

    // find old nodes that are connected to tiny elements and are not excluded
    // (most of the nodes of tiny elements are newly inserted by the remeshing library)
    for (int i=0; i<npoints; ++i) {
        // excluded nodes
        if (excl_func(old_bcflag[i])) continue;

        const double *p = old_points + i*NDIMS;
        for (int ee=0; ee<tiny_nelem; ++ee) {
            if (bary.is_inside_elem(p, ee)) {
                points_to_delete.push_back(i);
                break;
            }
        }
    }

    if (DEBUG) {
        std::cout << "points of tiny elements: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
    }
}


void delete_points(const int_vec &points_to_delete, int &npoints,
                   int nseg, double *points, int *segment)
{
    if (points_to_delete.size() == 0) return;

    if (DEBUG) {
        std::cout << "old points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
    }

    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = npoints - 1;

    // delete points from the end
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        // when a point is deleted, replace it with the last point
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        std::replace(segment, endsegment, end, *i);
        // std::cout << *i << " <- " << end << "\n";

        end --;
    }
    npoints -= points_to_delete.size();
}


void delete_facets(int &nseg, int *segment, int *segflag)
{
    // delete facets from the end
    for (int i=nseg-1; i>=0; i--) {
        if (segment[i*NODES_PER_FACET] == DELETED_FACET) {
            // safety check
            if (segment[i*NODES_PER_FACET + 1] != DELETED_FACET
#ifdef THREED
                || segment[i*NODES_PER_FACET + 2] != DELETED_FACET
#endif
                ) {
                std::cerr << "Error: segment array is corrupted before delete_facets()!\n";
                print(std::cerr, segment, nseg*NODES_PER_FACET);
                std::exit(11);
            }

            // replace deleted segment with the last segment
            for (int j=0; j<NODES_PER_FACET; ++j) {
                segment[i*NODES_PER_FACET + j] = segment[(nseg-1)*NODES_PER_FACET + j];
            }
            segflag[i] = segflag[nseg-1];
            nseg --;
        }
    }

    if (DEBUG) {
        std::cout << "segment: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }
}


void delete_points_and_merge_segments(const int_vec &points_to_delete, int &npoints,
                                      int nseg, double *points, int *segment,
                                      uint_vec &bcflag, double min_length)
{
#ifdef THREED
    std::cerr << "delete_points_and_merge_segments() doesn't work in 3D!\n";
    std::exit(12);
#endif

    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = npoints - 1;

    // delete points from the end
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        if (DEBUG) {
            std::cout << " deleting point " << *i << " replaced by point " << end << '\n';
        }

        // if the deleted point is a segment point, merge the neighboring two segments
        uint flag = bcflag[*i];
        if (is_boundary(flag)) {
            // segment and points:  aa------a = b-------bb

            int *a = std::find(segment, endsegment, *i);
            int *b = std::find(a+1, endsegment, *i);
            if (b == endsegment) {
                std::cerr << "Error: segment array is corrupted when merging segment!\n";
                std::exit(11);
            }

            // a could be the either first or second node of the segment,
            // which will affect the location of aa in the segment array.
            bool not_first;
            not_first = (a - segment) % NODES_PER_FACET;
            int *aa = not_first? (a - 1) : (a + 1);
            not_first = (b - segment) % NODES_PER_FACET;
            int *bb = not_first? (b - 1) : (b + 1);

            // if the length of the two segments are too long
            // do not delete this point
            double la2, lb2;
            la2 = dist2(points + (*a)*NDIMS, points + (*aa)*NDIMS);
            lb2 = dist2(points + (*b)*NDIMS, points + (*bb)*NDIMS);
            if (la2 > min_length*min_length && lb2 > min_length*min_length) {
                if (DEBUG) {
                    std::cout << " the segments of point " << *i << " have length^2 "
                              << la2 << ", " << lb2 << " -- skip deletion."<< '\n';
                }
                continue;
            }

            // merge the two segments
            *a = *bb;
            *bb = DELETED_FACET;
            *b = DELETED_FACET;

            if (DEBUG) {
                std::cout << "a: " << (a-segment) << "  b: " << (b-segment)
                          << " not_first? " << not_first << " bb: " << (bb-segment)
                          << " = " << *bb << '\n';
                std::cout << "segment after merging: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        // when a point is deleted, replace it with the last point
        flag = bcflag[end];
        bcflag[*i] = flag;
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        if (is_boundary(flag)) {
            std::replace(segment, endsegment, end, *i);
            // std::cout << *i << " <- " << end << "\n";
            if (DEBUG) {
                std::cout << "segment after replace: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        end --;
        npoints --;
    }

    /* PS: We don't check whether the merged segments belong to the same boundary or not,
     *     e.g., one segment is on the top boundary while the other segment is on the left
     *     boundary. The check is performed in delete_facets() instead.
     */
}


void delete_points_and_merge_facets(const int_vec &points_to_delete,
                                    const int_vec (&bnodes)[nbdrytypes],
                                    const int_vec (&bdry_polygons)[nbdrytypes],
                                    const int_vec (&bdrynode_deleting)[nbdrytypes],
                                    const double (&bnormals)[nbdrytypes][NDIMS],
                                    int &npoints,
                                    int &nseg, double *points,
                                    int *segment, int *segflag,
                                    uint_vec &bcflag, double min_length)
{
#ifndef THREED
    std::cerr << "delete_points_and_merge_facets() doesn't work in 2D!\n";
    std::exit(12);
#else

    int_vec inverse[nbdrytypes], nfacets;
    std::vector<int*> facet;

    // before deleting boundary points, create a new triangulation
    // TODO: skip boundary bdrynode_deleting.size()==0, but retaining its facets
    for (int i=0; i<nbdrytypes; ++i) {  // looping over all possible boundaries
        const int_vec& bdeleting = bdrynode_deleting[i];  // nodes to be deleted on this boundary, sorted
        const int_vec& bdry_nodes = bnodes[i];  // all nodes on this boundary, sorted
        if (bdry_nodes.size() == 0) continue;

        if (DEBUG) {
            std::cout << i << "-th boundary to be merged: ";
            print(std::cout, bdeleting);
            std::cout << '\n';
        }

        // 3d coordinate -> 2d
        Array2D<double,2> coord2d(bdry_nodes.size() - bdeleting.size());
        int_vec& inv = inverse[i];
        {
            if (DEBUG > 1) {
                std::cout << "bdry nodes:\n";
                print(std::cout, bdry_nodes);
                std::cout << "\n";
                std::cout << bdry_nodes.size() << ' ' << bdeleting.size() << '\n';
            }

            int major_direction = i / 2;  // largest component of normal vector
            if (i >= iboundn0) { // slant boundaries start at iboundn0
                double max = 0;
                for (int d=0; d<NDIMS; d++) {
                    if (std::abs(bnormals[i][d]) > max) {
                        max = std::abs(bnormals[i][d]);
                        major_direction = d;
                    }
                }
            }

            for (std::size_t j=0, k=0, n=0; j<bdry_nodes.size(); j++) {

                // is j belonged to bdeleting[]?
                if (k < bdeleting.size() && bdry_nodes[j] == bdeleting[k]) {
                    k++;
                    continue;
                }

                // record the node # for inverse mapping
                inv.push_back(bdry_nodes[j]);

                int dd = 0;
                for (int d=0; d<NDIMS; d++) {
                    if (d == major_direction) continue;
                    coord2d[n][dd] = points[bdry_nodes[j]*NDIMS + d];
                    dd++;
                }
                n++;
            }

            if (DEBUG > 1) {
                std::cout << i << "-th boundary to be remeshed: ";
                print(std::cout, coord2d);
                std::cout << '\n';
            }
        }

        // re-triangulate the affected boundary (connect only, not adding new points)
        {
            // converting polygon vertex to segment array
            const int_vec& polygon = bdry_polygons[i];
            int_vec surf_segflag(polygon.size()); // all 0, its value does not matter
            int_vec surf_segment(2 * polygon.size());
            std::size_t first = 0;
            int new_polygon_size = 0;
            for (std::size_t j=0; j<polygon.size(); j++) {
                auto search = std::find(inv.begin(), inv.end(), polygon[j]);
                if (search != inv.end()) {  // the vertex is not deleted
                    std::size_t ia = search - inv.begin();
                    if (new_polygon_size == 0) {
                        // start of the polygon
                        surf_segment[0] = ia;
                        first = ia;
                    }
                    else {
                        surf_segment[2*new_polygon_size-1] = surf_segment[2*new_polygon_size] = ia;
                    }
                    ++new_polygon_size;
                }
            }
            // end of the polygon
            surf_segment[2*new_polygon_size-1] = first;

            if (DEBUG) {
                std::cout << "inverse: \n";
                print(std::cout, inv);
                std::cout << '\n';
                std::cout << "polygon segments: ";
                print(std::cout, surf_segment, new_polygon_size*2);
                std::cout << '\n';
            }

            // temporary arrays, some will be allocated inside Triangle library
            int new_nnode, new_nelem, new_nseg;
            double *pcoord, *pregattr;
            int *pconnectivity, *psegment, *psegflag;

            Mesh mesh;
            mesh.min_angle = 0;
            mesh.meshing_verbosity = 0;
            points_to_new_surface(mesh, coord2d.size(), coord2d.data(),
                                  new_polygon_size, surf_segment.data(), surf_segflag.data(),
                                  0, NULL,
                                  0, 3,
                                  new_nnode, new_nelem, new_nseg,
                                  pcoord, pconnectivity, psegment, psegflag, pregattr);

            if (static_cast<std::size_t>(new_nnode) != coord2d.size()) {
                std::cerr << "Error: ponits_to_new_surface is adding new points!\n";
                std::cout << new_nnode << ' ' << coord2d.size() << '\n';
                std::cout << "old points: ";
                print(std::cout, coord2d);
                std::cout << '\n';
                std::cout << "new points: ";
                print(std::cout, pcoord, new_nnode*2);
                std::cout << '\n';
                std::cout << "new conn: ";
                print(std::cout, pconnectivity, new_nelem*3);
                std::cout << '\n';
                std::exit(12);
            }
            if (new_nseg != new_polygon_size) {
                std::cerr << "Error: points_to_new_surface is adding new segments!\n";
                std::exit(12);
            }

            delete [] pcoord;
            delete [] psegment;
            delete [] psegflag;
            delete [] pregattr;
            // remember to free pconnectivity later

            // translating index of local (per-boundary) node # to global (whole-mesh) node #
            for (int j=0; j<new_nelem; ++j) {
                for (int k=0; k<NODES_PER_FACET; ++k) {
                    int n = pconnectivity[NODES_PER_FACET*j + k];
                    pconnectivity[NODES_PER_FACET*j + k] = inv[n];
                }
            }

            // storing the new boundary facets for later
            // ownership of "pconnectivity" array is transferred to "facet"
            facet.push_back(pconnectivity);
            nfacets.push_back(new_nelem);
        }
    }

    int nseg2 = std::accumulate(nfacets.begin(), nfacets.end(), 0);
    if (nseg2 > nseg) {
        std::cerr << "Error: ponits_to_new_surface too many segments!\n";
        std::exit(12);
    }

    // appending facets of all boundaries into segment array
    for (int i=0, n=0; i<nbdrytypes; ++i) {
        if (bnodes[i].size() == 0) continue;

        for (int k=0; k<nfacets[i]; ++k, ++n) {
            for (int j=0; j<NODES_PER_FACET; ++j)
                segment[n*NODES_PER_FACET + j] = facet[i][k*NODES_PER_FACET + j];
            segflag[n] = 1 << i;
        }
        delete [] facet[i];
    }

    // mark deleted facets
    for (int i=nseg2; i<nseg; ++i) {
        for (int j=0; j<NODES_PER_FACET; ++j)
            segment[i*NODES_PER_FACET + j] = DELETED_FACET;
        segflag[i] = 0;
    }
    nseg = nseg2;

    // delete points from the end
    int *endsegment = segment + nseg * NODES_PER_FACET;  // last segment
    int end = npoints - 1;  // last point
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        if (DEBUG) {
            std::cout << " deleting point " << *i << " replaced by point " << end << '\n';
        }

        // when a point is deleted, replace it with the last point
        uint flag = bcflag[end];
        bcflag[*i] = flag;
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        if (is_boundary(flag)) {
            std::replace(segment, endsegment, end, *i);
            if (DEBUG > 1) {
                std::cout << *i << " <- " << end << "\n";
                std::cout << "segment after replace: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        end --;
        npoints --;
    }

#endif
}


void delete_points_on_boundary(int_vec &points_to_delete,
                               const int_vec (&bnodes)[nbdrytypes],
                               const int_vec (&bdry_polygons)[nbdrytypes],
                               const double (&bnormals)[nbdrytypes][NDIMS],
                               int &npoints,
                               int &nseg, double *points,
                               int *segment, int *segflag,
                               uint_vec &bcflag, double min_size)
{
    if (DEBUG > 1) {
        std::cout << "old points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
        std::cout << "segment before delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
    }

#ifdef THREED
    // are there any points in points_to_delete on the boundary? and on which boundary?
    // if the deleted point is a boundary point, store its index
    bool changed = 0;
    int_vec bdrynode_deleting[nbdrytypes];
    for (auto i=points_to_delete.begin(); i<points_to_delete.end(); ++i) {
        uint flag = bcflag[*i];
        for (int j=0; j<nbdrytypes; ++j) {
            uint bc = 1 << j;
            if (flag & bc) {
                bdrynode_deleting[j].push_back(*i);
                changed = 1;
            }
        }
    }

    // non-boundary points changed,
    if (! changed) {
        delete_points(points_to_delete, npoints, nseg,
                      points, segment);
        delete_facets(nseg, segment, segflag);
        return;
    }

    delete_points_and_merge_facets(points_to_delete, bnodes, bdry_polygons,
                                   bdrynode_deleting, bnormals, npoints, nseg,
                                   points, segment, segflag, bcflag, min_size);
    delete_facets(nseg, segment, segflag);
#else
    delete_points_and_merge_segments(points_to_delete, npoints, nseg,
                                     points, segment, bcflag, min_size);
    delete_facets(nseg, segment, segflag);
    // silence 'unused variable' complier warning
    (void) bnodes;
    (void) bdry_polygons;
#endif

    if (DEBUG > 1) {
        std::cout << "segment after  delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
    }
}


void new_mesh(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
    int_vec bdry_polygons[nbdrytypes];
    assemble_bdry_polygons(var, original_coord, original_connectivity, bdry_polygons);

    // create a copy of original mesh
    array_t old_coord(original_coord);
    conn_t old_connectivity(original_connectivity);
    segment_t old_segment(original_segment);
    segflag_t old_segflag(original_segflag);

    // raw pointers to old mesh
    double *qcoord = old_coord.data();
    int *qconn = old_connectivity.data();
    int *qsegment = old_segment.data();
    int *qsegflag = old_segflag.data();

    // size of old mesh
    int old_nnode = old_coord.size();
    int old_nelem = old_connectivity.size();
    int old_nseg = old_segment.size();

    // copying useful arrays of old mesh
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = var.bnodes[i];  // copying whole vector
    }

    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted
    switch (param.mesh.remeshing_option) {
    case 0:
    case 1:
    case 2:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 10:
    case 11:
        // DO NOT change the corners
        excl_func = &is_corner;
        break;
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        std::exit(1);
    }

    /* choosing which way to remesh the boundary */
    int_vec points_to_delete;
    const double min_dist = std::pow(param.mesh.smallest_size*sizefactor, 1./NDIMS) * param.mesh.resolution;
    switch (param.mesh.remeshing_option) {
    case 0:
    case 10:
        // Nothing
        break;
    case 1:
    case 11:
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 2:
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
        break;
    }

    if (bad_quality == 3) { // there is a tiny element
        // Marking (non-boundary) points of small elements, which will be deleted later
        int_vec tiny_elems;
        find_tiny_element(param, old_volume, tiny_elems);

        if (tiny_elems.size() > 0) {
            find_points_of_tiny_elem(old_coord, old_connectivity, old_volume,
                                     tiny_elems, old_nnode, qcoord, old_bcflag, points_to_delete,
                                     excl_func);
        }
    }

    // sort points_to_delete and remove duplicates
    {
        std::sort(points_to_delete.begin(), points_to_delete.end());
        auto last = std::unique(points_to_delete.begin(), points_to_delete.end());
        points_to_delete.resize(last - points_to_delete.begin());
    }

    // delete points
    switch (param.mesh.remeshing_option) {
    case 0:
    case 1:
    case 2:
        // deleting non-boundary points
        delete_points(points_to_delete, old_nnode, old_nseg,
                      qcoord, qsegment);
        delete_facets(old_nseg, qsegment, qsegflag);
        break;
    case 10:
    case 11:
        // deleting points, some of them might be on the boundary
        delete_points_on_boundary(points_to_delete, old_bnodes, bdry_polygons, var.bnormals,
                                  old_nnode, old_nseg,
                                  qcoord, qsegment, qsegflag, old_bcflag, min_dist);
        break;
    }

    int new_nnode, new_nelem, new_nseg;
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    int nloops = 0;
    Mesh mesh = param.mesh;  // temporary copy
    mesh.poly_filename = ""; // not used, but still need to be init'd
    while (1) {

        if (bad_quality == 3) {
            // lessen the quality constraint so that less new points got inserted to the mesh
            // and less chance of having tiny elements
            mesh.min_angle *= 0.9;
            mesh.max_ratio *= 0.9;
            mesh.min_tet_angle *= 1.1;
        }
#ifdef THREED
        if (nloops != 0 && bad_quality == 1) {
            // enable tetgen optimization for mesh with higher quality
            // tetgen might consume too much memory and crash!
            mesh.tetgen_optlevel = 3;
        }
#endif
        pregattr = NULL;

        //
        // new mesh
        //

        // We don't want to refine large elements during remeshing,
        // so using negative size as the max area
        const double max_elem_size = -1;
        const int vertex_per_polygon = 3;
        points_to_new_mesh(mesh, old_nnode, qcoord,
                           old_nseg, qsegment, qsegflag,
                           0, pregattr,
                           max_elem_size, vertex_per_polygon,
                           new_nnode, new_nelem, new_nseg,
                           pcoord, pconnectivity, psegment, psegflag, pregattr);

        array_t new_coord(pcoord, new_nnode);
        conn_t new_connectivity(pconnectivity, new_nelem);

        // deleting (non-boundary) nodes to avoid having tiny elements
        double_vec new_volume(new_nelem);
        compute_volume(new_coord, new_connectivity, new_volume);

        const double smallest_vol = param.mesh.smallest_size * sizefactor * std::pow(param.mesh.resolution, NDIMS);
        bad_quality = 0;
        for (int e=0; e<new_nelem; e++) {
            if (new_volume[e] < smallest_vol) {
                bad_quality = 3;
                break;
            }
        }
        int worst_elem;
        double q = worst_elem_quality(new_coord, new_connectivity,
                                      new_volume, worst_elem);
#ifdef THREED
        // normalizing q so that its magnitude is about the same in 2D and 3D
        q = std::pow(q, 1.0/3);
#endif
        if (q < param.mesh.min_quality) {
            bad_quality = 1;
        }

        new_coord.nullify();
        new_connectivity.nullify();
        if (! bad_quality) break;

        nloops ++;
        if (nloops > 5) {
            std::cout << "Warning: exceeding loop limit in remeshing. Proceeding with risks.\n";
            break;
        }

        delete [] pcoord;
        delete [] pconnectivity;
        delete [] psegment;
        delete [] psegflag;
        delete [] pregattr;
    }

    var.nnode = new_nnode;
    var.nelem = new_nelem;
    var.nseg = new_nseg;
    var.coord->reset(pcoord, new_nnode);
    var.connectivity->reset(pconnectivity, new_nelem);
    var.segment->reset(psegment, var.nseg);
    var.segflag->reset(psegflag, var.nseg);
}

void compute_metric_field(const Variables &var, const conn_t &connectivity, const double resolution0, double_vec &metric, double_vec &tmp_result_sg)
{
    /* dvoldt is the volumetric strain rate, weighted by the element volume,
     * lumped onto the nodes.
     */
    const double_vec& volume = *var.volume;
    const double_vec& volume_n = *var.volume_n;
    double resolution = resolution0;
    std::fill_n(metric.begin(), var.nnode, 0);

    #pragma omp parallel for default(none) shared(var,volume,connectivity,tmp_result_sg,resolution)
    for (int e=0;e<var.nelem;e++) {
        const int *conn = connectivity[e];
        double plstrain = resolution/(1.0+5.0*(*var.plstrain)[e]);
        tmp_result_sg[e] = plstrain * volume[e];
    }

    #pragma omp parallel for default(none) shared(var,metric,tmp_result_sg,volume_n)
    for (int n=0;n<var.nnode;n++) {
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e)
            metric[n] += tmp_result_sg[*e];
        metric[n] /= volume_n[n];
    }
}

#ifdef USEMMG
#ifdef THREED
void optimize_mesh(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
    // We don't want to refine large elements during remeshing,
    // so using negative size as the max area
    const double max_elem_size = -1;
    const int vertex_per_polygon = 3;
    const double min_dist = std::pow(param.mesh.smallest_size*sizefactor, 1./NDIMS) * param.mesh.resolution;
    Mesh mesh_param = param.mesh;
    mesh_param.poly_filename = "";

    int_vec bdry_polygons[nbdrytypes];    
    assemble_bdry_polygons(var, original_coord, original_connectivity,
                           bdry_polygons);

    // create a copy of original_coord and original_segment
    array_t old_coord(original_coord);
    conn_t old_connectivity(original_connectivity);
    conn_t old_connectivity_from_1(original_connectivity);
    segment_t old_segment(original_segment);
    segment_t old_segment_from_1(original_segment);
    segflag_t old_segflag(original_segflag);

    double *qcoord = old_coord.data();
    int *qconn = old_connectivity.data();
    int *qsegment = old_segment.data();
    int *qconn_from_1 = old_connectivity_from_1.data();
    int *qsegment_from_1 = old_segment_from_1.data();
    int *qsegflag = old_segflag.data();

    int old_nnode = old_coord.size();
    int old_nelem = old_connectivity.size();
    int old_nseg = old_segment.size();

    // copy
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = var.bnodes[i];
    }

    int_vec points_to_delete;
    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted

    /* choosing which way to remesh the boundary */
    switch (param.mesh.remeshing_option) {
    case 0:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 1:
        excl_func = &is_boundary;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 2:
        excl_func = &is_boundary;
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
        break;
    case 10:
        excl_func = &is_corner;
        break;
    case 11:
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        std::exit(1);
    }

    // --- STEP I: Initialization
    // 1) Initialisation of mesh and sol structures
    //   args of InitMesh:
    //     MMG5_ARG_start: we start to give the args of a variadic func
    //     MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
    //     &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
    //     MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
    //     &mmgSol: pointer toward your MMG5_pSol (that store your metric) */
    MMG5_pMesh      mmgMesh = NULL;
    MMG5_pSol       mmgSol  = NULL;
  
    MMG3D_Init_mesh(MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet,
                    &mmgSol, MMG5_ARG_end);

    // 2) Build mesh in MMG5 format
    // Manually set of the mesh 
    //  a) give the size of the mesh: vertices, tetra, prisms, triangles, quads, edges
    if ( MMG3D_Set_meshSize(mmgMesh, old_nnode, old_nelem,
                            0,old_nseg,0,0) != 1 )
        exit(EXIT_FAILURE);
    //   b) give the vertex coordinates. References are NULL but can be an integer array for boundary flag etc.
    if( MMG3D_Set_vertices(mmgMesh, qcoord, NULL) != 1)
        exit(EXIT_FAILURE);
    //   c) give the connectivity. References are NULL but can be an integer array for boundary flag etc.
    for (int i = 0; i < old_nelem*NODES_PER_ELEM; ++i)
        ++qconn_from_1[i];
    if( MMG3D_Set_tetrahedra(mmgMesh, qconn_from_1, NULL) != 1 )
        exit(EXIT_FAILURE);
    //   d) give the segments (i.e., boundary facet)
    for (int i = 0; i < old_nseg*NODES_PER_FACET; ++i)
        ++qsegment_from_1[i];
    if( MMG3D_Set_triangles(mmgMesh, qsegment_from_1, qsegflag) != 1 )
        exit(EXIT_FAILURE);

    // 3) Build sol in MMG5 format
    //      Here a 'solution' is a nodal field that becomes
    //      the basis for metric tensor for isotropic and anisotropic
    //      mesh adaptation.
    //
    //   a) give info for the sol structure
    if( MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, old_nnode, MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
    //   b) give solutions values and positions
    compute_metric_field(var, old_connectivity, param.mesh.resolution, *var.ntmp, *var.tmp_result_sg);
    //      i) If sol array is available:
    if( MMG3D_Set_scalarSols(mmgSol, (*var.ntmp).data()) != 1 )
        exit(EXIT_FAILURE);

    /*save init mesh*/
    //      ii) Otherwise, set a value node by node:
    // for (int i = 0; i < var.nnode; ++i) {
    //     if( MMG3D_Set_scalarSol(mmgSol, 100.0, i+1) != 1 )
    //         exit(EXIT_FAILURE);
    // }
    //      iii) If a metric field ('solution') is given, 
    //           optimization mode should be off.
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_optim, 0) != 1 ) 
        exit(EXIT_FAILURE);

    // 4) (not mandatory): check if the number of given entities match with mesh size
    if( MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1 ) exit(EXIT_FAILURE);

    //--- STEP  II: Remesh function
    /* debug mode ON (default value = OFF) */
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_debug, param.mesh.mmg_debug) != 1 )
        exit(EXIT_FAILURE);
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, param.mesh.mmg_verbose) != 1 )
        exit(EXIT_FAILURE);
 
    /* maximal memory size (default value = 50/100*ram) */
    //if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_mem, 600) != 1 )
    //exit(EXIT_FAILURE);

    // /* Maximal mesh size (default FLT_MAX)*/
    if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmax, param.mesh.mmg_hmax_factor*param.mesh.resolution) != 1 )
    exit(EXIT_FAILURE);

    /* Minimal mesh size (default 0)*/
    if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmin, param.mesh.mmg_hmin_factor*param.mesh.resolution) != 1 )
    exit(EXIT_FAILURE);

    /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
    if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hausd, param.mesh.mmg_hausd_factor*param.mesh.resolution) != 1 )
    exit(EXIT_FAILURE);

    // /* Gradation control*/
    // if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 3.0) != 1 )
    // exit(EXIT_FAILURE);

    // /* Gradation requirement */
    // if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgradreq, -1.0) != 1 )
    // exit(EXIT_FAILURE);

    const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);
    if ( ier == MMG5_STRONGFAILURE ) {
        fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
        exit(EXIT_FAILURE);
    } else if ( ier == MMG5_LOWFAILURE ) {
        fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");
        exit(EXIT_FAILURE);
    }

    //--- STEP III: Get results
    // 1) Preparations
    //   a) get the size of the mesh: vertices, tetra, triangles, edges */
    int na;
    if ( MMG3D_Get_meshSize(mmgMesh, &(var.nnode), &(var.nelem), NULL, &(var.nseg), NULL, &na) !=1 )
        exit(EXIT_FAILURE);
    std::cerr << "Updated mesh size\n";
    std::cerr << "New number of vertices:" << var.nnode << std::endl;
    std::cerr << "New number of elements:" << var.nelem << std::endl;
    std::cerr << "New number of segments:" << var.nseg << std::endl;

    //   b) Create mesh-defining arrays of new sizes
    array_t new_coord( var.nnode );
    conn_t new_connectivity( var.nelem );
    segment_t new_segment( var.nseg );
    segflag_t new_segflag( var.nseg );
    std::cerr << "Resized arrays\n";

    // double *ncoord = new_coord.data();
    // int *nconn = new_connectivity.data();
    // int *nsegment = new_segment.data();
    // int *nsegflag = new_segflag.data();

    // 2) Pupolate DES3D mesh-defining arrays
    //   a) Vertex recovering
    for (int i = 0; i < var.nnode; ++i) {
        if ( MMG3D_Get_vertex(mmgMesh, &(new_coord[i][0]), &(new_coord[i][1]), &(new_coord[i][2]), NULL, NULL, NULL) != 1 )
            exit(EXIT_FAILURE);
    }
    std::cerr << "New coordinates populated\n";

    //   b) Tetra recovering
    for (int i = 0; i < var.nelem; ++i) {
        if ( MMG3D_Get_tetrahedron(mmgMesh, &(new_connectivity[i][0]), &(new_connectivity[i][1]), &(new_connectivity[i][2]), &(new_connectivity[i][3]), NULL, NULL) != 1 )  
            exit(EXIT_FAILURE);
        for(std::size_t j = 0; j < NODES_PER_ELEM; ++j)
            new_connectivity[i][j] -= 1;
    }
    std::cerr << "New connectivity populated\n";

    //   c) segments recovering
    for (int i = 0; i < var.nseg; ++i) {
        if ( MMG3D_Get_triangle(mmgMesh, &(new_segment[i][0]), &(new_segment[i][1]), &(new_segment[i][2]),&(new_segflag.data()[i]), NULL) != 1 )
            exit(EXIT_FAILURE);
        for(int j = 0; j < NODES_PER_FACET; ++j)
            new_segment[i][j] -= 1;
    }     
    std::cerr << "New segments populated\n";

    //   d) Let the DES3D arrays point to the newly populated data 
    var.coord->steal_ref( new_coord );
    var.connectivity->steal_ref( new_connectivity );
    var.segment->steal_ref( new_segment );
    var.segflag->steal_ref( new_segflag );
    std::cerr << "Arrays transferred." << std::endl;

    // 3) Free the MMG3D5 structures
    MMG3D_Free_all(MMG5_ARG_start,
                    MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                    MMG5_ARG_end);
    std::cerr << "MMG3D freed." <<std::endl;

    std::cerr << "\nMesh optimization done" <<std::endl;
}

#else

void optimize_mesh_2d(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
    // We don't want to refine large elements during remeshing,
    // so using negative size as the max area
    const double max_elem_size = -1;
    const int vertex_per_polygon = 3;
    const double min_dist = std::pow(param.mesh.smallest_size*sizefactor, 1./NDIMS) * param.mesh.resolution;
    Mesh mesh_param = param.mesh;
    mesh_param.poly_filename = "";

    int_vec bdry_polygons[nbdrytypes];    
    assemble_bdry_polygons(var, original_coord, original_connectivity,
                           bdry_polygons);

    // create a copy of original_coord and original_segment
    array_t old_coord(original_coord);
    conn_t old_connectivity(original_connectivity);
    conn_t old_connectivity_from_1(original_connectivity);
    segment_t old_segment(original_segment);
    segment_t old_segment_from_1(original_segment);
    segflag_t old_segflag(original_segflag);

    double *qcoord = old_coord.data();
    int *qconn = old_connectivity.data();
    int *qsegment = old_segment.data();
    int *qconn_from_1 = old_connectivity_from_1.data();
    int *qsegment_from_1 = old_segment_from_1.data();
    int *qsegflag = old_segflag.data();

    int old_nnode = old_coord.size();
    int old_nelem = old_connectivity.size();
    int old_nseg = old_segment.size();

    // copy
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = var.bnodes[i];
    }

    int_vec points_to_delete;
    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted

    /* choosing which way to remesh the boundary */
    switch (param.mesh.remeshing_option) {
    case 0:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 1:
        excl_func = &is_boundary;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 2:
        excl_func = &is_boundary;
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
        break;
    case 10:
        excl_func = &is_corner;
        break;
    case 11:
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        std::exit(1);
    }

    // --- STEP I: Initialization
    // 1) Initialisation of mesh and sol structures
    //   args of InitMesh:
    //     MMG5_ARG_start: we start to give the args of a variadic func
    //     MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
    //     &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
    //     MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
    //     &mmgSol: pointer toward your MMG5_pSol (that store your metric) */
    MMG5_pMesh      mmgMesh = NULL;
    MMG5_pSol       mmgSol  = NULL;
  
    MMG2D_Init_mesh(MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet,
                    &mmgSol, MMG5_ARG_end);

    // 2) Build mesh in MMG5 format
    // Manually set of the mesh 
    //  a) give the size of the mesh: vertices, triangles, quads(=0), edges
    if ( MMG2D_Set_meshSize(mmgMesh, old_nnode, old_nelem, 0, old_nseg) != 1 )
        exit(EXIT_FAILURE);
    //   b) give the vertex coordinates. References are NULL but can be an integer array for boundary flag etc.
    if( MMG2D_Set_vertices(mmgMesh, qcoord, NULL) != 1)
        exit(EXIT_FAILURE);
    //   c) give the connectivity. References are NULL but can be an integer array for boundary flag etc.
    for (int i = 0; i < old_nelem*NODES_PER_ELEM; ++i)
        ++qconn_from_1[i];
    if( MMG2D_Set_triangles(mmgMesh, qconn_from_1, NULL) != 1 )
        exit(EXIT_FAILURE);
    //   d) give the segments (i.e., boundary edges)
    for (int i = 0; i < old_nseg*NODES_PER_FACET; ++i)
        ++qsegment_from_1[i];
    for (int i = 0; i < old_nseg; ++i)        
        if( MMG2D_Set_edge(mmgMesh, qsegment_from_1[i*NODES_PER_FACET], 
                qsegment_from_1[i*NODES_PER_FACET+1], qsegflag[i], i+1) != 1)
            exit(EXIT_FAILURE);
    // if( MMG2D_Set_edges(mmgMesh, qsegment_from_1, qsegflag) != 1 )
    //     exit(EXIT_FAILURE);

    // 3) Build sol in MMG5 format
    //      Here a 'solution' is a nodal field that becomes
    //      the basis for metric tensor for isotropic and anisotropic
    //      mesh adaptation.
    //
    //   a) give info for the sol structure
    if( MMG2D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, old_nnode, MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
    //   b) give solutions values and positions
    compute_metric_field(var, old_connectivity, param.mesh.resolution, *var.ntmp, *var.tmp_result_sg);
    //      i) If sol array is available:
    if( MMG2D_Set_scalarSols(mmgSol, (*var.ntmp).data()) != 1 )
        exit(EXIT_FAILURE);
    //      ii) Otherwise, set a value node by node:
    // for (int i = 0; i < var.nnode; ++i) {
    //     if( MMG2D_Set_scalarSol(mmgSol, 0.5, i+1) != 1 )
    //         exit(EXIT_FAILURE);
    // }
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_optim, 0) != 1 )
    exit(EXIT_FAILURE);

    // 4) (not mandatory): check if the number of given entities match with mesh size
    if( MMG2D_Chk_meshData(mmgMesh, mmgSol) != 1 ) exit(EXIT_FAILURE);

    //--- STEP  II: Remesh function
    /* debug mode ON (default value = OFF) */
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_debug, param.mesh.mmg_debug) != 1 )
    exit(EXIT_FAILURE);

    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose, param.mesh.mmg_verbose) != 1 )
    exit(EXIT_FAILURE);

    // if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_iso, 1) != 1 )
    // exit(EXIT_FAILURE);
 
    /* maximal memory size (default value = 50/100*ram) */
    //if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_mem, 600) != 1 )
    //exit(EXIT_FAILURE);

    /* Maximal mesh size (default FLT_MAX)*/
    if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax, param.mesh.mmg_hmax_factor*param.mesh.resolution) != 1 )
    exit(EXIT_FAILURE);

    /* Minimal mesh size (default 0)*/
    if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin, param.mesh.mmg_hmin_factor*param.mesh.resolution) != 1 )
    exit(EXIT_FAILURE);

    /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
    if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd, param.mesh.mmg_hausd_factor*param.mesh.resolution) != 1 )
    exit(EXIT_FAILURE);

    // /* Gradation control*/
    // if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hgrad, 3.0) != 1 )
    // exit(EXIT_FAILURE);

    // /* Gradation requirement */
    // if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hgradreq, 3.0) != 1 )
    // exit(EXIT_FAILURE);

    const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);
    if ( ier == MMG5_STRONGFAILURE ) {
        fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
        exit(EXIT_FAILURE);
    } else if ( ier == MMG5_LOWFAILURE ) {
        fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");
        exit(EXIT_FAILURE);
    }    

    //--- STEP III: Get results
    // 1) Preparations
    //   a) get the size of the mesh: vertices, tetra, triangles, edges */
    if ( MMG2D_Get_meshSize(mmgMesh, &(var.nnode), &(var.nelem), NULL, &(var.nseg)) !=1 )
        exit(EXIT_FAILURE);
    std::cerr << "Updated mesh size\n";
    std::cerr << "New number of vertices:" << var.nnode << std::endl;
    std::cerr << "New number of elements:" << var.nelem << std::endl;
    std::cerr << "New number of segments:" << var.nseg << std::endl;

    //   b) Create mesh-defining arrays of new sizes
    array_t new_coord( var.nnode );
    conn_t new_connectivity( var.nelem );
    segment_t new_segment( var.nseg );
    segflag_t new_segflag( var.nseg );
    std::cerr << "Resized arrays\n";

    // 2) Pupolate DES3D mesh-defining arrays
    //   a) Vertexes recovering
    for (int i = 0; i < var.nnode; ++i) {
        if ( MMG2D_Get_vertex(mmgMesh, &(new_coord[i][0]), &(new_coord[i][1]), NULL, NULL, NULL) != 1 )
            exit(EXIT_FAILURE);
    }
    std::cerr << "New coordinates populated\n";

    //   b) Triangles recovering
    for (int i = 0; i < var.nelem; ++i) {
        if ( MMG2D_Get_triangle(mmgMesh, &(new_connectivity[i][0]), &(new_connectivity[i][1]), &(new_connectivity[i][2]), NULL, NULL) != 1 )  
            exit(EXIT_FAILURE);
        for(std::size_t j = 0; j < NODES_PER_ELEM; ++j)
            new_connectivity[i][j] -= 1;
    }
    std::cerr << "New connectivity populated\n";

    //   c) segments recovering
    for (int i = 0; i < var.nseg; ++i) {
        if ( MMG2D_Get_edge(mmgMesh, &(new_segment[i][0]), &(new_segment[i][1]), &(new_segflag.data()[i]), NULL, NULL) != 1 )
            exit(EXIT_FAILURE);
        for(std::size_t j = 0; j < NODES_PER_FACET; ++j)
            new_segment[i][j] -= 1;
    }
    std::cerr << "New segments populated\n";

    //   d) Let the DES3D arrays point to the newly populated data 
    var.coord->steal_ref( new_coord );
    var.connectivity->steal_ref( new_connectivity );
    var.segment->steal_ref( new_segment );
    var.segflag->steal_ref( new_segflag );
    std::cerr << "Arrays transferred." << std::endl;

    // 3) Free the MMG3D5 structures
    MMG2D_Free_all(MMG5_ARG_start,
                    MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                    MMG5_ARG_end);
    std::cerr << "MMG2D freed." <<std::endl;

    std::cerr << "\nMesh optimization done" <<std::endl;
}
#endif  // end of if THREED
#endif // end of if USEMMG

#if defined THREED && defined ADAPT
void optimize_mesh(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
    // We don't want to refine large elements during remeshing,
    // so using negative size as the max area
    const double max_elem_size = -1;
    const int vertex_per_polygon = 3;
    const double min_dist = std::pow(param.mesh.smallest_size*sizefactor, 1./NDIMS) * param.mesh.resolution;
    Mesh mesh_param = param.mesh;
    mesh_param.poly_filename = "";

    int_vec bdry_polygons[nbdrytypes];    
    assemble_bdry_polygons(var, original_coord, original_connectivity,
                           bdry_polygons);

    // create a copy of original_coord and original_segment
    array_t old_coord(original_coord);
    conn_t old_connectivity(original_connectivity);
    segment_t old_segment(original_segment);
    segflag_t old_segflag(original_segflag);

    double *qcoord = old_coord.data();
    int *qconn = old_connectivity.data();
    int *qsegment = old_segment.data();
    int *qsegflag = old_segflag.data();

    int old_nnode = old_coord.size();
    int old_nelem = old_connectivity.size();
    int old_nseg = old_segment.size();

    // copy
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = var.bnodes[i];
    }

    int_vec points_to_delete;
    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted

    /* choosing which way to remesh the boundary */
    switch (param.mesh.remeshing_option) {
    case 0:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 1:
        excl_func = &is_boundary;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 2:
        excl_func = &is_boundary;
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
        break;
    case 10:
        excl_func = &is_corner;
        break;
    case 11:
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        std::exit(1);
    }

    ////// Optimize mesh using libadaptivity

    // 1. Create a vtkUnstructuredGrid (UG) object.
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

    //// 1.1 add old points to the UG object.
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (std::size_t i = 0; i < old_nnode; ++i) {
        points->InsertNextPoint( &qcoord[i*NDIMS] );
    }
    ug->SetPoints(points);

    //// 1.2 add old connectivity to the UG object.
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    for (std::size_t i = 0; i < old_nelem; ++i) {
        vtkSmartPointer<vtkTetra> tetra =  vtkSmartPointer<vtkTetra>::New();
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            tetra->GetPointIds()->SetId(j, qconn[i*NODES_PER_ELEM + j]);
        cells->InsertNextCell( tetra );
    }
    ug->SetCells(VTK_TETRA, cells);

#if VTK_MAJOR_VERSION <= 5
    ug->Update();
#endif

    // These should be populated only once at the outset of a simulation
    // and be maintained thereafter.
    std::vector<int> SENList, sids;
    {
        DiscreteGeometryConstraints constraints;
        constraints.verbose_off();
        constraints.set_coplanar_tolerance(0.9999999);
        constraints.set_volume_input(ug);

        sids.assign(qsegflag, qsegflag + old_nseg);
        SENList.assign(qsegment, qsegment + old_nseg * NODES_PER_FACET);
        //constraints.get_coplanar_ids(sids);
        //constraints.get_surface(SENList);
        //assert(sids.size()*3==SENList.size());
        //cout<<"Found "<<sids.size()<<" surface elements\n";
    }

    DiscreteGeometryConstraints constraints;
    constraints.verbose_off();
    constraints.set_surface_input(ug, SENList, sids);

    std::vector<double> max_len;
    constraints.get_constraints(max_len);

    // Prepare the field to be used for error analysis: e.g., plastic strain or strain rate.
    // Setup data for the triangle. Attach a value of 1.45.
    // This can be anything you wish to store with it)
    vtkSmartPointer<vtkDoubleArray> cellData = vtkSmartPointer<vtkDoubleArray>::New();
    cellData->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    cellData->SetName("plasticStrain"); //set the name of the value
    for (std::size_t i = 0; i < old_nelem; ++i)
        cellData->InsertNextValue( (*var.plstrain)[i] ); //set the actual value
    ug->GetCellData()->AddArray(cellData);

    vtkSmartPointer<vtkCellDataToPointData> cell2point = vtkSmartPointer<vtkCellDataToPointData>::New();
#if VTK_MAJOR_VERSION <= 5
    cell2point->SetInput(ug);
#else
    cell2point->SetInputData(ug);
#endif
    cell2point->PassCellDataOff();
    cell2point->Update();

    ug->GetPointData()->AddArray(cell2point->GetUnstructuredGridOutput()->GetPointData()->GetArray("plasticStrain"));

    /* Test merging of metrics.
     */
    ErrorMeasure error;
    error.verbose_off();
    error.set_input(ug);
    //error.add_field("plasticStrain", 0.05, false, 0.01); // 1.0, false, 1.0
    error.add_field_simple("plasticStrain", 2.0*param.mesh.resolution, false, 0.01); // 1.0, false, 1.0
    //error.set_max_length(param.mesh.resolution*5.0);
    //error.set_max_length(&(max_len[0]), ug->GetNumberOfPoints()); // 1 or NumofPoints.
    error.set_min_length(param.mesh.resolution * param.mesh.smallest_size);
    error.apply_gradation(1.3);
    error.set_max_nodes(5*old_nnode);

    // // For debugging
    // constraints.write_vtk(std::string("sids_before.vtu"));

    // error.diagnostics();

    // vtkXMLUnstructuredGridWriter *metric_writer = vtkXMLUnstructuredGridWriter::New();
    // metric_writer->SetFileName("metric_new.vtu");
    // metric_writer->SetInput(ug);
    // metric_writer->Write();
    // metric_writer->Delete();
    // ug->GetPointData()->RemoveArray("mean_desired_lengths");
    // ug->GetPointData()->RemoveArray("desired_lengths");

    // vtkXMLUnstructuredGridWriter *ug_writer = vtkXMLUnstructuredGridWriter::New();
    // ug_writer->SetFileName("before_adapted.vtu");
    // ug_writer->SetInput(ug);
    // ug_writer->Write();
    // /////////////////

    Adaptivity adapt;
    adapt.verbose_on();
    adapt.set_from_vtk(ug, true);
    adapt.set_adapt_sweeps(5);
    adapt.set_surface_mesh(SENList);
    adapt.set_surface_ids(sids);
    adapt.adapt();
    adapt.get_surface_ids(sids);
    adapt.get_surface_mesh(SENList);
    vtkSmartPointer<vtkUnstructuredGrid> adapted_ug = adapt.get_adapted_vtu();
    // ug->Delete();

    // // For debugging
    // ug_writer->SetFileName("after_adapted.vtu");
    // ug_writer->SetInput(adapted_ug);
    // ug_writer->Write();

    // constraints.set_surface_input(adapted_ug, SENList, sids);
    // constraints.write_vtk(std::string("sids_after.vtu"));
    ////////////////


    // update mesh info.
    var.nnode = adapted_ug->GetNumberOfPoints();
    var.nelem = adapted_ug->GetNumberOfCells();
    var.nseg = sids.size();
    std::cerr << "Updated mesh size\n";

    array_t new_coord( var.nnode );
    conn_t new_connectivity( var.nelem );
    segment_t new_segment( var.nseg );
    segflag_t new_segflag( var.nseg );
    std::cerr << "Resized arrays\n";

    for (std::size_t i = 0; i < var.nnode; ++i) {
        double *x = adapted_ug->GetPoints()->GetPoint(i);
        for(int j=0; j < NDIMS; j++ )
            new_coord[i][j] = x[j];
    }
    std::cerr << "New coordinates populated\n";

    for (std::size_t i = 0; i < var.nelem; ++i) {
        vtkSmartPointer<vtkTetra> tetra = (vtkTetra *)adapted_ug->GetCell(i);
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            new_connectivity[i][j] = tetra->GetPointId(j);
    }
    std::cerr << "New connectivity populated\n";

    // copy optimized surface triangle connectivity to dynearthsol3d.
    for (std::size_t i = 0; i < var.nseg; ++i) {
        for(int j=0; j < NODES_PER_FACET; j++ )
            new_segment[i][j] = SENList[i*NODES_PER_FACET + j];
        new_segflag.data()[i] = sids[i];
    }
    std::cerr << "New segments populated\n";

    var.coord->steal_ref( new_coord );
    var.connectivity->steal_ref( new_connectivity );
    var.segment->steal_ref( new_segment );
    var.segflag->steal_ref( new_segflag );
    std::cerr << "Arrays transferred. Mesh optimization done \n";
}
#endif

} // anonymous namespace


int bad_mesh_quality(const Param &param, const Variables &var, int &index)
{
    /* Check the quality of the mesh, return 0 if the mesh quality (by several
     * measures) is good. Non-zero returned values indicate --
     * 1: an element has bad quality (too acute / narrow / flat).
     * 2: a bottom node has moved too far away from the flat bottom.
     * 3: an element is smaller than (mesh.smallest_size * [volume of a equilateral triangle/tetrahedron
     *    of side = mesh.resolution]).
     */

    // check tiny elements
    const double smallest_vol = param.mesh.smallest_size * sizefactor * std::pow(param.mesh.resolution, NDIMS);
    for (int e=0; e<var.nelem; e++) {
        if ((*var.volume)[e] < smallest_vol) {
            index = e;
            std::cout << "    The size of element #" << index << " is too small.\n";
            return 3;
        }
    }

    // check if any bottom node is too far away from the bottom depth
    if (param.mesh.remeshing_option == 1 ||
        param.mesh.remeshing_option == 2 ||
        param.mesh.remeshing_option == 11) {
        double bottom = - param.mesh.zlength;
        const double dist = param.mesh.max_boundary_distortion * param.mesh.resolution;
        for (int i=0; i<var.nnode; ++i) {
            if (is_bottom((*var.bcflag)[i])) {
                double z = (*var.coord)[i][NDIMS-1];
                if (std::fabs(z - bottom) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the bottm: z = " << z << "\n";
                    return 2;
                }
            }
        }
    }

    // check element distortion
    int worst_elem;
    double q = worst_elem_quality(*var.coord, *var.connectivity,
                                  *var.volume, worst_elem);
#ifdef THREED
    // normalizing q so that its magnitude is about the same in 2D and 3D
    q = std::pow(q, 1.0/3);
#endif
    if (q < param.mesh.min_quality) {
        index = worst_elem;
        std::cout << "    Element #" << worst_elem << " has mesh quality = " << q << ".\n";
        return 1;
    }
    return 0;
}


void remesh(const Param &param, Variables &var, int bad_quality)
{
    std::cout << "  Remeshing starts...\n";

    {
        // creating a "copy" of mesh pointer so that they are not deleted
        array_t old_coord;
        conn_t old_connectivity;
        segment_t old_segment;
        segflag_t old_segflag;
        old_coord.steal_ref(*var.coord);
        old_connectivity.steal_ref(*var.connectivity);
        old_segment.steal_ref(*var.segment);
        old_segflag.steal_ref(*var.segflag);

#ifdef THREED
#if defined ADAPT || defined USEMMG
        optimize_mesh(param, var, bad_quality, old_coord, old_connectivity,
                 old_segment, old_segflag);
#else
        new_mesh(param, var, bad_quality, old_coord, old_connectivity,
                 old_segment, old_segflag);
#endif
#else  // if 2d
#if defined USEMMG
        optimize_mesh_2d(param, var, bad_quality, old_coord, old_connectivity,
                 old_segment, old_segflag);
#else
        new_mesh(param, var, bad_quality, old_coord, old_connectivity,
                 old_segment, old_segflag);
#endif
#endif
        renumbering_mesh(param, *var.coord, *var.connectivity, *var.segment, NULL);        

        {
            std::cout << "    Interpolating fields.\n";
            Barycentric_transformation bary(old_coord, old_connectivity, *var.volume);

            // interpolating fields defined on elements
            nearest_neighbor_interpolation(var, bary, old_coord, old_connectivity);

            // interpolating fields defined on nodes
            barycentric_node_interpolation(var, bary, old_coord, old_connectivity);
        }

        delete var.support;
        create_support(var);

        std::cout << "    Remapping markers.\n";
        // remap markers. elemmarkers are updated here, too.
        remap_markers(param, var, old_coord, old_connectivity);
  
        // old_coord et al. are destroyed before exiting this block
    }

    // memory for new fields
    reallocate_variables(param, var);

    // updating other arrays
    create_boundary_flags(var);
    for (int i=0; i<nbdrytypes; ++i) {
        var.bnodes[i].clear();
        var.bfacets[i].clear();
    }
    create_boundary_nodes(var);
    create_boundary_facets(var);
    /* // moved before remap_markers()
     * delete var.support;
     * create_support(var);
     */

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    // TODO: using edvoldt and volume to get volume_old
    std::copy(var.volume->begin(), var.volume->end(), var.volume_old->begin());
    compute_mass(param, var, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.tmp_result,
                 *var.support);
    compute_shape_fn(var, *var.coord, *var.connectivity, *var.volume,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    if (param.mesh.remeshing_option==1 ||
        param.mesh.remeshing_option==2 ||
        param.mesh.remeshing_option==11) {
        /* Reset coord0 of the bottom nodes */
        for (auto i=var.bnodes[iboundz0].begin(); i<var.bnodes[iboundz0].end(); ++i) {
            int n = *i;
            (*var.coord0)[n][NDIMS-1] = -param.mesh.zlength;
        }
    }

    if (param.sim.has_output_during_remeshing) {
        // the following variables need to be re-computed only when we are
        // outputing right after remeshing
        update_strain_rate(var, *var.strain_rate);
        update_force(param, var, *var.force, *var.tmp_result);
    }

    std::cout << "  Remeshing finished.\n";
}


