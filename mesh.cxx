#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <sstream>

#ifdef THREED

#define TETLIBRARY
#include "tetgen/tetgen.h"
#undef TETLIBRARY

#endif // THREED

#ifdef USEEXODUS
#include "netcdf.h"
#include "exodusII.h"
#endif // USEEXODUS

#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS


#include "constants.hpp"
#include "parameters.hpp"
#include "sortindex.hpp"
#include "utils.hpp"
#include "mesh.hpp"
#include "markerset.hpp"

#ifdef WIN32
#ifndef _MSC_VER
namespace std {
  static std::string to_string(long double t)
  {
    char temp[32];
    sprintf(temp,"%f",double(t));
    return std::string(temp);
  }
}
#endif //_MSC_VER
#endif // WIN32


namespace { // anonymous namespace

void set_verbosity_str(std::string &verbosity, int meshing_verbosity)
{
    switch (meshing_verbosity) {
    case -1:
        verbosity = "Q";
        break;
    case 0:
        verbosity = "";
        break;
    case 1:
        verbosity = "V";
        break;
    case 2:
        verbosity = "VV";
        break;
    case 3:
        verbosity = "VVV";
        break;
    default:
        verbosity = "";
        break;
    }
}


void set_volume_str(std::string &vol, double max_volume)
{
    vol.clear();
    if (max_volume > 0) {
        vol += 'a';
        vol += std::to_string((long double)max_volume);
    }
    else if (max_volume == 0) {
        // max_volume is set inside regional attributes (usu. read from poly files)
        vol += 'a';
    }
}


void set_2d_quality_str(std::string &quality, double min_angle)
{
    quality.clear();
    if (min_angle > 0) {
        quality += 'q';
        quality += std::to_string((long double)min_angle);
    }
}


void triangulate_polygon
(double min_angle, double max_area,
 int meshing_verbosity,
 int npoints, int nsegments,
 const double *points, const int *segments, const int *segflags,
 const int nregions, const double *regionattributes,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags, double **outregattr)
{
    char options[255];
    triangulateio in, out;

    std::string verbosity, vol, quality;
    set_verbosity_str(verbosity, meshing_verbosity);
    set_volume_str(vol, max_area);
    set_2d_quality_str(quality, min_angle);

    if( nregions > 0 )
        std::sprintf(options, "%s%spjz%sA", verbosity.c_str(), quality.c_str(), vol.c_str());
    else
        std::sprintf(options, "%s%spjz%s", verbosity.c_str(), quality.c_str(), vol.c_str());

    if( meshing_verbosity >= 0 )
        std::cout << "The meshing option is: " << options << '\n';

    in.pointlist = const_cast<double*>(points);
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.numberofpoints = npoints;
    in.numberofpointattributes = 0;

    in.trianglelist = NULL;
    in.triangleattributelist = NULL;
    in.trianglearealist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 3;
    in.numberoftriangleattributes = 0;

    in.segmentlist = const_cast<int*>(segments);
    in.segmentmarkerlist = const_cast<int*>(segflags);
    in.numberofsegments = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.numberofregions = nregions;
    if( nregions > 0 )
        in.regionlist = const_cast<double*>(regionattributes);
    else
        in.regionlist = NULL;

    out.pointlist = NULL;
    out.pointattributelist = NULL;
    out.pointmarkerlist = NULL;
    out.trianglelist = NULL;
    out.triangleattributelist = NULL;
    out.neighborlist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.edgelist = NULL;
    out.edgemarkerlist = NULL;

    /*******************************/
    triangulate(options, &in, &out, NULL);
    /*******************************/

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;

    *ntriangles = out.numberoftriangles;
    *triangles = out.trianglelist;

    *noutsegments = out.numberofsegments;
    *outsegments = out.segmentlist;
    *outsegflags = out.segmentmarkerlist;
    *outregattr = out.triangleattributelist;

    trifree(out.pointmarkerlist);
}


void set_3d_quality_str(std::string &quality, double max_ratio,
                        double min_dihedral_angle, double max_dihedral_angle)
{
    quality.clear();
    if (max_ratio > 0) {
        quality += 'q';
        quality += std::to_string((long double)max_ratio);
        quality += "qq";
        quality += std::to_string((long double)min_dihedral_angle);
        quality += "qqq";
        quality += std::to_string((long double)max_dihedral_angle);
    }
}


#ifdef THREED
void tetrahedralize_polyhedron
(double max_ratio, double min_dihedral_angle, double max_volume,
 int vertex_per_polygon, int meshing_verbosity, int optlevel,
 int npoints, int nsegments,
 const double *points, const int *segments, const int *segflags,
 const tetgenio::facet *facets,
 const int nregions, const double *regionattributes,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags, double **outregattr)
{
    //
    // Setting Tetgen options.
    //
    char options[255];
    double max_dihedral_angle = 180 - 3 * min_dihedral_angle;

    std::string verbosity, vol, quality;
    set_verbosity_str(verbosity, meshing_verbosity);
    set_volume_str(vol, max_volume);
    set_3d_quality_str(quality, max_ratio, min_dihedral_angle, max_dihedral_angle);

    if( nregions > 0 )
        std::sprintf(options, "%s%s%spzs%dA", verbosity.c_str(), quality.c_str(), vol.c_str(), optlevel);
    else
        std::sprintf(options, "%s%s%spzs%d", verbosity.c_str(), quality.c_str(), vol.c_str(), optlevel);

    if( meshing_verbosity >= 0 )
        std::cout << "The meshing option is: " << options << '\n';

    //
    // Setting input arrays to tetgen
    //

    // NOTE: all tetgenio pointers will need to be
    // reset to NULL to prevent double-free error
    tetgenio in;
    in.pointlist = const_cast<double*>(points);
    in.numberofpoints = npoints;

    tetgenio::polygon *polys;
    tetgenio::facet *fl;
    if (facets == NULL) {
        polys = new tetgenio::polygon[nsegments];
        for (int i=0; i<nsegments; ++i) {
            polys[i].vertexlist = const_cast<int*>(&segments[i*vertex_per_polygon]);
            polys[i].numberofvertices = vertex_per_polygon;
        }

        fl = new tetgenio::facet[nsegments];
        for (int i=0; i<nsegments; ++i) {
            fl[i].polygonlist = &polys[i];
            fl[i].numberofpolygons = 1;
            fl[i].holelist = NULL;
            fl[i].numberofholes = 0;
        }
    }
    else {
        fl = const_cast<tetgenio::facet*>(facets);
    }

    in.facetlist = fl;
    in.facetmarkerlist = const_cast<int*>(segflags);
    in.numberoffacets = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.numberofregions = nregions;
    if( nregions > 0 )
        in.regionlist = const_cast<double*>(regionattributes);
    else
        in.regionlist = NULL;

    tetgenio out;
    /*******************************/
    tetrahedralize(options, &in, &out, NULL, NULL);
    /*******************************/

    // the destructor of tetgenio will free any non-NULL pointer
    // set in.pointers to NULL to prevent double-free
    in.pointlist = NULL;
    in.facetmarkerlist = NULL;
    in.facetlist = NULL;
    in.regionlist = NULL;
    if (facets == NULL) {
        delete [] polys;
        delete [] fl;
    }

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;
    out.pointlist = NULL;

    *ntriangles = out.numberoftetrahedra;
    *triangles = out.tetrahedronlist;
    out.tetrahedronlist = NULL;

    *noutsegments = out.numberoftrifaces;
    *outsegments = out.trifacelist;
    *outsegflags = out.trifacemarkerlist;
    *outregattr = out.tetrahedronattributelist;
    out.trifacelist = NULL;
    out.trifacemarkerlist = NULL;
    out.tetrahedronattributelist = NULL;

}
#endif


void points_to_mesh(const Param &param, Variables &var,
                    int npoints, const double *points,
                    int n_init_segments, const int *init_segments, const int *init_segflags,
                    int nregions, const double *regattr,
                    double max_elem_size, int vertex_per_polygon)
{
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    points_to_new_mesh(param.mesh, npoints, points,
                       n_init_segments, init_segments, init_segflags,
                       nregions, regattr,
                       max_elem_size, vertex_per_polygon,
                       var.nnode, var.nelem, var.nseg,
                       pcoord, pconnectivity, psegment, psegflag, pregattr);

    var.coord = new array_t(pcoord, var.nnode);
    var.connectivity = new conn_t(pconnectivity, var.nelem);
    var.segment = new segment_t(psegment, var.nseg);
    var.segflag = new segflag_t(psegflag, var.nseg);
    var.regattr = new regattr_t(pregattr, var.nelem);
}


void new_mesh_uniform_resolution(const Param& param, Variables& var)
{
    int npoints = 4 * (NDIMS - 1); // 2D:4;  3D:8
    double *points = new double[npoints*NDIMS];

    int n_init_segments = 2 * NDIMS; // 2D:4;  3D:6
    int n_segment_nodes = 2 * (NDIMS - 1); // 2D:2; 3D:4
    int *init_segments = new int[n_init_segments*n_segment_nodes];
    int *init_segflags = new int[n_init_segments];

    const int attr_ndata = NDIMS+2;
    const int nregions = (param.ic.mattype_option == 0) ? param.mat.nmat : 1;
    double *regattr = new double[nregions*attr_ndata]; // each region has (NDIMS+2) data fields: x, (y,) z, region marker (mat type), and volume.

    double elem_size;  // size of a typical element
    int vertex_per_polygon = 4;

#ifndef THREED
    {
	/* Define 4 corner points of the rectangle, with this order:
         *            BOUNDZ1
         *          0 ------- 3
         *  BOUNDX0 |         | BOUNDX1
         *          1 ------- 2
         *            BOUNDZ0
         */
        // corner 0
        points[0] = 0;
	points[1] = 0;
        // corner 1
	points[2] = 0;
	points[3] = -param.mesh.zlength;
        // corner 2
	points[4] = param.mesh.xlength;
	points[5] = -param.mesh.zlength;
        // corner 3
	points[6] = param.mesh.xlength;
	points[7] = 0;

	for (int i=0; i<n_init_segments; ++i) {
            // 0th node of the i-th segment
	    init_segments[2*i] = i;
            // 1st node of the i-th segment
	    init_segments[2*i+1] = i+1;
	}
	// the 1st node of the last segment is connected back to the 0th node
	init_segments[2*n_init_segments-1] = 0;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDZ0;
        init_segflags[2] = BOUNDX1;
        init_segflags[3] = BOUNDZ1;

        elem_size = 1.5 * param.mesh.resolution * param.mesh.resolution;

        // Need to modify when nregions > 1.
        for (int i = 0; i < nregions; i++) {
            regattr[i * attr_ndata] = 0.5*param.mesh.xlength;
            regattr[i * attr_ndata + 1] = -0.5*param.mesh.zlength;
            regattr[i * attr_ndata + 2] = 0;
            regattr[i * attr_ndata + 3] = -1;
        }
    }
#else
    {
	/* Define 8 corner points of the box, with this order:
         *         4 ------- 7
         *        /         /|
         *       /         / 6
         *      0 ------- 3 /
         *      |         |/
         *      1 ------- 2
         *
         * Cut-out diagram with boundary flag:
         *             4 ------- 7
         *             | BOUNDZ1 |
         *   4 ------- 0 ------- 3 ------- 7 ------- 4
         *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
         *   5 ------- 1 ------- 2 ------- 6 ------- 5
         *             | BOUNDZ0 |
         *             5 ------- 6
         */

        // corner 0
        points[0] = 0;
	points[1] = 0;
	points[2] = 0;
        // corner 1
	points[3] = 0;
	points[4] = 0;
	points[5] = -param.mesh.zlength;
        // corner 2
	points[6] = param.mesh.xlength;
	points[7] = 0;
	points[8] = -param.mesh.zlength;
        // corner 3
	points[9] = param.mesh.xlength;
	points[10] = 0;
	points[11] = 0;
        // corner 4
	points[12] = 0;
	points[13] = param.mesh.ylength;
	points[14] = 0;
        // corner 5
	points[15] = 0;
	points[16] = param.mesh.ylength;
	points[17] = -param.mesh.zlength;
        // corner 6
	points[18] = param.mesh.xlength;
	points[19] = param.mesh.ylength;
	points[20] = -param.mesh.zlength;
        // corner 7
	points[21] = param.mesh.xlength;
	points[22] = param.mesh.ylength;
	points[23] = 0;

        // BOUNDX0
        init_segments[0] = 0;
        init_segments[1] = 1;
        init_segments[2] = 5;
        init_segments[3] = 4;
        // BOUNDY0
        init_segments[4] = 0;
        init_segments[5] = 3;
        init_segments[6] = 2;
        init_segments[7] = 1;
        // BOUNDZ0
        init_segments[8] = 1;
        init_segments[9] = 2;
        init_segments[10] = 6;
        init_segments[11] = 5;
        // BOUNDX1
        init_segments[12] = 3;
        init_segments[13] = 7;
        init_segments[14] = 6;
        init_segments[15] = 2;
        // BOUNDY1
        init_segments[16] = 7;
        init_segments[17] = 4;
        init_segments[18] = 5;
        init_segments[19] = 6;
        // BOUNDZ1
        init_segments[20] = 0;
        init_segments[21] = 4;
        init_segments[22] = 7;
        init_segments[23] = 3;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDY0;
        init_segflags[2] = BOUNDZ0;
        init_segflags[3] = BOUNDX1;
        init_segflags[4] = BOUNDY1;
        init_segflags[5] = BOUNDZ1;

        elem_size = 0.7 * param.mesh.resolution
            * param.mesh.resolution * param.mesh.resolution;

        // Need to modify when nregions > 1. Using .poly file might be better.
        for (int i = 0; i < nregions; i++) {
            regattr[i * attr_ndata] = 0.5*param.mesh.xlength;
            regattr[i * attr_ndata + 1] = 0.5*param.mesh.ylength;
            regattr[i * attr_ndata + 2] = -0.5*param.mesh.zlength;
            regattr[i * attr_ndata + 3] = 0;
            regattr[i * attr_ndata + 4] = -1;
        }
    }
#endif

    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags, nregions, regattr,
                   elem_size, vertex_per_polygon);

    delete [] points;
    delete [] init_segments;
    delete [] init_segflags;
    delete [] regattr;
}


void new_mesh_refined_zone(const Param& param, Variables& var)
{
    const Mesh& m = param.mesh;

    // bounding box
    const double Lx = m.xlength;
#ifdef THREED
    const double Ly = m.ylength;
#endif
    const double Lz = m.zlength;

    // typical distance between nodes in the refined zone
    const double d = m.resolution;

    // adjust the bounds of the refined zone so that nodes are not on the boundary
    const double x0 = std::max(m.refined_zonex.first, d / Lx);
    const double x1 = std::min(m.refined_zonex.second, 1 - d / Lx);
#ifdef THREED
    const double y0 = std::max(m.refined_zoney.first, d / Ly);
    const double y1 = std::min(m.refined_zoney.second, 1 - d / Ly);
#endif
    const double z0 = std::max(m.refined_zonez.first, d / Lz);
    const double z1 = std::min(m.refined_zonez.second, 1 - d / Lz);

    const int npoints = 2 * 4 * (NDIMS - 1);
    double *points = new double[npoints*NDIMS];

    const int nsegments = 2 * (2 * NDIMS); // 2 x 2D:4;  3D:6
    const int nnodes = 2 * (NDIMS - 1);    // 2D:2; 3D:4
    const int vertex_per_polygon = 4;
    int *segments = new int[nsegments*nnodes];
    int *segflags = new int[nsegments];

#ifndef THREED
    {
      /* Define 4 corner points of the rectangle, with this order:
       *            BOUNDZ1
       *          0 ------- 3
       *  BOUNDX0 |         | BOUNDX1
       *          1 ------- 2
       *            BOUNDZ0
       */

      // outer box            // inner box
      // corner 0             // corner 0
      points[0] = 0.0;        points[8+0] = +x0*Lx;
      points[1] = 0.0;        points[8+1] = -z0*Lz;
      // corner 1             // corner 1
      points[2] = 0.0;        points[8+2] = +x0*Lx;
      points[3] = -Lz;        points[8+3] = -z1*Lz;
      // corner 2             // corner 2
      points[4] = +Lx;        points[8+4] = +x1*Lx;
      points[5] = -Lz;        points[8+5] = -z1*Lz;
      // corner 3             // corner 3
      points[6] = +Lx;        points[8+6] = +x1*Lx;
      points[7] = 0.0;        points[8+7] = -z0*Lz;

      // outer segments       // inner segments
      // BOUNDX0              // BOUNDX0
      segments[0] = 0;        segments[8+0] = 4+0;
      segments[1] = 1;        segments[8+1] = 4+1;
      // BOUNDZ0              // BOUNDZ0
      segments[2] = 1;        segments[8+2] = 4+1;
      segments[3] = 2;        segments[8+3] = 4+2;
      // BOUNDX1              // BOUNDX1
      segments[4] = 2;        segments[8+4] = 4+2;
      segments[5] = 3;        segments[8+5] = 4+3;
      // BOUNDZ1              // BOUNDZ1
      segments[6] = 3;        segments[8+6] = 4+3;
      segments[7] = 0;        segments[8+7] = 4+0;

      // outer boundary       // inner boundary
      segflags[0] = BOUNDX0;  segflags[4+0] = 0;
      segflags[1] = BOUNDZ0;  segflags[4+1] = 0;
      segflags[2] = BOUNDX1;  segflags[4+2] = 0;
      segflags[3] = BOUNDZ1;  segflags[4+3] = 0;
    }
#else
    {
      /* Define 8 corner points of the box, with this order:
       *         4 ------- 7
       *        /         /|
       *       /         / 6
       *      0 ------- 3 /
       *      |         |/
       *      1 ------- 2
       *
       * Cut-out diagram with boundary flag:
       *             4 ------- 7
       *             | BOUNDZ1 |
       *   4 ------- 0 ------- 3 ------- 7 ------- 4
       *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
       *   5 ------- 1 ------- 2 ------- 6 ------- 5
       *             | BOUNDZ0 |
       *             5 ------- 6
       */

      // outer box            // inner box
      // corner 0             // corner 0
      points[ 0] = 0.0;       points[24+ 0] = +x0*Lx;
      points[ 1] = 0.0;       points[24+ 1] = +y0*Ly;
      points[ 2] = 0.0;       points[24+ 2] = -z0*Lz;
      // corner 1             // corner 1
      points[ 3] = 0.0;       points[24+ 3] = +x0*Lx;
      points[ 4] = 0.0;       points[24+ 4] = +y0*Ly;
      points[ 5] = -Lz;       points[24+ 5] = -z1*Lz;
      // corner 2             // corner 2
      points[ 6] = +Lx;       points[24+ 6] = +x1*Lx;
      points[ 7] = 0.0;       points[24+ 7] = +y0*Ly;
      points[ 8] = -Lz;       points[24+ 8] = -z1*Lz;
      // corner 3             // corner 3
      points[ 9] = +Lx;       points[24+ 9] = +x1*Lx;
      points[10] = 0.0;       points[24+10] = +y0*Ly;
      points[11] = 0.0;       points[24+11] = -z0*Lz;
      // corner 4             // corner 4
      points[12] = 0.0;       points[24+12] = +x0*Lx;
      points[13] = +Ly;       points[24+13] = +y1*Ly;
      points[14] = 0.0;       points[24+14] = -z0*Lz;
      // corner 5             // corner 5
      points[15] = 0.0;       points[24+15] = +x0*Lx;
      points[16] = +Ly;       points[24+16] = +y1*Ly;
      points[17] = -Lz;       points[24+17] = -z1*Lz;
      // corner 6             // corner 6
      points[18] = +Lx;       points[24+18] = +x1*Lx;
      points[19] = +Ly;       points[24+19] = +y1*Ly;
      points[20] = -Lz;       points[24+20] = -z1*Lz;
      // corner 7             // corner 7
      points[21] = +Lx;       points[24+21] = +x1*Lx;
      points[22] = +Ly;       points[24+22] = +y1*Ly;
      points[23] = 0.0;       points[24+23] = -z0*Lz;

      // outer segments       // inner segments
      // BOUNDX0              // BOUNDX0
      segments[ 0] = 0;       segments[24+ 0] = 8+0;
      segments[ 1] = 1;       segments[24+ 1] = 8+1;
      segments[ 2] = 5;       segments[24+ 2] = 8+5;
      segments[ 3] = 4;       segments[24+ 3] = 8+4;
      // BOUNDY0              // BOUNDY0
      segments[ 4] = 0;       segments[24+ 4] = 8+0;
      segments[ 5] = 3;       segments[24+ 5] = 8+3;
      segments[ 6] = 2;       segments[24+ 6] = 8+2;
      segments[ 7] = 1;       segments[24+ 7] = 8+1;
      // BOUNDZ0              // BOUNDZ0
      segments[ 8] = 1;       segments[24+ 8] = 8+1;
      segments[ 9] = 2;       segments[24+ 9] = 8+2;
      segments[10] = 6;       segments[24+10] = 8+6;
      segments[11] = 5;       segments[24+11] = 8+5;
      // BOUNDX1              // BOUNDX1
      segments[12] = 3;       segments[24+12] = 8+3;
      segments[13] = 7;       segments[24+13] = 8+7;
      segments[14] = 6;       segments[24+14] = 8+6;
      segments[15] = 2;       segments[24+15] = 8+2;
      // BOUNDY1              // BOUNDY1
      segments[16] = 7;       segments[24+16] = 8+7;
      segments[17] = 4;       segments[24+17] = 8+4;
      segments[18] = 5;       segments[24+18] = 8+5;
      segments[19] = 6;       segments[24+19] = 8+6;
      // BOUNDZ1              // BOUNDZ1
      segments[20] = 0;       segments[24+20] = 8+0;
      segments[21] = 4;       segments[24+21] = 8+4;
      segments[22] = 7;       segments[24+22] = 8+7;
      segments[23] = 3;       segments[24+23] = 8+3;

      // outer boundary       // inner boundary
      segflags[0] = BOUNDX0;  segflags[6+0] = 0;
      segflags[1] = BOUNDY0;  segflags[6+1] = 0;
      segflags[2] = BOUNDZ0;  segflags[6+2] = 0;
      segflags[3] = BOUNDX1;  segflags[6+3] = 0;
      segflags[4] = BOUNDY1;  segflags[6+4] = 0;
      segflags[5] = BOUNDZ1;  segflags[6+5] = 0;
    }
#endif

    const int nregions = 2;
    const int ndata = NDIMS+2; // each region has (NDIMS+2) data fields: x, (y,) z, region marker (mat type), area/volume.
    double *regattr = new double[nregions*ndata];
    double *region0 = regattr;
    double *region1 = regattr + ndata;
#ifndef THREED
    const double area1 = 1.5 * (d*d);
#else
    const double area1 = 0.7 * (d*d*d);
#endif
    const double area0 = area1 * m.largest_size;
    int pos = 0;
    region0[pos] = +d/2;  region1[pos++] = +x0*Lx+d/2; // x
#if THREED
    region0[pos] = +d/2;  region1[pos++] = +y0*Ly+d/2; // y
#endif
    region0[pos] = -d/2;  region1[pos++] = -z0*Lz-d/2; // z
    region0[pos] = 0.0;   region1[pos++] = 0.0;        // attribute
    region0[pos] = area0; region1[pos++] = area1;      // area/volume

    const double max_elem_size = 0;

    points_to_mesh(param, var, npoints, points,
                   nsegments, segments, segflags, nregions, regattr,
                   max_elem_size, vertex_per_polygon);

    delete [] points;
    delete [] segments;
    delete [] segflags;
    delete [] regattr;
}


void my_fgets(char *buffer, std::size_t size, std::FILE *fp,
              int &lineno, const std::string &filename)
{
    char *s;
    while (1) {
        ++ lineno;
        s = std::fgets(buffer, size, fp);
        if (! s) {
            std::cerr << "Error: reading line " << lineno
                      << " of '" << filename << "'\n";
            std::exit(2);
        }
        if (std::strlen(buffer) == size-1 && buffer[size-2] != '\n') {
            std::cerr << "Error: reading line " << lineno
                      << " of '" << filename << "', line is too long.\n";
            std::exit(2);
        }

        // check for blank lines and comments
        if (buffer[0] != '\n' && buffer[0] != '#') break;
    }
}


void new_mesh_from_polyfile(const Param& param, Variables& var)
{
    /* The format specifiction for the poly file can be found in:
     *   2D:  http://www.cs.cmu.edu/~quake/triangle.poly.html
     *   3D:  http://wias-berlin.de/software/tetgen/fformats.poly.html
     *
     * Note:
     * 1. the poly file in 3D has a more complicated format.
     * 2. The last row in the regional attributes is the regional constraint
     *    on the maximum triangle area or tetrahedra volume, which is the max
     *    size of an element.
     *    when meshing_option=90, the unit is meter^dim.
     *    when meshing_option=91, the unit is mesh.resolution^dim.
     */

#ifdef THREED
    const double std_elem_size = 0.7 * param.mesh.resolution
        * param.mesh.resolution * param.mesh.resolution;
#else
    const double std_elem_size = 1.5 * param.mesh.resolution * param.mesh.resolution;
#endif

    std::FILE *fp = std::fopen(param.mesh.poly_filename.c_str(), "r");
    if (! fp) {
        std::cerr << "Error: Cannot open poly_filename '" << param.mesh.poly_filename << "'\n";
        std::exit(2);
    }

    int lineno = 0;
    int n;
    char buffer[2550];

    // get header of node list
    int npoints;
    {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        int dim, nattr, nbdrym;
        n = std::sscanf(buffer, "%d %d %d %d", &npoints, &dim, &nattr, &nbdrym);
        if (n != 4) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }

        if (dim != NDIMS ||
            nattr != 0 ||
            nbdrym != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get node list
    double *points = new double[npoints * NDIMS];
    for (int i=0; i<npoints; i++) {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        int k;
        double *x = &points[i*NDIMS];
#ifdef THREED
        n = std::sscanf(buffer, "%d %lf %lf %lf", &k, x, x+1, x+2);
#else
        n = std::sscanf(buffer, "%d %lf %lf", &k, x, x+1);
#endif
        if (n != NDIMS+1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (k != i) {
            std::cerr << "Error: node number is continuous from 0 at line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get header of segment (facet) list
    int n_init_segments;
    {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        int has_bdryflag;
        n = std::sscanf(buffer, "%d %d", &n_init_segments, &has_bdryflag);
        if (n != 2) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }

        if (has_bdryflag != 1) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get segment (facet) list
#ifdef THREED
    auto facets = new tetgenio::facet[n_init_segments];
    int *init_segflags = new int[n_init_segments];
    for (int i=0; i<n_init_segments; i++) {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        auto &f = facets[i];
        int npolygons, nholes, bdryflag;
        n = std::sscanf(buffer, "%d %d %d", &npolygons, &nholes, &bdryflag);
        if (n != 3) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (npolygons <= 0 ||
            nholes != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (bdryflag == 0) goto flag_ok;
        for (int j=0; j<nbdrytypes; j++) {
            if (bdryflag == 1 << j) goto flag_ok;
        }
        std::cerr << "Error: bdry_flag has multiple bits set in line " << lineno
                  << " of '" << param.mesh.poly_filename << "'\n";
        std::exit(1);
    flag_ok:
        init_segflags[i] = bdryflag;

        f.polygonlist = new tetgenio::polygon[npolygons];
        f.numberofpolygons = npolygons;
        f.holelist = NULL;
        f.numberofholes = 0;

        for (int j=0; j<npolygons; j++) {
            my_fgets(buffer, 4096, fp, lineno, param.mesh.poly_filename);

            std::istringstream inbuf(std::string(buffer, 4096));
            int nvertex;
            inbuf >> nvertex;
            if (nvertex < NODES_PER_FACET || nvertex > 9999) {
                std::cerr << "Error: unsupported number of polygon points in line " << lineno
                          << " of '" << param.mesh.poly_filename << "'\n";
                std::exit(1);
            }

            f.polygonlist[j].vertexlist = new int[nvertex];
            f.polygonlist[j].numberofvertices = nvertex;
            for (int k=0; k<nvertex; k++) {
                inbuf >> f.polygonlist[j].vertexlist[k];
                if (f.polygonlist[j].vertexlist[k] < 0 ||
                    f.polygonlist[j].vertexlist[k] >= npoints) {
                    std::cerr << "Error: segment contains out-of-range node # [0-" << npoints
                              <<"] in line " << lineno << " of '"
                              << param.mesh.poly_filename << "'\n";
                    std::exit(1);
                }
            }
        }
    }
#else
    int *init_segments = new int[n_init_segments * NODES_PER_FACET];
    int *init_segflags = new int[n_init_segments];
    for (int i=0; i<n_init_segments; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int *x = &init_segments[i*NODES_PER_FACET];
        int junk, bdryflag;
        n = std::sscanf(buffer, "%d %d %d %d", &junk, x, x+1, &bdryflag);
        if (n != NODES_PER_FACET+2) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (bdryflag == 0) goto flag_ok;
        for (int j=0; j<nbdrytypes; j++) {
            if (bdryflag == 1 << j) goto flag_ok;
        }
        std::cerr << "Error: bdry_flag has multiple bits set in line " << lineno
                  << " of '" << param.mesh.poly_filename << "'\n";
        std::exit(1);
    flag_ok:
        init_segflags[i] = bdryflag;
    }
    for (int i=0; i<n_init_segments; i++) {
        int *x = &init_segments[i*NODES_PER_FACET];
        for (int j=0; j<NODES_PER_FACET; j++) {
            if (x[j] < 0 || x[j] >= npoints) {
                std::cerr << "Error: segment contains out-of-range node # [0-" << npoints
                          <<"] in line " << lineno << " of '"
                          << param.mesh.poly_filename << "'\n";
                std::exit(1);
            }
        }
    }
#endif

    // get header of hole list
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int nholes;
        n = std::sscanf(buffer, "%d", &nholes);
        if (n != 1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }

        if (nholes != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get header of region list
    int nregions;
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        n = std::sscanf(buffer, "%d", &nregions);
        if (n != 1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
        if (nregions <= 0) {
            std::cerr << "Error: nregions <= 0, at line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            std::exit(1);
        }
    }

    // get region list
    double *regattr = new double[nregions * (NDIMS+2)]; // each region has these data fields: x, (y,) z, region marker (mattype), and volume.
    bool has_max_size = false;
    for (int i=0; i<nregions; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int junk;
        double *x = &regattr[i*(NDIMS+2)];
#ifdef THREED
        n = std::sscanf(buffer, "%d %lf %lf %lf %lf %lf", &junk, x, x+1, x+2, x+3, x+4);
#else
        n = std::sscanf(buffer, "%d %lf %lf %lf %lf", &junk, x, x+1, x+2, x+3);
#endif
        if (n != NDIMS+3) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'. "<<NDIMS+3<<" values should be given but only "<<n<<" found.\n";
            std::exit(1);
        }

        if ( x[NDIMS] < 0 || x[NDIMS] >= param.mat.nmat ) {
            std::cerr << "Error: "<<NDIMS+2<<"-th value in line "<<lineno<<" should be >=0 and < "<<param.mat.nmat<<" (=mat.num_materials) but is "<<x[NDIMS]<<"\n";
            std::cerr << "Note that this parameter is directly used as the index of mat. prop. arrays.\n";
            std::exit(1);
        }

        if ( x[NDIMS+1] > 0 ) {
            has_max_size = true; // max area is set for this region

            // normalize the element size w.r.t mesh.resolution
            if (param.mesh.meshing_option == 91) {
                x[NDIMS+1] *= std_elem_size;
            }
        }
    }
    double max_elem_size = std_elem_size;
    if ( has_max_size ) max_elem_size = 0; // special value, see set_volume_str() above.

#ifdef THREED
    // shortcutting points_to_mesh() for polygon boundary
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    tetrahedralize_polyhedron(param.mesh.max_ratio,
                              param.mesh.min_tet_angle, max_elem_size,
                              0,
                              param.mesh.meshing_verbosity,
                              param.mesh.tetgen_optlevel,
                              npoints, n_init_segments, points,
                              NULL, init_segflags,
                              facets,
                              nregions, regattr,
                              &var.nnode, &var.nelem, &var.nseg,
                              &pcoord, &pconnectivity,
                              &psegment, &psegflag, &pregattr);

    var.coord = new array_t(pcoord, var.nnode);
    var.connectivity = new conn_t(pconnectivity, var.nelem);
    var.segment = new segment_t(psegment, var.nseg);
    var.segflag = new segflag_t(psegflag, var.nseg);
    var.regattr = new regattr_t(pregattr, var.nelem);
#else
    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags, nregions, regattr,
                   max_elem_size, NODES_PER_FACET);
#endif

    delete [] points;
#ifdef THREED
    for (int i=0; i<n_init_segments; i++) {
        auto f = facets[i];
        for (int j=0; j<f.numberofpolygons; j++) {
            delete [] f.polygonlist[j].vertexlist;
        }
        delete [] f.polygonlist;
    }
    delete [] facets;
#else
    delete [] init_segments;
#endif
    delete [] init_segflags;
    delete [] regattr;
}

#ifdef USEEXODUS
void new_mesh_from_exofile(const Param& param, Variables& var)
{

#ifndef THREED
    std::cerr << "Error: Importing an exofile currently works in 3D only.\n";
    std::exit(2);
#endif

    // Open an .exo file.
    int CPU_word_size = 0;
    int IO_word_size = 0;
    float version = 0.0;
    int exoid = ex_open( param.mesh.exo_filename.c_str(), /* filename path */
                         EX_READ, /* access mode = READ */
                         &CPU_word_size, /* CPU word size */
                         &IO_word_size, /* IO word size */
                         &version); /* ExodusII library version */
    if (exoid < 0) {
        std::cerr << "Error: Cannot open exo_filename '" << param.mesh.exo_filename << "'\n";
        std::exit(2);
    }

    // Read database parameters.
    int error;
    char title[MAX_LINE_LENGTH+1];
    int num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;

    error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem,
                         &num_elem_blk, &num_node_sets, &num_side_sets);
    if( error != 0 ) {
        std::cerr << "Error: Unable to read database parameters from '" << param.mesh.exo_filename << std::endl;
        std::exit(2);
    }
    var.nnode = num_nodes;
    var.nelem = num_elem;
    std::cerr <<" Reading " << param.mesh.exo_filename <<"."<<std::endl;
    std::cerr <<" Numbers of nodes and elements are " << var.nnode <<" and "<< var.nelem <<"."<<std::endl;
    std::cerr <<" Number of element blocks is " << num_elem_blk <<"."<<std::endl;
    std::cerr <<" Number of node sets is " << num_node_sets <<"."<<std::endl;
    std::cerr <<" Number of side sets is " << num_side_sets <<"."<<std::endl;
    if( param.mat.nmat != num_elem_blk) {
        std::cerr <<"param.mat.nmat is not equal to # of element blocks in this exo file!"<<std::endl;
        std::cerr <<"Check if your material parameters are properly set!!"<<std::endl;
        std::exit(2);
    }

    // Assign node coordinates.
    float *x = (float *) calloc(var.nnode, sizeof(float));
    float *y = (float *) calloc(var.nnode, sizeof(float));
    float *z = (float *) calloc(var.nnode, sizeof(float));

    error = ex_get_coord (exoid, x, y, z);
    if( error != 0 ) {
        std::cerr << "Error: Unable to read coordinates from '" << param.mesh.exo_filename << "'\n";
        std::exit(2);
    }

    // assign node coordinates to var.coord.
    var.coord = new array_t(var.nnode);
    double* coord = var.coord->data();
    for(int i=0; i<var.nnode; i++) {
        coord[i*NDIMS]   = static_cast<double>(x[i]);
        coord[i*NDIMS+1] = static_cast<double>(y[i]);
        coord[i*NDIMS+2] = static_cast<double>(z[i]);
    }
    free(x);
    free(y);
    free(z);

    // Process element blocks.
    //
    // Note that blocks are the basic units for
    // organizing connectivity and material id.
    //
    int *ids = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_elem_in_block = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_nodes_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_edges_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_faces_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_attr = (int *) calloc(num_elem_blk, sizeof(int));
    char elem_type[MAX_STR_LENGTH+1];

    // - Read element block ids.
    error = ex_get_ids (exoid, EX_ELEM_BLOCK, ids);
    if( error != 0 ) {
        std::cerr << "Error: Unable to get element block ids." << std::endl;
        std::exit(2);
    }
    // - Read element block parameters.
    for (int i=0; i<num_elem_blk; i++) {
        error = ex_get_block (exoid, EX_ELEM_BLOCK, ids[i], elem_type,
                                &(num_elem_in_block[i]), &(num_nodes_per_elem[i]),
                                &(num_edges_per_elem[i]), &(num_faces_per_elem[i]), 
                                &(num_attr[i]));
        if( error != 0 ) {
            std::cerr << "Error: Unable to element block " << ids[i] ;
            std::cerr << " out of " << num_elem_blk << " blocks." << std::endl;
            std::exit(2);
        }
        if( NODES_PER_ELEM != num_nodes_per_elem[i] ) {
            std::cerr << "Error: Element has " << num_nodes_per_elem[i] << " nodes per element but should have "<<NODES_PER_ELEM<<" because element type should be uniformly tetrahedral."<< std::endl;
            std::exit(2);
        }
    }

    // Create the array of region attributes.
    // Element blocks must be numbered to be consistent
    // with the following simplistic mapping between block id and mat. id:
    //
    // block id: 1 2 3 4 ...
    // mat.  id: 0 1 2 3 ...
    //
    // Will need a more general mapping from block id to mat id.
    // One way is to provide a map from block names to mat id. TBD.
    //
    var.regattr = new regattr_t(var.nelem);
    double *attr = var.regattr->data();
    int start = 0;
    for (int i=0; i<num_elem_blk; i++) {
        for(int j=0; j<num_elem_in_block[i]; j++)
            attr[start+j] = static_cast<double>(ids[i]-1);
        start += num_elem_in_block[i];
    }

    // - Read and append element connectivity
    int* connect[num_elem_blk];
    for (int i=0; i<num_elem_blk; i++) {
        connect[i] = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), sizeof(int));
        error = ex_get_conn (exoid, EX_ELEM_BLOCK, ids[i], connect[i], 0, 0);
        if( error != 0 ) {
            std::cerr << "Error: Unable to connectivity for element block " << ids[i] ;
            std::cerr << " out of " << num_elem_blk << " blocks." << std::endl;
            std::exit(2);
        }
    }

    { // To render 'int *conn' local to this block.
        var.connectivity = new conn_t(var.nelem);
        int *conn = var.connectivity->data();
        start = 0;
        for (int i=0; i<num_elem_blk; i++) {
            for(int j=0; j<num_elem_in_block[i]; j++) {
                const int elem_num = start + j;
                for(int k=0; k < NODES_PER_ELEM; k++ )
                    conn[ NODES_PER_ELEM*elem_num + k ] = connect[i][ NODES_PER_ELEM*j + k ]-1;
            }
            // can the j and k loops be replaced with memcpy?
            // std::memcpy( &(conn[NODES_PER_ELEM*start]), &(connect[i]), num_elem_in_block[i]*NODES_PER_ELEM*sizeof(int) );
            start += num_elem_in_block[i];
        }
    } // code block to localize conn

    for (int i=0; i<num_elem_blk; i++) free(connect[i]);
    free (ids);
    free (num_elem_in_block);
    free (num_nodes_per_elem);
    free (num_edges_per_elem);
    free (num_faces_per_elem);
    free (num_attr);

    //
    // Read side sets (i.e., boundaries) one by one.
    //
    ids = (int *) calloc(num_side_sets, sizeof(int));
    int *num_sides_in_set = (int *) calloc(num_side_sets, sizeof(int));
    int *num_df_in_set = (int *) calloc(num_side_sets, sizeof(int));

    // - Read side set ids.
    error = ex_get_ids (exoid, EX_SIDE_SET, ids);
    if( error != 0 ) {
            std::cerr << "Error: Unable to get side set ids." << std::endl;
            std::exit(2);
    }
    var.nseg = 0;
    for(int i=0; i<num_side_sets; i++) {
        error = ex_get_set_param (exoid, EX_SIDE_SET, ids[i], &(num_sides_in_set[i]),
                                       &(num_df_in_set[i]));
        if( error != 0 ) {
            std::cerr << "Error: Unable to read "<<i<<"-th side set parameters." << std::endl;
            std::exit(2);
        }
        var.nseg += num_sides_in_set[i];
    }
    std::cerr <<" Numbers of segments are " << var.nseg <<"."<<std::endl;

    // list of elements in the side set
    int *elem_list[num_side_sets];
    // list of facets (sides) of a corresponding element in the side set
    int *side_list[num_side_sets];
    // list of the number of nodes of a facet. Uniformly 3 in DES3D.
    int *node_cnt_list[num_side_sets];
    // list of node IDs of the facet-composing nodes
    int *node_list[num_side_sets];
    // list of distribution factors
    float *dist_fact[num_side_sets];

    for(int i=0; i<num_side_sets; i++) {
        elem_list[i]     = (int *) calloc(num_sides_in_set[i], sizeof(int));
        side_list[i]     = (int *) calloc(num_sides_in_set[i], sizeof(int));
        node_cnt_list[i] = (int *) calloc(num_sides_in_set[i], sizeof(int));
        node_list[i]     = (int *) calloc(num_sides_in_set[i]*FACETS_PER_ELEM, sizeof(int));
        dist_fact[i]     = (float *) calloc(num_df_in_set[i], sizeof(float));

        // list of elements of which facet belongs to a side set
        // Note: The # of elements is same as # of sides!
        error = ex_get_set (exoid, EX_SIDE_SET, ids[i], elem_list[i], side_list[i]);
        if( error != 0 ) {
            std::cerr << "Error: Unable to read "<< i <<"-th side set." << std::endl;
            std::exit(2);
        }
        if (*num_df_in_set > 0) {
            error = ex_get_set_dist_fact (exoid, EX_SIDE_SET, ids[i], dist_fact[i]);
            if( error != 0 ) {
                std::cerr << "Error: Unable to read "<< i;
                std::cerr <<"-th side set's distribution factor list." << std::endl;
                std::exit(2);
            }
        }
    }

    // Assign memory to segments (boundary facets)
    var.segment = new segment_t(var.nseg, 0);
    int *segments = var.segment->data();

    // Assign memory to segment flags (which boundary a segments belong to).
    var.segflag = new segflag_t(var.nseg, 0);
    int *segflags = var.segflag->data();

    // Populate segment and segflag arrays
    // For sideset node ordering,
    // refer Table 4.2 on p. 30, "Exodus: A finite element data model" by Sjaardema et al.
    // Retrieved on 2019/09/28 from gsjaardema.github.io/seacas/exodusII-new.pdf.
    std::vector< std::vector<int> > local_node_list{{1,2,4},{2,3,4},{1,4,3},{1,3,2}};
    start = 0;
    const int *conn = var.connectivity->data();
    for (int i=0; i<num_side_sets; i++) {
        for (int j=0; j<num_sides_in_set[i]; j++) {
            const int elem_num = elem_list[i][j] - 1;
            const int side_num = side_list[i][j] - 1;
            for (int k=0; k<NODES_PER_FACET; k++) { // 3 nodes per triangular facet
                const int local_node_number = local_node_list[side_num][k] - 1;
                segments[ (start + j)*NODES_PER_FACET + k] = conn[ elem_num*NODES_PER_ELEM + local_node_number ];
            }
            segflags[start + j] = ids[i];
        }
        start += num_sides_in_set[i];
    }

    // free up memory.
    for (int i=0; i<num_side_sets; i++) {
        free (elem_list[i]);
        free (side_list[i]);
        free (node_cnt_list[i]);
        free (node_list[i]);
        free (dist_fact[i]);
    }
    free (num_sides_in_set);
    free (num_df_in_set);
    free (ids);

} // end of function new_mesh_from_exofile
#endif // end of USEEXODUS

} // end of anonymous namespace


void points_to_new_mesh(const Mesh &mesh, int npoints, const double *points,
                        int n_init_segments, const int *init_segments, const int *init_segflags,
                        int n_regions, const double *regattr,
                        double max_elem_size, int vertex_per_polygon,
                        int &nnode, int &nelem, int &nseg, double *&pcoord,
                        int *&pconnectivity, int *&psegment, int *&psegflag, double *&pregattr)
{
#ifdef THREED

    tetrahedralize_polyhedron(mesh.max_ratio,
                              mesh.min_tet_angle, max_elem_size,
                              vertex_per_polygon,
                              mesh.meshing_verbosity,
                              mesh.tetgen_optlevel,
                              npoints, n_init_segments, points,
                              init_segments, init_segflags,
                              NULL,
                              n_regions, regattr,
                              &nnode, &nelem, &nseg,
                              &pcoord, &pconnectivity,
                              &psegment, &psegflag, &pregattr);

#else

    triangulate_polygon(mesh.min_angle, max_elem_size,
                        mesh.meshing_verbosity,
                        npoints, n_init_segments, points,
                        init_segments, init_segflags,
                        n_regions, regattr,
                        &nnode, &nelem, &nseg,
                        &pcoord, &pconnectivity,
                        &psegment, &psegflag, &pregattr);

#endif

    if (nelem <= 0) {
#ifdef THREED
        std::cerr << "Error: tetrahedralization failed\n";
#else
        std::cerr << "Error: triangulation failed\n";
#endif
        std::exit(10);
    }

}


void points_to_new_surface(const Mesh &mesh, int npoints, const double *points,
                           int n_init_segments, const int *init_segments, const int *init_segflags,
                           int n_regions, const double *regattr,
                           double max_elem_size, int vertex_per_polygon,
                           int &nnode, int &nelem, int &nseg, double *&pcoord,
                           int *&pconnectivity, int *&psegment, int *&psegflag, double *&pregattr)
{
#ifdef THREED
    /* For triangulation of boundary surfaces in 3D */

    triangulate_polygon(mesh.min_angle, max_elem_size,
                        mesh.meshing_verbosity,
                        npoints, n_init_segments, points,
                        init_segments, init_segflags,
                        n_regions, regattr,
                        &nnode, &nelem, &nseg,
                        &pcoord, &pconnectivity,
                        &psegment, &psegflag, &pregattr);

    if (nelem <= 0) {
        std::cerr << "Error: surface triangulation failed\n";
        std::exit(10);
    }
#endif
}


void discard_internal_segments(int &nseg, segment_t &segment, segflag_t &segflag)
{
    int n = 0;
    while (n < nseg) {
        if (segflag[n][0] & BOUND_ANY) {
            // This is a boundary segment, go to next.
            n++;
        }
        else {
            // This is an internal boundary segment, replace it with the last segment.
            nseg--;
            segflag[n][0] = segflag[nseg][0];
            for (int i=0; i<NDIMS; i++) {
                segment[n][i] = segment[nseg][i];
            }
        }
    }

    segment.resize(nseg);
    segflag.resize(nseg);

}


void renumbering_mesh(const Param& param, array_t &coord, conn_t &connectivity,
                      segment_t &segment, regattr_t *regattr)
/* Note: the last argument is a pointer to regional attribute. If the pointer
 * is not NULL (when the initial mesh is created, before creating markers),
 * regattr is renumbered.
 */
{
    /* Renumbering nodes and elements to enhance cache coherance and better parallel performace. */

    const int nnode = coord.size();
    const int nelem = connectivity.size();
    const int nseg = segment.size();

    //
    // find the longest/shortest dimension
    //
    double_vec lengths = {param.mesh.xlength,
#ifdef THREED
                          param.mesh.ylength,
#endif
                          param.mesh.zlength};
    std::vector<std::size_t> idx(NDIMS);
    sortindex(lengths, idx);
    const int dmin = idx[0];
    const int dmid = idx[1];
    const int dmax = idx[NDIMS-1];

    //
    // sort coordinate of nodes and element centers
    //
    double_vec wn(nnode);
    const double f = 1e-3;
    for(int i=0; i<nnode; i++) {
        wn[i] = coord[i][dmax]
#ifdef THREED
            + f * coord[i][dmid]
#endif
            + f * f * coord[i][dmin];
    }

    double_vec we(nelem);
    for(int i=0; i<nelem; i++) {
        const int *conn = connectivity[i];
        we[i] = wn[conn[0]] + wn[conn[1]]
#ifdef THREED
            + wn[conn[2]]
#endif
            + wn[conn[NODES_PER_ELEM-1]];
    }

    // arrays to store the result of sorting
    std::vector<int> nd_idx(nnode);
    std::vector<int> el_idx(nelem);
    sortindex(wn, nd_idx);
    sortindex(we, el_idx);

    // inverse permutation
    std::vector<int> nd_inv(nnode);
    for(int i=0; i<nnode; i++)
      nd_inv[nd_idx[i]] = i;

    //
    // renumbering
    //
    array_t coord2(nnode);
    for(int i=0; i<nnode; i++) {
        int n = nd_idx[i];
        for(int j=0; j<NDIMS; j++)
            coord2[i][j] = coord[n][j];
    }
    coord.steal_ref(coord2);

    conn_t conn2(nelem);
    for(int i=0; i<nelem; i++) {
        int n = el_idx[i];
        for(int j=0; j<NODES_PER_ELEM; j++) {
            int k = connectivity[n][j];
            conn2[i][j] = nd_inv[k];
        }
    }
    connectivity.steal_ref(conn2);

    segment_t seg2(nseg);
    for(int i=0; i<nseg; i++) {
        for(int j=0; j<NDIMS; j++) {
            int k = segment[i][j];
            seg2[i][j] = nd_inv[k];
        }
    }
    segment.steal_ref(seg2);

    if (regattr != NULL) {
        regattr_t regattr2(nelem);
        for(int i=0; i<nelem; i++) {
            int n = el_idx[i];
            regattr2[i][0] = (*regattr)[n][0];
        }
        regattr->steal_ref(regattr2);
    }
}


void create_boundary_flags2(uint_vec &bcflag, int nseg,
                            const int *psegment, const int *psegflag)
{
    for (int i=0; i<nseg; ++i) {
        uint flag = static_cast<uint>(psegflag[i]);
        const int *n = psegment + i * NODES_PER_FACET;
        for (int j=0; j<NODES_PER_FACET; ++j) {
            bcflag[n[j]] |= flag;
        }
    }
}


void create_boundary_flags(Variables& var)
{
    // allocate and init to 0
    if (var.bcflag) delete var.bcflag;
    var.bcflag = new uint_vec(var.nnode);

    create_boundary_flags2(*var.bcflag, var.segment->size(),
                           var.segment->data(), var.segflag->data());
}


void create_boundary_nodes(Variables& var)
{
    /* var.bnodes[i] contains a list of nodes on the i-th boundary.
     * (See constants.hpp for the order of boundaries.)
     */
    for (std::size_t i=0; i<var.bcflag->size(); ++i) {
        uint f = (*var.bcflag)[i];
        for (int j=0; j<nbdrytypes; ++j) {
            if (f & (1<<j)) {
                // this node belongs to a boundary
                (var.bnodes[j]).push_back(i);
            }
        }
    }

    // for (int j=0; j<nbdrytypes; ++j) {
    //     std::cout << "boundary " << j << '\n';
    //     print(std::cout, var.bnodes[j]);
    //     std::cout << '\n';
    // }
}


namespace {

    struct OrderedInt
    {
#ifdef THREED
        int a, b, c;
        OrderedInt(int x, int y, int z)
        {
            // ensure that a < b < c
            if (x < y) {
                if (y < z)
                    a = x, b = y, c = z;
                else {
                    if (x < z)
                        a = x, b = z, c = y;
                    else
                        a = z, b = x, c = y;
                }
            }
            else {
                if (x < z)
                    a = y, b = x, c = z;
                else {
                    if (y < z)
                        a = y, b = z, c = x;
                    else
                        a = z, b = y, c = x;
                }
            }
        }

        bool operator==(OrderedInt &rhs)
        {
            return a==rhs.a && b==rhs.b && c==rhs.c;
        }

#else

        int a, b;
        OrderedInt(int x, int y)
        {
            // ensure that a < b
            if (x < y)
                a = x, b = y;
            else
                a = y, b = x;
        }

        bool operator==(OrderedInt &rhs)
        {
            return a==rhs.a && b==rhs.b;
        }
#endif
    };
}


void create_boundary_facets(Variables& var)
{
    /* var.bfacets[i] contains a list of facets (or segments in 2D)
     * on the i-th boundary. (See constants.hpp for the order of boundaries.)
     */

    // Looping through var.segment
    for (int i=0; i<var.nseg; ++i) {
        uint flag = static_cast<uint>((*var.segflag)[i][0]);
        if ((flag & BOUND_ANY) == 0) continue; // not a boundary facet
        // Nodes of this facet
        OrderedInt af((*var.segment)[i][0], (*var.segment)[i][1]
#ifdef THREED
                      , (*var.segment)[i][2]
#endif
                      );

        // Finding the corresponding element and facet #
        for (int e=0; e<var.nelem; ++e) {
            const int *conn = (*var.connectivity)[e];
            for (int f=0; f<FACETS_PER_ELEM; ++f) {
                if ((flag & (*var.bcflag)[conn[NODE_OF_FACET[f][0]]]
                    & (*var.bcflag)[conn[NODE_OF_FACET[f][1]]]
#ifdef THREED
                    & (*var.bcflag)[conn[NODE_OF_FACET[f][2]]]
#endif
                     ) == 0U) continue; // skip

                OrderedInt bf(conn[NODE_OF_FACET[f][0]], conn[NODE_OF_FACET[f][1]]
#ifdef THREED
                              , conn[NODE_OF_FACET[f][2]]
#endif
                              );
                if (af == bf) {
                    for (int k=0; k<nbdrytypes; ++k) {
                        if (flag == (1U << k)) {
                            var.bfacets[k].push_back(std::make_pair(e,f));
                            goto found_facet; // break out of nested loops
                        }
                    }
                }
            }
        }
        // not found
        std::cerr << "Error: " << i << "-th segment is not on any element\n";
        std::exit(12);

    found_facet:
        continue;
    }

    // for (int n=0; n<nbdrytypes; ++n) {
    //     std::cout << "boundary facet " << n << ":\n";
    //     print(std::cout, var.bfacets[n]);
    //     std::cout << '\n';
    //     for (int i=0; i<var.bfacets[n].size(); ++i) {
    //         int e = var.bfacets[n][i].first;
    //         int f = var.bfacets[n][i].second;
    //         const int *conn = (*var.connectivity)[e];
    //         std::cout << i << ", " << e << ":";
    //         for (int j=0; j<NODES_PER_FACET; ++j) {
    //             std::cout << " " << conn[NODE_OF_FACET[f][j]];
    //         }
    //         std::cout << '\n';
    //     }
    //     std::cout << '\n';
    // }
}


void create_support(Variables& var)
{
    var.support = new std::vector<int_vec>(var.nnode);

    // create the inverse mapping of connectivity
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            (*var.support)[conn[i]].push_back(e);
        }
    }
    // std::cout << "support:\n";
    // print(std::cout, *var.support);
    // std::cout << "\n";
}


void create_elemmarkers(const Param& param, Variables& var)
{
    var.elemmarkers = new int_vec2D( var.nelem, std::vector<int>(param.mat.nmat, 0) );
    if (param.control.has_hydration_processes)
        var.hydrous_elemmarkers = new Array2D<int,1>( var.nelem, 0 );

}


void create_markers(const Param& param, Variables& var)
{
    var.markersets.push_back(new MarkerSet(param, var, std::string("markerset")));
    if (param.control.has_hydration_processes) {
        var.hydrous_marker_index = var.markersets.size();
        var.markersets.push_back(new MarkerSet(std::string("hydrous-markerset")));
    }
}


void create_new_mesh(const Param& param, Variables& var)
{
    switch (param.mesh.meshing_option) {
    case 1:
        new_mesh_uniform_resolution(param, var);
        break;
    case 2:
        new_mesh_refined_zone(param, var);
        break;
    case 90:
    case 91:
        new_mesh_from_polyfile(param, var);
        break;
    case 95:
#ifdef USEEXODUS
        new_mesh_from_exofile(param, var);
#else
        std::cout << "Error: Install Exodus library and rebuild with 'useexo' turned on in Makefile." << std::endl;
        std::exit(1);
#endif
        break;
    default:
        std::cout << "Error: unknown meshing option: " << param.mesh.meshing_option << '\n';
        std::exit(1);
    }

    if (param.mesh.is_discarding_internal_segments)
        discard_internal_segments(var.nseg, *var.segment, *var.segflag);

    renumbering_mesh(param, *var.coord, *var.connectivity, *var.segment, var.regattr);

    // std::cout << "segment:\n";
    // print(std::cout, *var.segment);
    // std::cout << '\n';
    // std::cout << "segflag:\n";
    // print(std::cout, *var.segflag);
    // std::cout << '\n';
}


double** elem_center(const array_t &coord, const conn_t &connectivity)
{
    /* Returns the centroid of the elements.
     * Note: center[0] == tmp
     * The caller is responsible to delete [] center[0] and center!
     */
    int nelem = connectivity.size();
    double *tmp = new double[nelem*NDIMS];
    double **center = new double*[nelem];
    #pragma omp parallel for default(none)          \
        shared(nelem, tmp, coord, connectivity, center)
    for(int e=0; e<nelem; e++) {
        const int* conn = connectivity[e];
        center[e] = tmp + e*NDIMS;
        for(int d=0; d<NDIMS; d++) {
            double sum = 0;
            for(int k=0; k<NODES_PER_ELEM; k++) {
                sum += coord[conn[k]][d];
            }
            center[e][d] = sum / NODES_PER_ELEM;
        }
    }
    return center;
}


