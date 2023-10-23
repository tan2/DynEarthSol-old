#ifndef DYNEARTHSOL3D_MESH_HPP
#define DYNEARTHSOL3D_MESH_HPP

void points_to_new_mesh(const Mesh &mesh, int npoints, const double *points,
                        int n_init_segments, const int *init_segments, const int *init_segflags,
                        int n_regions, const double *regattr,
                        double max_elem_size, int vertex_per_polygon,
                        int &nnode, int &nelem, int &nseg,
                        double *&pcoord, int *&pconnectivity,
                        int *&psegment, int *&psegflag, double *&pregattr);
void points_to_new_surface(const Mesh &mesh, int npoints, const double *points,
                           int n_init_segments, const int *init_segments, const int *init_segflags,
                           int n_regions, const double *regattr,
                           double max_elem_size, int vertex_per_polygon,
                           int &nnode, int &nelem, int &nseg,
                           double *&pcoord, int *&pconnectivity,
                           int *&psegment, int *&psegflag, double *&pregattr);
void renumbering_mesh(const Param& param, array_t &coord, conn_t &connectivity,
                      segment_t &segment, regattr_t *regattr);
void create_boundary_flags2(uint_vec &bcflag, int nseg,
                            const int *psegment, const int *psegflag);
void create_boundary_flags(Variables& var);
void create_boundary_nodes(Variables& var);
void create_boundary_facets(Variables& var);
void create_support(Variables& var);
void create_elemmarkers(const Param&, Variables&);
void create_markers(const Param&, Variables&);
void create_new_mesh(const Param&, Variables&);
double** elem_center(const array_t &coord, const conn_t &connectivity);

#endif
