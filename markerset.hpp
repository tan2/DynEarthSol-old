#ifndef DYNEARTHSOL3D_MARKERSET_HPP
#define DYNEARTHSOL3D_MARKERSET_HPP

#include <string>
#include "barycentric-fn.hpp"


// forward declaration
class BinaryOutput;
class BinaryInput;

class MarkerSet
{

public:
    MarkerSet( const std::string& name );
    MarkerSet( const Param& param, Variables& var, const std::string& name );
    MarkerSet( const Param& param, Variables& var, BinaryInput& bin, const std::string& name );
    ~MarkerSet()
    {
        delete _slope;
        delete _distance;
        delete _z;
        delete _time;
        delete _id;
        delete _eta; 
        delete _elem;
        delete _mattype;
    }

    static void random_eta( double* ); // class method
//    void create_marker_in_elem(Variables& var);
//    void update_marker_in_elem(Variables& var);
    void create_melt_markers(const int mat, int_vec& melt_markers);
    void correct_surface_marker(const Variables& var, array_t& coord0s, Barycentric_transformation &bary);
    void set_surface_marker(const Variables& var, const double smallest_size, const int mattype_sed, array_t& edhacc, int_vec2D& elemmarkers);
    void remap_marker(const Variables &var, const double *m_coord, const int e, int &new_elem, double *new_eta, int &inc);
    void append_random_marker_in_elem( int el, int mt);
    void append_marker( const double *eta, int el, int mt, double time, double z, double distance, double slope);
    void remove_marker(int i);
    void resize(const int);
    void write_chkpt_file(BinaryOutput &bin) const;
    void read_chkpt_file(Variables &var, BinaryInput &bin);
    void write_save_file(const Variables &var, BinaryOutput &bin) const;

    inline int get_nmarkers() const { return _nmarkers; }
    inline void set_nmarkers(int n) { _nmarkers = n; }

    inline bool if_melt(const int mat) const { return (std::find((*_mattype).begin(), (*_mattype).end(), mat) != (*_mattype).end()); }

    inline int get_id(int m) const { return (*_id)[m]; }
    inline void set_id(const int m, const int i) { (*_id)[m] = i; }

    inline int get_elem(int m) const { return (*_elem)[m]; }
    inline void set_elem(const int m, const int e) { (*_elem)[m] = e; }

    inline int get_mattype(int m) const { return (*_mattype)[m]; }
    inline void set_mattype(const int m, const int mt) { (*_mattype)[m] = mt; }

    inline double get_time(int m) const { return (*_time)[m]; }
    inline void set_time(const int m, const double ti) { (*_time)[m] = ti; }

    inline double get_z(int m) const { return (*_z)[m]; }
    inline void set_z(const int m, const double z) { (*_z)[m] = z; }

    inline double get_distance(int m) const { return (*_distance)[m]; }
    inline void set_distance(const int m, const double d) { (*_distance)[m] = d; }

    inline double get_slope(int m) const { return (*_slope)[m]; }
    inline void set_slope(const int m, const double s) { (*_slope)[m] = s; }

    inline const double *get_eta(int m) const { return (*_eta)[m]; }
    inline void set_eta( const int i, const double r[NDIMS] );

private:
    const std::string _name;

    // Didn't create a data type for an individual marker to follow the "structure of arrays" concept.
    // Number of markers (may change during remeshing)

    int _nmarkers;
    int _reserved_space;
    int _last_id;

    // Barycentric (local) coordinate within the reference element
    shapefn *_eta;
    // Containing element
    int_vec *_elem;
    // Material type
    int_vec *_mattype;
    // Unique id
    int_vec *_id;
    // Cearte time
    double_vec *_time;
    // Cearte z
    double_vec *_z;
    // Distance to coastline
    double_vec *_distance;
    // Slope of surface
    double_vec *_slope;

    void random_markers( const Param&, Variables& );
    void regularly_spaced_markers( const Param&, Variables& );
    void allocate_markerdata( const int );

    int initial_mattype( const Param&, const Variables&,
                         int elem, const double eta[NODES_PER_ELEM],
                         const double *x=NULL );
    int layered_initial_mattype( const Param& param, const Variables &var,
                                 int elem, const double eta[NODES_PER_ELEM],
                                 const double *x);
    int custom_initial_mattype( const Param& param, const Variables &var,
                                int elem, const double eta[NODES_PER_ELEM],
                                const double *x );

};

void remap_markers(const Param&, Variables &, 
                   const array_t &, const conn_t &);
void advect_hydrous_markers(const Param &, const Variables &, double,
                            MarkerSet &, Array2D<int,1> &);

#endif
