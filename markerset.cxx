#ifdef USE_NPROF
#include <nvToolsExt.h>
#endif
#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>
#include <unordered_map>

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
#include "binaryio.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "geometry.hpp"
#include "utils.hpp"


namespace {

    const int DEBUG = 0;
    const double over_alloc_ratio = 2.0;  // how many extra space to allocate for future expansion

}


MarkerSet::MarkerSet(const std::string& name) :
    _name(name)
{
    _last_id = _nmarkers = 0;
    allocate_markerdata(4 * 1024); // pre-allocate a small amount of markers
}


MarkerSet::MarkerSet(const Param& param, Variables& var, const std::string& name) :
    _name(name)
{
    _last_id = _nmarkers = 0;

    switch ( param.markers.init_marker_option ) {
    case 1:
        random_markers(param, var);
        break;
    case 2:
        regularly_spaced_markers(param, var);
        break;
    default:
        std::cerr << "Error: unknown init_marker_option: " << param.markers.init_marker_option << ". The only valid option is '1'.\n";
        std::exit(1);
        break;
    }

    for( int e = 0; e < var.nelem; e++ ) {
        int num_markers_in_elem = 0;
        for( int i = 0; i < param.mat.nmat; i++ )
            num_markers_in_elem += (*(var.elemmarkers))[e][i];

        if (num_markers_in_elem <= 0) {
            std::cerr << "Error: no marker in element #" << e
                      << ". Please increase the number of markers.\n";
            std::exit(1);
        }
    }
}


MarkerSet::MarkerSet(const Param& param, Variables& var, BinaryInput& bin, const std::string& name) :
    _name(name)
{
    // init from checkpoint file
    read_chkpt_file(var, bin);
}


void MarkerSet::allocate_markerdata( const int max_markers )
{
    _reserved_space = max_markers;
    _eta = new shapefn( max_markers );
    _elem = new int_vec( max_markers );
    _mattype = new int_vec( max_markers );
    _id = new int_vec( max_markers );
}


void MarkerSet::random_eta( double *eta )
{
    // eta for randomly scattered markers within an element.
    // An alternative would be to fix barycentric coordinates and add random perturbations.
    //
    while (1) {
        // sum(eta) == 1 && all components of eta are greater than zero
        double sum = 0;
        for( int n = 0; n < NDIMS; n++ ) {
            eta[n] = (rand()/(double)RAND_MAX);
            sum += eta[n];
        }
        if (sum < 1) {
            eta[NODES_PER_ELEM - 1] = 1 - sum;
            break;
        }
    }
}


void MarkerSet::append_marker( const double *eta, int el, int mt )
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // Ensure sufficient array size
    if( _nmarkers == _reserved_space ) {
        // Resize the marker-related arrays if necessary.
        const int newsize = _nmarkers * over_alloc_ratio;
        resize( newsize );
    }

    int m = _nmarkers;
    std::memcpy((*_eta)[m], eta, NODES_PER_ELEM*sizeof(double));
    (*_elem)[m] = el;
    (*_mattype)[m] = mt;
    (*_id)[m] = _last_id;

    if(DEBUG > 1) {
        std::cout << el << " " << m << " "
                  << _nmarkers << " " << (*_mattype)[_nmarkers] << " "
                  << eta[0] << "+" << eta[1] << "+" << eta[2];
#ifdef THREED
        std::cout << "+" << eta[3]
                  << "=" << (eta[0]+eta[1]+eta[2]+eta[3]) << "\n";
#else
        std::cout << "=" << (eta[0]+eta[1]+eta[2]) << "\n";
#endif
    }

    ++_nmarkers;
    ++_last_id;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::append_random_marker_in_elem( int el, int mt )
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::random_markers( const Param& param, Variables &var )
{
    const int ne = var.nelem;
    const int mpe = param.markers.markers_per_element;
    const int num_markers = ne*mpe;
    const int max_markers = num_markers * over_alloc_ratio;

    // allocate memory for data members.
    allocate_markerdata( max_markers );
    
    // initialize random seed:
    if (param.markers.random_seed)
        srand(param.markers.random_seed);
    else
        srand(time(NULL));

    // Generate particles in each element.
    for( int e = 0; e < ne; e++ )
        for( int m = 0; m < mpe; m++ ) {
            // random barycentric coordinate
            double eta[NODES_PER_ELEM];
            random_eta(eta);

            // decide the mattype of markers
            int mt = initial_mattype(param, var, e, eta);
            append_marker(eta, e, mt);
            ++(*var.elemmarkers)[e][mt];
        }
}


void MarkerSet::regularly_spaced_markers( const Param& param, Variables &var )
{
    const int d = param.markers.init_marker_spacing * param.mesh.resolution;

    double domain_min[NDIMS], domain_max[NDIMS];
    {
        for (int d=0; d<NDIMS; d++)
            domain_min[d] = domain_max[d] = (*var.coord)[0][d];

        for (int i=1; i<var.nnode; i++) {
            for (int d=0; d<NDIMS; d++) {
                domain_min[d] = std::min(domain_min[d], (*var.coord)[i][d]);
                domain_max[d] = std::max(domain_max[d], (*var.coord)[i][d]);
            }
        }
        // print(std::cout, domain_min, NDIMS);
        // print(std::cout, domain_max, NDIMS);
    }

    const double xlength = domain_max[0] - domain_min[0];
    const int nx = xlength / d + 1;
    const double x0 = domain_min[0] + 0.5 * (xlength - (nx-1)*d);
    const double zlength = domain_max[NDIMS-1] - domain_min[NDIMS-1];
    const int nz = zlength / d + 1;
    const double z0 = domain_min[NDIMS-1] + 0.5 * (zlength - (nz-1)*d);
#ifdef THREED
    const double ylength = domain_max[1] - domain_min[1];
    const int ny = ylength / d + 1;
    const double y0 = domain_min[1] + 0.5 * (ylength - (ny-1)*d);
#else
    const int ny = 1;
#endif

    const int num_markers = nx * ny * nz;
    const int max_markers = num_markers * over_alloc_ratio;

    allocate_markerdata( max_markers );

    // nearest-neighbor search structure
    double **centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements
    ANNkd_tree kdtree(centroid, var.nelem, NDIMS);
    const int k = std::min(20, var.nelem);  // how many nearest neighbors to search?
    const double eps = 0.001; // tolerance of distance error
    int nn_idx[20];
    double dd[20];

    double_vec new_volume( var.nelem );
    compute_volume( *var.coord, *var.connectivity, new_volume );
    Barycentric_transformation bary(*var.coord, *var.connectivity, new_volume);

    for (int n=0; n< num_markers; ++n) {
        int ix = n % nx;
        int iy = (n / nx) % ny;
        int iz = n / (nx * ny);

        // Physical coordinate of new marker
        double x[NDIMS] = { x0 + ix*d,
#ifdef THREED
                            y0 + iy*d,
#endif
                            z0 + iz*d };
        bool found = false;

        // Look for nearby elements.
        // Note: kdtree.annkSearch() is not thread-safe, cannot use openmp in this loop
        kdtree.annkSearch(x, k, nn_idx, dd, eps);

        for( int j=0; j<k; j++ ) {
            int e = nn_idx[j];
            double eta[NODES_PER_ELEM];

            bary.transform(x, e, eta);

            // Compute the last component of eta, with the constraint sum(eta)==1
            double tmp = 1;
            for( int d=0; d<NDIMS; ++d) {
                tmp -= eta[d];
            }
            eta[NDIMS] = tmp;

            if (bary.is_inside(eta)) {
                int mt = initial_mattype(param, var, e, eta, x);
                append_marker(eta, e, mt);
                ++(*var.elemmarkers)[e][mt];
                found = true;
                break;
            }
        }

        if (! found) {
            // x is outside the domain (ex: the domain is not rectangular)
            continue;
        }
    }

    delete [] centroid[0];
    delete [] centroid;
}


int MarkerSet::initial_mattype( const Param& param, const Variables &var,
                                int elem, const double eta[NODES_PER_ELEM],
                                const double *x)
{
    int mt;
    if (param.ic.mattype_option == 0) {
        mt = (*var.regattr)[elem][0]; // mattype should take a reginal attribute assigned during meshing.
    }
    else {
        double p[NDIMS] = {0};
        if (x) {
            std::memcpy(p, x, NDIMS*sizeof(double));
        }
        else {
            const int *conn = (*var.connectivity)[elem];
            for(int i=0; i<NDIMS; i++) {
                for(int j=0; j<NODES_PER_ELEM; j++)
                    p[i] += (*var.coord)[ conn[j] ][i] * eta[j];
            }
        }
        // modify mt according to the marker coordinate p
        switch (param.ic.mattype_option) {
        case 1:
            mt = layered_initial_mattype(param, var, elem, eta, p);
            break;
        case 101:
            mt = custom_initial_mattype(param, var, elem, eta, p);
            break;
        default:
            std::cerr << "Error: unknown ic.mattype_option: " << param.ic.mattype_option << '\n';
            std::exit(1);
        }
    }
    return mt;
}


int MarkerSet::layered_initial_mattype( const Param& param, const Variables &var,
                                        int elem, const double eta[NODES_PER_ELEM],
                                        const double *x)
{
    int mt = param.ic.layer_mattypes[param.ic.layer_mattypes.size() - 1];
    const double_vec &layers = param.ic.mattype_layer_depths;
    for (std::size_t i=0; i<layers.size(); ++i) {
        if (x[NDIMS-1] >= -param.mesh.zlength * layers[i]) {
            mt = param.ic.layer_mattypes[i];
            break;
        }
    }
    return mt;
}


int MarkerSet::custom_initial_mattype( const Param& param, const Variables &var,
                                       int elem, const double eta[NODES_PER_ELEM],
                                       const double *x )
{
    /* User defined function */
    int mt = 0;

    return mt;
}


void MarkerSet::remove_marker(int i)
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // Replace marker i by the last marker.
    --_nmarkers;
    std::memcpy( (*_eta)[i], (*_eta)[_nmarkers], sizeof(double)*(NODES_PER_ELEM) );
    (*_id)[i] = (*_id)[_nmarkers];
    (*_elem)[i] = (*_elem)[_nmarkers];
    (*_mattype)[i] = (*_mattype)[_nmarkers];
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::set_eta( const int i, const double r[NDIMS] ) {
    double sum = 0.0;
    for( int j = 0; j < NDIMS; j++ ) {
        (*_eta)[i][j] = r[j];
        sum += r[j];
    }
    (*_eta)[i][NDIMS] = 1.0 - sum;
}


void MarkerSet::resize( const int newsize )
{
    if( newsize > _reserved_space ) {
        // enlarge arrays
        std::cout << "  Increasing marker arrays size to " << newsize << " markers.\n";
        _reserved_space = newsize;

        double *new_eta = new double[NODES_PER_ELEM * newsize];
        std::copy( (*_eta)[0], (*_eta)[_nmarkers], new_eta );
        _eta->reset( new_eta, newsize );

        _elem->resize( newsize );
        _mattype->resize( newsize );
        _id->resize( newsize );
    }
    // else if( nmarkers_new < _reserved_space ) {
    //     // TBD: shrink arrays
    //     _reserved_space = newsize;
    // }
    else {
        // New size is too close to old size, don't do anything.
    }
}


void MarkerSet::write_chkpt_file(BinaryOutput &bin) const
{
    int_vec itmp(2);
    itmp[0] = _nmarkers;
    itmp[1] = _last_id;
    bin.write_array(itmp, (_name + " size").c_str(), itmp.size());

    bin.write_array(*_eta, (_name + ".eta").c_str(), _nmarkers);
    bin.write_array(*_elem, (_name + ".elem").c_str(), _nmarkers);
    bin.write_array(*_mattype, (_name + ".mattype").c_str(), _nmarkers);
    bin.write_array(*_id, (_name + ".id").c_str(), _nmarkers);

}


void MarkerSet::read_chkpt_file(Variables &var, BinaryInput &bin)
{
    int_vec itmp(2);
    bin.read_array(itmp, (_name + " size").c_str());
    _nmarkers = itmp[0];
    _last_id = itmp[1];

    allocate_markerdata(_nmarkers);

    if (_nmarkers != 0) {
        bin.read_array(*_eta, (_name + ".eta").c_str());
        bin.read_array(*_elem, (_name + ".elem").c_str());
        bin.read_array(*_mattype, (_name + ".mattype").c_str());
        bin.read_array(*_id, (_name + ".id").c_str());
    }

    if (_name == "markerset")
        for( int i = 0; i < _nmarkers; i++ ) {
            int e = (*_elem)[i];
            int mt = (*_mattype)[i];
            ++(*var.elemmarkers)[e][mt];
        }
    else if (_name == "hydrous-markerset")
        for( int i = 0; i < _nmarkers; i++ ) {
            int e = (*_elem)[i];
            ++(*var.hydrous_elemmarkers)[e][0];
        }
}


void MarkerSet::write_save_file(const Variables &var, BinaryOutput &bin) const
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    int_vec itmp(1);
    itmp[0] = _nmarkers;
    bin.write_array(itmp, (_name + " size").c_str(), itmp.size());

    array_t mcoord(_nmarkers, 0);
    const array_t &coord = *var.coord;

    for (int i=0; i<_nmarkers; ++i) {
        int e = (*_elem)[i];
        const int *conn = (*var.connectivity)[e];
        const double *eta = (*_eta)[i];
        double *x = mcoord[i];
        for (int j = 0; j < NDIMS; j++)
            for (int k = 0; k < NODES_PER_ELEM; k++)
                x[j] += eta[k] * coord[ conn[k] ][j];

        // std::cout << i << '\t';
        // print(std::cout, x, NDIMS);
        // std::cout << "\n";
    }

    bin.write_array(mcoord, (_name + ".coord").c_str(), _nmarkers);
    bin.write_array(*_elem, (_name + ".elem").c_str(), _nmarkers);
    bin.write_array(*_mattype, (_name + ".mattype").c_str(), _nmarkers);
    bin.write_array(*_id, (_name + ".id").c_str(), _nmarkers);

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

    template <class T>
    void find_markers_in_element(MarkerSet& ms, T& elemmarkers,
                                 ANNkd_tree& kdtree, const Barycentric_transformation &bary,
                                 const array_t& old_coord, const conn_t &old_connectivity)
    {
        const int k = std::min((std::size_t) 20, old_connectivity.size());  // how many nearest neighbors to search?
        const double eps = 0.001; // tolerance of distance error
        int nn_idx[20];
        double dd[20];


        // Loop over all the old markers and identify a containing element in the new mesh.
        int last_marker = ms.get_nmarkers();
        int i = 0;
        while (i < last_marker) {
            bool found = false;

            // 1. Get physical coordinates, x, of an old marker.
            int eold = ms.get_elem(i);
            double x[NDIMS] = {0};
            for (int j = 0; j < NDIMS; j++)
                for (int k = 0; k < NODES_PER_ELEM; k++)
                    x[j] += ms.get_eta(i)[k]*
                        old_coord[ old_connectivity[eold][k] ][j];

            if (DEBUG) {
                std::cout << "marker #" << i << " old_elem " << eold << " x: ";
                print(std::cout, x, NDIMS);
            }

            // 2. look for nearby elements.
            // Note: kdtree.annkSearch() is not thread-safe, cannot use openmp in this loop
            kdtree.annkSearch(x, k, nn_idx, dd, eps);

            for( int j = 0; j < k; j++ ) {
                int e = nn_idx[j];
                double r[NDIMS];

                bary.transform(x, e, r);

                // change this if-condition to (i == N) to debug the N-th marker
                if (0) {
                    std::cout << '\n' << j << " check elem #" << e << ' ';
                    print(std::cout, r, NDIMS);
                }

                if (bary.is_inside(r)) {
                    ms.set_eta(i, r);
                    ms.set_elem(i, e);
                    ++elemmarkers[e][ms.get_mattype(i)];
            
                    found = true;
                    ++i;
                    if (DEBUG) {
                        std::cout << " in element " << e << '\n';
                    }
                    break;
                }
            }

            if( found ) continue;

            if (DEBUG) {
                std::cout << " not in any element" << '\n';
            }

            /* not found */
            {
                // Since no containing element has been found, delete this marker.
                // Note i is not inc'd.
                --last_marker;
                ms.remove_marker(i);
            }
        }
    }


    void replenish_markers_with_mattype_0(const Param& param, Variables &var,
                                          int e, int num_marker_in_elem)
    {
#ifdef USE_NPROF
        nvtxRangePush(__FUNCTION__);
#endif
        while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
            const int mt = 0;
            var.markersets[0]->append_random_marker_in_elem(e, mt);
            if (DEBUG) {
                std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
            }

            ++(*var.elemmarkers)[e][mt];
            ++num_marker_in_elem;
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_cpdf(const Param& param, Variables &var,
                                                  int e, int num_marker_in_elem)
    {
#ifdef USE_NPROF
        nvtxRangePush(__FUNCTION__);
#endif
        // cummulative probability density function of mattype
        double_vec cpdf(param.mat.nmat, 0);

        if (num_marker_in_elem > 0) {
            // cpdf of this element
            cpdf[0] = (*(var.elemmarkers))[e][0] / double(num_marker_in_elem);
            for( int i = 1; i < param.mat.nmat - 1; i++ )
                cpdf[i] = cpdf[i-1] + (*(var.elemmarkers))[e][i] / double(num_marker_in_elem);
        }
        else {
            // No markers in this element.
            // Construct cpdf from neighboring elements
            int num_markers_in_nbr_elems = 0;
            // Looping over all neighboring elements (excluding self)
            for( int kk = 0; kk < NODES_PER_ELEM; kk++) {
                int n = (*var.connectivity)[e][kk]; // node of this element
                for( auto ee = (*var.support)[n].begin(); ee < (*var.support)[n].end(); ++ee) {
                    if (*ee == e) continue;
                    // Note: some (NODES_PER_ELEM) elements will be iterated
                    // more than once (NDIMS times). These elements are direct neighbors,
                    // i.e. they share facets (3D) or edges (2D) with element e.
                    // So they are counted multiple times to reprensent a greater weight.
                    for( int i = 0; i < param.mat.nmat; i++ ) {
                        cpdf[i] += (*(var.elemmarkers))[*ee][i];
                        num_markers_in_nbr_elems += (*(var.elemmarkers))[*ee][i];
                    }
                }
            }
            for( int i = 1; i < param.mat.nmat - 1; i++ )
                cpdf[i] += cpdf[i-1];
            for( int i = 0; i < param.mat.nmat - 1; i++ )
                cpdf[i] = cpdf[i] / double(num_markers_in_nbr_elems);
        }

        cpdf[param.mat.nmat - 1] = 1; // fix to 1 to avoid round-off error
        if (DEBUG > 1) {
            std::cout << num_marker_in_elem << " markers in element " << e << '\n'
                      << "  cpdf: ";
            print(std::cout, cpdf);
            std::cout << '\n';
        }

        while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
            // Determine new marker's matttype based on cpdf
            auto upper = std::upper_bound(cpdf.begin(), cpdf.end(), rand()/(double)RAND_MAX);
            const int mt = upper - cpdf.begin();
            var.markersets[0]->append_random_marker_in_elem(e, mt);
            if (DEBUG) {
                std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
            }

            ++(*var.elemmarkers)[e][mt];
            ++num_marker_in_elem;
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_nn_preparation(const Param& param, const Variables &var,
                                                            ANNkd_tree *&kdtree, double **&points)
    {
#ifdef USE_NPROF
        nvtxRangePush(__FUNCTION__);
#endif
        const MarkerSet &ms = *var.markersets[0];
        const int nmarkers = ms.get_nmarkers();

        double *tmp = new double[nmarkers*NDIMS];
        points = new double*[nmarkers];
        // The caller is responsible to delete [] points[0] and points!

        for (int n=0; n<nmarkers; n++) {
            const int e = ms.get_elem(n);
            const double* eta = ms.get_eta(n);
            const int* conn = (*var.connectivity)[e];

            points[n] = tmp + n*NDIMS;
            for(int d=0; d<NDIMS; d++) {
                double sum = 0;
                for(int k=0; k<NODES_PER_ELEM; k++) {
                    sum += (*var.coord)[ conn[k] ][d] * eta[k];
                }
                points[n][d] = sum;
            }
        }

        kdtree = new ANNkd_tree(points, nmarkers, NDIMS);
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_nn(const Param& param, Variables &var,
                                                ANNkd_tree &kdtree,
                                                int e, int num_marker_in_elem)
    {
#ifdef USE_NPROF
        nvtxRangePush(__FUNCTION__);
#endif
        MarkerSet &ms = *var.markersets[0];
        while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
            double eta[NODES_PER_ELEM];
            ms.random_eta(eta);

            double x[NDIMS] = {0};
            const int *conn = (*var.connectivity)[e];
            for (int d=0; d<NDIMS; d++) {
                for (int i=0; i<NODES_PER_ELEM; i++) {
                    x[d] += (*var.coord)[ conn[i] ][d] * eta[i];
                }
            }

            const int k = 1;  // how many nearest neighbors to search?
            const double eps = 0; // tolerance of distance error
            int nn_idx[k];
            double dd[k];

            // Look for nearest marker.
            // Note: kdtree.annkSearch() is not thread-safe, cannot use openmp in this loop
            kdtree.annkSearch(x, k, nn_idx, dd, eps);

            int m = nn_idx[0]; // nearest marker
            const int mt = ms.get_mattype(m);

            ms.append_marker(eta, e, mt);
            if (DEBUG) {
                std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
            }

            ++(*var.elemmarkers)[e][mt];
            ++num_marker_in_elem;
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

} // anonymous namespace


void remap_markers(const Param& param, Variables &var, const array_t &old_coord,
                   const conn_t &old_connectivity)
{
#ifdef USE_NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // Re-create elemmarkers
    delete var.elemmarkers;
    if (param.control.has_hydration_processes)
        delete var.hydrous_elemmarkers;
    create_elemmarkers( param, var );

    // Locating markers in new elements
    {
        double_vec new_volume( var.nelem );
        compute_volume( *var.coord, *var.connectivity, new_volume );

        Barycentric_transformation bary( *var.coord, *var.connectivity, new_volume );

        // nearest-neighbor search structure
        double **centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements
        ANNkd_tree kdtree(centroid, var.nelem, NDIMS);

        find_markers_in_element(*var.markersets[0], *var.elemmarkers,
                                kdtree, bary, old_coord, old_connectivity);
        if (param.control.has_hydration_processes)
            find_markers_in_element(*var.markersets[var.hydrous_marker_index], *var.hydrous_elemmarkers,
                                    kdtree, bary, old_coord, old_connectivity);

        delete [] centroid[0];
        delete [] centroid;
    }

    // If any new element has too few markers, generate markers in them.
    ANNkd_tree *kdtree = NULL;
    double **points = NULL; // coordinate of markers
    if (param.markers.replenishment_option == 2) {
        replenish_markers_with_mattype_from_nn_preparation(param, var, kdtree, points);
    }

    for( int e = 0; e < var.nelem; e++ ) {
        int num_marker_in_elem = 0;
        for( int i = 0; i < param.mat.nmat; i++ )
            num_marker_in_elem += (*(var.elemmarkers))[e][i];

        if (num_marker_in_elem < param.markers.min_num_markers_in_element) {
            switch (param.markers.replenishment_option) {
            case 0:
                replenish_markers_with_mattype_0(param, var, e, num_marker_in_elem);
                break;
            case 1:
                replenish_markers_with_mattype_from_cpdf(param, var, e, num_marker_in_elem);
                break;
            case 2:
                replenish_markers_with_mattype_from_nn(param, var, *kdtree, e, num_marker_in_elem);
                break;
            default:
                std::cerr << "Error: unknown markers.replenishment_option: " << param.markers.replenishment_option << '\n';
                std::exit(1);
            }
        }
    }

    if (param.markers.replenishment_option == 2) {
        delete kdtree;
        delete [] points[0];
        delete [] points;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

    Barycentric_transformation* get_bary_from_cache(std::unordered_map<int, Barycentric_transformation*> &cache,
                                                    int el, const array_t &coordinate, const int *conn,
                                                    double_vec &volume)
    {
        Barycentric_transformation *bary;
        auto search = cache.find(el);
        if (search == cache.end()) {
            const double *coord[NODES_PER_ELEM];
            for(int j=0; j<NODES_PER_ELEM; j++) {
                coord[j] = coordinate[ conn[j] ];
            }
            bary = new Barycentric_transformation(coord, volume[el]);
            cache[el] =  bary;
        }
        else {
            bary = search->second;
        }
        return bary;
    }
}


void advect_hydrous_markers(const Param& param, const Variables& var, double dt_subtotal,
                            MarkerSet& hydms, Array2D<int,1>& hydem)
{
    std::unordered_map<int, Barycentric_transformation*> cache;
    Barycentric_transformation *bary;

    int last_marker = hydms.get_nmarkers();
    int m = 0;
    while (m < last_marker) {
        // Find physical coordinate of the marker
        int el = hydms.get_elem(m);
        double x[NDIMS] = {0};
        const int *conn = (*var.connectivity)[el];
        for(int i=0; i<NDIMS; i++) {
            for(int j=0; j<NODES_PER_ELEM; j++)
                x[i] += hydms.get_eta(m)[j] * (*var.coord)[ conn[j] ][i];
        }

        // Advect the marker
        x[NDIMS-1] += dt_subtotal * param.control.hydration_migration_speed;

        // Transform back to barycentric coordinate
        double r[NDIMS];
        bary = get_bary_from_cache(cache, el, *var.coord, conn, *var.volume);
        bary->transform(x, 0, r); // always (local) 0-th element for bary

        if (bary->is_inside(r)) {
            hydms.set_eta(m, r);
            ++m;
            goto next;
        }
        else {
            // Marker has moved out of el. Find the new containing element.
            for(int j=0; j<NODES_PER_ELEM; j++) {
                const int_vec& supp = (*var.support)[ conn[j] ];
                for (std::size_t k=0; k<supp.size(); k++) {
                    int ee = supp[k];
                    const int *conn2 = (*var.connectivity)[ee];
                    bary = get_bary_from_cache(cache, ee, *var.coord, conn2, *var.volume);
                    bary->transform(x, 0, r);
                    if (bary->is_inside(r)) {
                        hydms.set_elem(m, ee);
                        hydms.set_eta(m, r);
                        --hydem[el][0];
                        ++hydem[ee][0];
                        ++m;
                        goto next;
                    }
                }
            }
            // Since no containing element has been found, delete this marker.
            // Note m is not inc'd.
            if (DEBUG) {
                std::cout << m << "-th hydrous marker not in any element\n";
            }
            --last_marker;
            hydms.remove_marker(m);
            --hydem[el][0];
        }
    next:;
    }

    // clean up all Barycentric_transformation instances
    for (auto i=cache.begin(); i!=cache.end(); ++i) {
        delete i->second;
    }
}


