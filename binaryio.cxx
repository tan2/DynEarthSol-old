#include <cstdio>
#include <cstring>
#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

/*****************************************************************************
 * The format of the binary file:
 * 1  The first 'headerlen' bytes are ASCII text.
 *   1.1  The 1st line in the header is the revision string. Starting with
 *        "# DynEarthSol ndims=%1 revision=%2", with %1 equal to 2 or 3
 *        (indicating 2D or 3D simulation) and %2 an integer.
 *   1.2  The following lines are 'name', 'position' pairs, separated by a
 *        TAB character. This line tells the name of the data and the
 *        starting position (in bytes) of the data in this file.
 * 2  The rests are binary data.
 ****************************************************************************/

namespace {
    const std::size_t headerlen = 4096;
    const char revision_str[] = "# DynEarthSol ndims="
#ifdef THREED
        "3"
#else
        "2"
#endif
        " revision=3\n";
}


/* Not using C++ stream IO for bulk file io since it can be much slower than C stdio. */

BinaryOutput::BinaryOutput(const char *filename)
{
    f = std::fopen(filename, "wb");
    if (f == NULL) {
        std::cerr << "Error: cannot open file: " << filename << '\n';
        std::exit(2);
    }

    header = new char[headerlen]();
    hd_pos = std::strcat(header, revision_str);
    eof_pos = headerlen;

    std::fseek(f, eof_pos, SEEK_SET);
}


BinaryOutput::~BinaryOutput()
{
    close();
}


void BinaryOutput::close()
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    if (f) {
        /* write header buffer to the beginning of file */
        std::fseek(f, 0, SEEK_SET);
        std::fwrite(header, sizeof(char), headerlen, f);
        std::fclose(f);
        f = NULL;
    }

    delete [] header;
    header = NULL;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void BinaryOutput::write_header(const char *name)
{
    /* write to header buffer */
    const std::size_t bsize = 256;
    char buffer[bsize];
    std::size_t len = std::snprintf(buffer, bsize, "%s\t%ld\n", name, eof_pos);
    if (len >= bsize) {
        std::cerr << "Error: exceeding buffer length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        std::exit(12);
    }
    if (len >= headerlen - (hd_pos - header)*sizeof(char)) {
        std::cerr << "Error: exceeding header length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        std::exit(12);
    }
    hd_pos = std::strncat(hd_pos, buffer, len);
}

// XXX: when A is *var.bcflag, i.e. T is uint, g++ cannot instantiate the template
template <typename T>
void BinaryOutput::write_array(const std::vector<T>& A, const char *name, std::size_t size)
{
    write_header(name);
    std::size_t n = std::fwrite(A.data(), sizeof(T), size, f);
    eof_pos += n * sizeof(T);
}


// specialize for uint
void BinaryOutput::write_array(const std::vector<uint>& A, const char *name, std::size_t size)
{
    write_header(name);
    std::size_t n = std::fwrite(A.data(), sizeof(uint), size, f);
    eof_pos += n * sizeof(uint);
}


template <typename T, int N>
void BinaryOutput::write_array(const Array2D<T,N>& A, const char *name, std::size_t size)
{
    write_header(name);
    std::size_t n = std::fwrite(A.data(), sizeof(T), size*N, f);
    eof_pos += n * sizeof(T);
}


// explicit instantiation
template
void BinaryOutput::write_array<int>(const std::vector<int>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<double>(const std::vector<double>& A, const char *name, std::size_t);

template
void BinaryOutput::write_array<double,NDIMS>(const Array2D<double,NDIMS>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<double,NSTR>(const Array2D<double,NSTR>& A, const char *name, std::size_t);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void BinaryOutput::write_array<double,NODES_PER_ELEM>(const Array2D<double,NODES_PER_ELEM>& A, const char *name, std::size_t);
#endif
template
void BinaryOutput::write_array<double,1>(const Array2D<double,1>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<int,NODES_PER_ELEM>(const Array2D<int,NODES_PER_ELEM>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<int,NDIMS>(const Array2D<int,NDIMS>& A, const char *name, std::size_t);
template
void BinaryOutput::write_array<int,1>(const Array2D<int,1>& A, const char *name, std::size_t);

//////////////////////////////////////////////////////////////////////////////

BinaryInput::BinaryInput(const char *filename)
{
    f = std::fopen(filename, "r");
    if (f == NULL) {
        std::cerr << "Error: cannot open file: " << filename << '\n';
        std::exit(2);
    }
    read_header();
}


BinaryInput::~BinaryInput()
{
    std::fclose(f);
}


void BinaryInput::read_header()
{
    /* Read into header buffer */
    std::fseek(f, 0, SEEK_SET);
    char *header = new char[headerlen]();
    std::size_t n = std::fread(header, sizeof(char), headerlen, f);
    if (n != headerlen) {
        std::cerr << "Error: error reading file header\n";
        std::exit(2);
    }

    /* Parse the content of header buffer */
    char *line = header;

    // Compare revision string (excluding the trailing new line)
    line = std::strtok(header, "\n");
    if (strncmp(line, revision_str, strlen(revision_str)-1) != 0) {
        std::cerr << "Error: mismatching revision string in header\n"
                  << "  Expect: " << revision_str
                  << "  Got: "<< line << '\n';
        std::exit(1);
    }

    line = std::strtok(NULL, "\n");
    while (line != NULL) {
        /* Each line is a string (might contain space), a tab, and an integer */
        char *tab = std::strchr(line, '\t');
        if (tab == NULL) {
            std::cerr << "Error: error parsing file header\n"
                      << " Line is:" << line << '\n';
            std::exit(1);
        }
        std::string name(line, tab-line);
        std::size_t loc;
        std::sscanf(tab, "%zu", &loc);

        offset[name] = loc;
        line = std::strtok(NULL, "\n");
    }

    delete [] header;
}


void BinaryInput::seek_to_array(const char *name)
{
    std::string name2(name);
    auto it = offset.find(name);
    if (it == offset.end()) {
        std::cerr << "Error: no array with a name: " << name << '\n';
        std::exit(1);
    }
    std::size_t loc = it->second;
    //std::cout << name << ' ' << loc << '\n';
    std::fseek(f, loc, SEEK_SET);
}


template <typename T>
void BinaryInput::read_array(std::vector<T>& A, const char *name)
{
    /* The caller must ensure A is of right size to hold the array */

    int size = A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    seek_to_array(name);
    int n = std::fread(A.data(), sizeof(T), size, f);

    if (n != size) {
        std::cerr << "Error: cannot read array: " << name << '\n';
        std::exit(1);
    }
}


template <typename T, int N>
void BinaryInput::read_array(Array2D<T,N>& A, const char *name)
{
    /* The caller must ensure A is of right size to hold the array */

    int size = A.size();
    if (A.size() == 0) {
        std::cerr << "Error: array size is 0: " << name << '\n';
        std::exit(1);
    }

    seek_to_array(name);
    int n = std::fread(A.data(), sizeof(T), size*N, f);

    if (n != N*size) {
        std::cerr << "Error: cannot read array: " << name << '\n';
        std::exit(1);
    }
}


// explicit instantiation
template
void BinaryInput::read_array<double>(std::vector<double>& A, const char *name);
template
void BinaryInput::read_array<int>(std::vector<int>& A, const char *name);
template
void BinaryInput::read_array<double,NDIMS>(Array2D<double,NDIMS>& A, const char *name);
template
void BinaryInput::read_array<double,NSTR>(Array2D<double,NSTR>& A, const char *name);
#ifdef THREED // when 2d, NSTR == NODES_PER_ELEM == 3
template
void BinaryInput::read_array<double,NODES_PER_ELEM>(Array2D<double,NODES_PER_ELEM>& A, const char *name);
#endif
template
void BinaryInput::read_array<double,1>(Array2D<double,1>& A, const char *name);
template
void BinaryInput::read_array<int,NDIMS>(Array2D<int,NDIMS>& A, const char *name);
template
void BinaryInput::read_array<int,NODES_PER_ELEM>(Array2D<int,NODES_PER_ELEM>& A, const char *name);
template
void BinaryInput::read_array<int,1>(Array2D<int,1>& A, const char *name);

