#pragma once

// This will likely move around a bit.  Also note that we have an odd
// choice of directory here.
#include <monty/random/random.hpp>
#include <monty/random/density.hpp>
#include <dust2/array.hpp>
#include <dust2/packing.hpp>
#include <dust2/tools.hpp>
#include <dust2/zero.hpp>
#include <cpp11/list.hpp>

// In an odd place, we might update that later too
#include <dust2/interpolate/interpolate.hpp>

namespace dust2 {

struct no_data {};
struct no_internal_state {};

namespace r {

// The actual definitions of these are elsewhere, but these are
// functions that systems may use so we declare them here.
inline double read_real(cpp11::list pars, const char * name);
inline double read_real(cpp11::list pars, const char * name,
                        double default_value);

inline int read_int(cpp11::list pars, const char * name);
inline int read_int(cpp11::list pars, const char * name,
                    int default_value);

inline size_t read_size(cpp11::list pars, const char * name);
inline size_t read_size(cpp11::list pars, const char * name,
                        size_t default_value);

inline bool read_bool(cpp11::list pars, const char * name);
inline bool read_bool(cpp11::list pars, const char * name,
                      bool default_value);

// Read a vector into allocated storage, if it is found
template <typename real_type>
void read_real_vector(cpp11::list pars, size_t len, real_type * dest,
                             const char *name, bool required);
template <typename real_type, size_t rank>
void read_real_array(cpp11::list args, const dust2::array::dimensions<rank>& dim,
                     real_type * dest, const char *name, bool required);
template <size_t rank>
void read_int_array(cpp11::list args, const dust2::array::dimensions<rank>& dim,
                    int * dest, const char *name, bool required);
template <size_t rank>
dust2::array::dimensions<rank> read_dimensions(cpp11::list args, const char * name);
template <size_t rank>
void read_dimensions(cpp11::list args, const char * name, dust2::array::dimensions<rank>& dest);
template <>
void read_dimensions(cpp11::list args, const char * name, dust2::array::dimensions<1>& dest);

}

}
