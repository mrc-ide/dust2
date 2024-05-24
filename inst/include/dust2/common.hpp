#pragma once

// This will likely move around a bit.  Also note that we have an odd
// choice of directory here.
#include <mcstate/random/random.hpp>
#include <mcstate/random/density.hpp>
#include <dust2/tools.hpp>
#include <cpp11/list.hpp>

namespace dust2 {

struct no_data {};
struct no_internal_state {};
struct no_shared_state {};

namespace r {

// The actual definitions of these are elsewhere, but these are
// functions that models may use so we declare them here.
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

}

}
