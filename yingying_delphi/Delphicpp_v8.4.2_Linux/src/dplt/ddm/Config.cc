
#include "../ddm/internal/Config.h"
#include "../ddm/internal/Logging.h"

#include "../ddm/util/Config.h"
#include "../ddm/util/Config.h"

#include <unordered_map>
#include <string>
#include <algorithm>
#include <cstring>
#include <unistd.h>
#include <limits>


// Environment variables as array of strings, terminated by null pointer.
extern char ** environ;

namespace ddm {
namespace util {

std::unordered_map<std::string, Config::callback_fun>  Config::callbacks_;
std::unordered_map<std::string, std::string>           Config::config_values_;

void Config::init()
{
  int    i          = 1;
  char * env_var_kv = *environ;

  Config::callbacks_["DDM_ENABLE_LOGGING"] =
    &Config::ddm_enable_logging_callback;
  Config::callbacks_["DDM_ENABLE_LOGGING_BOOL"] =
    &Config::ddm_enable_logging_callback;

  for (; env_var_kv != 0; ++i) {
    // Split into key and value:
    char * flag_name_cstr  = env_var_kv;
    char * flag_value_cstr = strstr(env_var_kv, "=");
    int    flag_name_len   = flag_value_cstr - flag_name_cstr;
    std::string flag_name(flag_name_cstr, flag_name_cstr + flag_name_len);
    std::string flag_value(flag_value_cstr+1);
    // Flag name has prefix "DDM_":
    if (strstr(env_var_kv, "DDM_") == env_var_kv) {
      set(flag_name, flag_value);
    }
    else if (strstr(env_var_kv, "MP_RDMA_MTU") == env_var_kv) {
      set("DDM_RDMA_MTU_SIZE", flag_value);
    }
    else if (strstr(env_var_kv, "MP_FIFO_MTU") == env_var_kv) {
      set("DDM_FIFO_MTU_SIZE", flag_value);
    }
    else if (strstr(env_var_kv, "MP_BULK_MIN_MSG_SIZE") == env_var_kv) {
      set("DDM_BULK_MIN_MSG_SIZE", flag_value);
    }
    env_var_kv = *(environ + i);
  }
#ifndef DART_MPI_DISABLE_SHARED_WINDOWS
  set("DDM_ENABLE_MPI_SHWIN", true);
#endif
#ifdef DDM_ENABLE_OPENMP
  set("DDM_ENABLE_OPENMP",    true);
#endif
#ifdef DDM_ENABLE_PAPI
  set("DDM_ENABLE_PAPI",      true);
#endif
#ifdef DDM_ENABLE_HWLOC
  set("DDM_ENABLE_HWLOC",     true);
#endif
#ifdef DDM_ENABLE_NUMA
  set("DDM_ENABLE_NUMA",      true);
#endif
#ifdef DDM_ENABLE_MKL
  set("DDM_ENABLE_MKL",       true);
#endif
#ifdef DDM_ENABLE_BLAS
  set("DDM_ENABLE_BLAS",      true);
#endif
#ifdef DDM_ENABLE_LAPACK
  set("DDM_ENABLE_LAPACK",    true);
#endif
#ifdef DDM_ENABLE_SCALAPACK
  set("DDM_ENABLE_SCALAPACK", true);
#endif
#ifdef DDM_ENABLE_PLASMA
  set("DDM_ENABLE_PLASMA",    true);
#endif
#ifdef DDM_ENABLE_HDF5
  set("DDM_ENABLE_HDF5",      true);
#endif

#ifdef MPI_IMPL_ID
   set("DART_MPI_IMPL", ddm__toxstr(MPI_IMPL_ID));
#endif
}

void Config::set(
  const std::string & key,
  std::string         value)
{
  DDM_LOG_TRACE("util::Config::set(string,string)", key, value);
  Config::config_values_[key] = value;
  Config::on_change(key, value);

  // Parse boolean values:
  std::string value_lowercase = value;
  std::transform(value.begin(), value.end(), value_lowercase.begin(),
                 ::tolower);
  if (value_lowercase == "true" || value_lowercase == "yes" ||
      value_lowercase == "t"    || value_lowercase == "y"   ||
      value_lowercase == "on") {
    set_str(key + "_BOOL", "1");
    return;
  }
  else if (value_lowercase == "false" || value_lowercase == "no" ||
           value_lowercase == "f"     || value_lowercase == "n"  ||
           value_lowercase == "off") {
    set_str(key + "_BOOL", "0");
    return;
  }

  // Parse sizes from human-readable format into number of bytes,
  // e.g. 2K -> 2048.
  // Parse config keys that end in '_SIZE' with values ending in
  // 'K', 'M' or 'G'.
  std::string size_suffix = "_SIZE";
  const char * key_suffix = key.c_str()
                            + (key.length() - size_suffix.length());
  if (size_suffix == key_suffix) {
    std::vector<char> value_cstr(value.begin(), value.end());
    value_cstr.push_back('\0');
    char * value_cstr_end = &value_cstr[0] + value.length();
    int    sh = 0;
    errno = 0;
    auto   value_bytes = strtoull(value.c_str(), &value_cstr_end, 10);
    if (!errno && value_cstr_end != value.c_str()) {
      switch(*value_cstr_end) {
        case 'K': sh = 10; break;
        case 'M': sh = 20; break;
        case 'G': sh = 30; break;
        case 0:   sh = 0;  break;
        default:  sh = -1; break;
      }
      if (sh > 0 && value_bytes < (std::numeric_limits<size_t>::max()) >> sh) {
        value_bytes <<= sh;
      }
    }
    std::string key_name_bytes = key + "_BYTES";
    set(key_name_bytes, value_bytes);
  }
}

template<>
std::string Config::get<std::string>(
  const std::string & key)
{
  std::string value = get_str(key);
  return value;
}

} // namespace util
} // namespace ddm
