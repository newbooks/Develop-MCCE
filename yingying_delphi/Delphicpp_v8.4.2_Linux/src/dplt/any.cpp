#include "any.h"

#ifdef __GNUG__ // gnu C++ compiler

#include <cxxabi.h>

std::string demangle(const char* mangled_name) {

    std::size_t len = 0;
    int status = 0;
    std::unique_ptr< char, decltype(&std::free) > ptr(
        __cxxabiv1::__cxa_demangle(mangled_name, nullptr, &len, &status), &std::free);
    return ptr.get();
}

#else

std::string demangle(const char* name) { return name; }

#endif // _GNUG_