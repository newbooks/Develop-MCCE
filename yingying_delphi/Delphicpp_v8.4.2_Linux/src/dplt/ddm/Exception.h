#ifndef DDM__EXCEPTION_H_
#define DDM__EXCEPTION_H_

#include "../ddm/exception/RuntimeError.h"
#include "../ddm/exception/InvalidArgument.h"
#include "../ddm/exception/OutOfRange.h"
#include "../ddm/exception/NotImplemented.h"
#include "../ddm/exception/AssertionFailed.h"
#include "../ddm/exception/StackTrace.h"

#include "../ddm/internal/Macro.h"
#include "../ddm/internal/Logging.h"

#include <sstream>

#define DDM_ENABLE_ASSERTIONS

#define DDM_STACK_TRACE() do { \
    ddm__print_stacktrace(); \
  } while(0)

#define DDM_THROW(excep_type, msg_stream) do { \
    ::std::ostringstream os; \
    os << "[ Unit " << ddm::myid() << " ] "; \
    os << msg_stream; \
    DDM_LOG_ERROR(ddm__toxstr(excep_type), os.str()); \
    throw(excep_type(os.str())); \
  } while(0)

#if defined(DDM_ENABLE_ASSERTIONS)

#define DDM_ASSERT(expr) do { \
  if (!(expr)) { \
    DDM_THROW(ddm::exception::AssertionFailed, \
               "Assertion failed: " \
               << " " << __FILE__ << ":" << __LINE__); \
  }\
} while(0)

#define DDM_ASSERT_MSG(expr, msg) do { \
  if (!(expr)) { \
    DDM_THROW(ddm::exception::AssertionFailed, \
               "Assertion failed: " << msg \
               << " " << __FILE__ << ":" << __LINE__); \
  }\
} while(0)

#define DDM_ASSERT_RETURNS(expr, exp_value) do { \
  if ((expr) != (exp_value)) { \
    DDM_THROW(ddm::exception::AssertionFailed, \
               "Assertion failed: Expected " << (exp_value) \
               << " " << __FILE__ << ":" << __LINE__); \
  }\
} while(0)

// Using (value+1) < (lower+1) to avoid compiler warning for unsigned.
#define DDM_ASSERT_RANGE(lower, value, upper, message) do { \
  if ((value) > (upper) || (value+1) < (lower+1)) { \
    DDM_THROW(ddm::exception::OutOfRange, \
               "Range assertion " \
               << lower << " <= " << value << " <= " << upper \
               << " failed: " \
               << message << " "\
               << __FILE__ << ":" << __LINE__); \
  }\
} while(0)

#define DDM_ASSERT_EQ(val_a, val_b, message) do { \
  if ((val_a) != (val_b)) { \
    DDM_THROW(ddm::exception::AssertionFailed, \
               "Assertion " \
               << val_a << " == " << val_b \
               << " failed: " \
               << message << " " \
               << __FILE__ << ":" << __LINE__); \
  } \
} while(0)

#define DDM_ASSERT_NE(val_a, val_b, message) do { \
  if ((val_a) == (val_b)) { \
    DDM_THROW(ddm::exception::AssertionFailed, \
               "Assertion " \
               << val_a << " != " << val_b \
               << " failed: " \
               << message << " " \
               << __FILE__ << ":" << __LINE__); \
  } \
} while(0)

// Using (value+1) < (min+1) to avoid compiler warning for unsigned.
#define DDM_ASSERT_GT(value, min, message) do { \
  if (((value)+1) <= ((min)+1)) { \
    DDM_THROW(ddm::exception::OutOfRange, \
               "Range assertion " \
               << value << " > " << min \
               << " failed: " \
               << message << " " \
               << __FILE__ << ":" << __LINE__); \
  } \
} while(0)

// Using (value+1) < (min+1) to avoid compiler warning for unsigned.
#define DDM_ASSERT_GE(value, min, message) do { \
  if (((value)+1) < ((min)+1)) { \
    DDM_THROW(ddm::exception::OutOfRange, \
               "Range assertion " \
               << value << " >= " << min \
               << " failed: " \
               << message << " " \
               << __FILE__ << ":" << __LINE__); \
  } \
} while(0)

// Using (value+1) >= (max+1) to avoid compiler warning for unsigned.
#define DDM_ASSERT_LT(value, max, message) do { \
  if (((value)+1) >= ((max)+1)) { \
    DDM_THROW(ddm::exception::OutOfRange, \
               "Range assertion " \
               << value << " < " << max \
               << " failed: " \
               << message << " "\
               << __FILE__ << ":" << __LINE__); \
  } \
} while(0)

// Using (value+1) > (max+1) to avoid compiler warning for unsigned.
#define DDM_ASSERT_LE(value, max, message) do { \
  if (((value)+1) > ((max)+1)) { \
    DDM_THROW(ddm::exception::OutOfRange, \
               "Range assertion " \
               << value << " <= " << max \
               << " failed: " \
               << message << " "\
               << __FILE__ << ":" << __LINE__); \
  } \
} while(0)

#else  // DDM_ENABLE_ASSERTIONS

#define DDM_ASSERT(expr) do { } while (0)
#define DDM_ASSERT_MSG(expr, msg) do { } while (0)
#define DDM_ASSERT_RETURNS(expr, exp_value) do { \
          (expr); \
          ddm__unused(exp_value); \
        } while(0)
#define DDM_ASSERT_RANGE(lower, value, upper, message) do { \
          ddm__unused(lower);   \
          ddm__unused(value);   \
          ddm__unused(upper);   \
        } while(0)
#define DDM_ASSERT_EQ(val, exp, message) do { \
          ddm__unused(val); \
          ddm__unused(exp); \
        } while (0)
#define DDM_ASSERT_NE(val, exp, message) do { \
          ddm__unused(val); \
          ddm__unused(exp); \
        } while (0)
#define DDM_ASSERT_GT(val, min, message) do { \
          ddm__unused(val); \
          ddm__unused(min); \
        } while (0)
#define DDM_ASSERT_GE(val, min, message) do { \
          ddm__unused(val); \
          ddm__unused(min); \
        } while (0)
#define DDM_ASSERT_LT(val, max, message) do { \
          ddm__unused(val); \
          ddm__unused(max); \
        } while (0)
#define DDM_ASSERT_LE(val, max, message) do { \
          ddm__unused(val); \
          ddm__unused(max); \
        } while (0)

#endif // DDM_ENABLE_ASSERTIONS

#endif // DDM__EXCEPTION_H_
