#ifndef DDM__UTIL__TIME_MEASURE_H_
#define DDM__UTIL__TIME_MEASURE_H_

namespace ddm {
namespace util {

class TimeMeasure {
 public: 
   typedef enum {
     Counter = 0,
     Clock   = 1
   } MeasureMode; 

   typedef enum {
     Cycles = 0,
     Time   = 1
   } MeasureDomain; 
};

} // namespace util
} // namespace ddm

#endif // DDM__UTIL__TIME_MEASURE_H_
