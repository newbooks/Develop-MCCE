#ifndef DDM__FUTURE_H__INCLUDED
#define DDM__FUTURE_H__INCLUDED

#include <cstddef>
#include <functional>
#include <sstream>
#include <iostream>

#include "../ddm/Exception.h"
#include "../ddm/internal/Logging.h"


namespace ddm {

template<typename ResultT>
class Future
{
private:
  typedef Future<ResultT>               self_t;
  typedef std::function<ResultT (void)> func_t;

private:
  func_t    _func;
  ResultT   _value;
  bool      _ready     = false;
  bool      _has_func  = false;

public:
  // For ostream output
  template<typename ResultT_>
  friend std::ostream & operator<<(
      std::ostream & os,
      const Future<ResultT_> & future);

public:
  Future()
  : _ready(false),
    _has_func(false)
  { }

  Future(const func_t & func)
  : _func(func),
    _ready(false),
    _has_func(true)
  { }

  Future(
    const self_t & other)
  : _func(other._func),
    _value(other._value),
    _ready(other._ready),
    _has_func(other._has_func)
  { }

  Future<ResultT> & operator=(const self_t & other)
  {
    if (this != &other) {
      _func      = other._func;
      _value     = other._value;
      _ready     = other._ready;
      _has_func  = other._has_func;
    }
    return *this;
  }

  void wait()
  {
    DDM_LOG_TRACE_VAR("Future.wait()", _ready);
    if (_ready) {
      return;
    }
    if (!_has_func) {
      DDM_LOG_ERROR("Future.wait()", "No function");
      DDM_THROW(
        ddm::exception::RuntimeError,
        "Future not initialized with function");
    }
    _value = _func();
    _ready = true;
    DDM_LOG_TRACE_VAR("Future.wait >", _ready);
  }

  bool test() const
  {
    return _ready;
  }

  ResultT & get()
  {
    DDM_LOG_TRACE_VAR("Future.get()", _ready);
    wait();
    DDM_LOG_TRACE_VAR("Future.get >", _value);
    return _value;
  }

}; // class Future

template<typename ResultT>
std::ostream & operator<<(
  std::ostream & os,
  const Future<ResultT> & future)
{
  std::ostringstream ss;
  ss << "ddm::Future<" << typeid(ResultT).name() << ">(";
  if (future._ready) {
    ss << future._value;
  } else {
    ss << "not ready";
  }
  ss << ")";
  return operator<<(os, ss.str());
}

}  // namespace ddm

#endif // DDM__FUTURE_H__INCLUDED
