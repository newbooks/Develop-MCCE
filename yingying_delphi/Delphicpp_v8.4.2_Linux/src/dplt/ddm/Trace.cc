#include "../ddm/util/Trace.h"
#include "../ddm/util/Config.h"
#include "../ddm/Team.h"

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include <unistd.h>

std::map<std::string, ddm::util::TraceStore::trace_events_t>
ddm::util::TraceStore::_traces
  = {{ }};

bool ddm::util::TraceStore::_trace_enabled
  = false;

bool ddm::util::TraceStore::on()
{
  _trace_enabled = ddm::util::Config::get<bool>("DDM_ENABLE_TRACE");
  return _trace_enabled;
}

void ddm::util::TraceStore::off()
{
  _trace_enabled = false;
}

bool ddm::util::TraceStore::enabled()
{
  return _trace_enabled == true;
}

void ddm::util::TraceStore::clear()
{
  _traces.clear();
}

void ddm::util::TraceStore::clear(const std::string & context)
{
  _traces[context].clear();
}

void ddm::util::TraceStore::add_context(const std::string & context)
{
  if (_traces.count(context) == 0) {
    _traces[context] = trace_events_t();
  }
}

ddm::util::TraceStore::trace_events_t &
ddm::util::TraceStore::context_trace(const std::string & context)
{
  return _traces[context];
}

void ddm::util::TraceStore::write(std::ostream & out)
{
  if (!ddm::util::Config::get<bool>("DDM_ENABLE_TRACE")) {
    return;
  }

  std::ostringstream os;
  auto unit   = ddm::Team::GlobalUnitID();
  auto nunits = ddm::size();
  for (auto context_traces : _traces) {
    std::string      context = context_traces.first;
    trace_events_t & events  = context_traces.second;

    // Master prints CSV headers:
    if (unit == 0) {
      os << "-- [TRACE] "
         << std::setw(15) << "context"  << ","
         << std::setw(5)  << "unit"     << ","
         << std::setw(15) << "start"    << ","
         << std::setw(15) << "end"      << ","
         << std::setw(12) << "state"
         << std::endl;
    }

    ddm::barrier();

    for (auto state_timespan : events) {
      auto   start    = state_timespan.start;
      auto   end      = state_timespan.end;
      auto   state    = state_timespan.state;
      os << "-- [TRACE] "
         << std::setw(15) << std::fixed << context  << ","
         << std::setw(5)  << std::fixed << unit     << ","
         << std::setw(15) << std::fixed << start    << ","
         << std::setw(15) << std::fixed << end      << ","
         << std::setw(12) << std::fixed << state
         << std::endl;
    }
  }
  // Print trace events of units sequentially:
  for (size_t trace_unit = 0; trace_unit < nunits; ++trace_unit) {
    if (trace_unit == static_cast<size_t>(unit)) {
      out << os.str();
    }
    ddm::barrier();
  }

  // To help synchronization of writes:
  sleep(1);

  ddm::barrier();
}

void ddm::util::TraceStore::write(
  const std::string & filename,
  const std::string & path)
{
    if (!ddm::util::Config::get<bool>("DDM_ENABLE_TRACE")) {
    return;
  }

  std::string trace_log_dir;
  if (ddm::util::Config::is_set("DDM_TRACE_LOG_PATH")) {
    trace_log_dir = ddm::util::Config::get<std::string>(
                      "DDM_TRACE_LOG_PATH");
    if (path.length() > 0) {
      trace_log_dir += "/";
    }
  }
  trace_log_dir += path;

  auto unit = ddm::Team::GlobalUnitID();
  std::ostringstream fn;
  fn << trace_log_dir << "/"
     << "u" << std::setfill('0') << std::setw(5) << unit
     << "." << filename;
  std::string trace_file = fn.str();
  std::ofstream out(trace_file);
  write(out);
  out.close();
}
