
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstddef>
#include <sstream>
#include "dplt.h"

DPLT::DPLT(int argc, char* argv[]) { init(argc, argv); std::cout << "***DPLT is constructed***" << std::endl; }

DPLT::DPLT() { char** empty;  init(0, empty); }

void DPLT::init(int argc, char* argv[])
{
    ddm::init(&argc, &argv);
    team_size = ddm::Team::All().size();
    teamspec_3d = ddm::TeamSpec<3>(team_size, 1, 1);
    teamspec_1d = ddm::TeamSpec<1>(team_size);
    teamspec_3d.balance_extents();
    teamspec_1d.balance_extents();
}

DPLT::~DPLT(void)
{
    finalize();
}


bool DPLT::createGlobalChar(const std::string objKey, size_t size)
{
    bool success = dataMap.insert(std::make_pair(objKey, createGlobalMem<char>(size))).second;
    DDM_ASSERT_MSG(success, ("Failed to create string object " + objKey + "in global memory map!"));

    return true;
}

bool DPLT::writeGlobalChar(const std::string objKey, const std::string &send)
{
    ddm::GlobMem<char>& gString = getGlobalMemObj<char>(objKey);

    size_t i = 0;
    for (std::string::const_iterator it = send.begin(); it != send.end(); it++)
    {
        gString.put_value(*it, i);
        i++;
    }
    gString.flush_all();
    return true;
}


std::string DPLT::readGlobalChar(const std::string objKey)
{
    ddm::GlobMem<char>& gString = getGlobalMemObj<char>(objKey);
    std::string receive;
    size_t j = gString.size();
    char temp_char;
    for (size_t i = 0; i<j; i++)
    {
        gString.get_value(&temp_char, i);
        receive += temp_char;
    }
    return receive;
}