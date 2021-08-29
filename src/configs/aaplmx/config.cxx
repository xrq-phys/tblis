#include "config.hpp"

#include "util/cpuid.hpp"

namespace tblis
{

int aaplmx_check()
{
#ifndef __APPLE__
    if (get_verbose() >= 1) printf("tblis: aaplmx: Not on Apple platform.\n");
    return -1;
#endif

    return 5;
}

}
