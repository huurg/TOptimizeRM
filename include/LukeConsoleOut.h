#ifndef LUKE_CONSOLE_OUT_HEADER
#define LUKE_CONSOLE_OUT_HEADER

#include <iostream>
using namespace std;

#include <ostream>

#include <ctime>

extern int dout_n;

namespace LukeConsoleOut {
    extern int LOut_Pad;

    ostream& LOut();
    void dout();

    double secs(clock_t tic, clock_t toc);
}


#endif
