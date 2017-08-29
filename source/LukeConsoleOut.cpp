#include "LukeConsoleOut.h"

#include <iostream>
using namespace std;

#include <ostream>
#include <ctime>

int LukeConsoleOut::LOut_Pad = 0;
int dout_n = 0;

ostream& LukeConsoleOut::LOut() {
    for(int i = 0; i < LOut_Pad; i++) {
        cout << "\t";
    }
    return cout;
}

void LukeConsoleOut::dout() {
    LOut() << dout_n << endl;
    dout_n++;
}

double LukeConsoleOut::secs(clock_t tic, clock_t toc) {
    double out = 0.0;

    out = ((double)toc-(double)tic)/(double)CLOCKS_PER_SEC;

    return out;
}

void LukeConsoleOut::warning(const char* message, const char* function_name, const char* class_name) {
    LOut() << "WARNING! ";
    if(function_name||class_name) {
        LOut() << "In `";
        if(class_name) {
            LOut() << class_name;
            if(function_name) {
                LOut() << "::";
            }
        }
        if(function_name) {
            LOut() << function_name;
        }
        LOut() << "'. ";
    }
    LOut() << message << endl;
}

void LukeConsoleOut::error(const char* message, const char* function_name, const char* class_name) {
    LOut() << "ERROR! ";
    if(function_name||class_name) {
        LOut() << "In `";
        if(class_name) {
            LOut() << class_name;
            if(function_name) {
                LOut() << "::";
            }
        }
        if(function_name) {
            LOut() << function_name;
        }
        LOut() << "'. ";
    }
    LOut() << message << endl;
}
