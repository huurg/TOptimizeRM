#include "LukeConsoleOut.h"

#include <iostream>
using namespace std;

#include <ostream>

int LukeConsoleOut::LOut_Pad = 0;

ostream& LukeConsoleOut::LOut() {
    for(int i = 0; i < LOut_Pad; i++) {
        cout << "\t";
    }
    return cout;
}
