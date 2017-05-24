#include "LukeConsoleOut.h"

#include <iostream>
using namespace std;

#include <ostream>

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
