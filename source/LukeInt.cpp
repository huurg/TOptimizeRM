#include "LukeInt.h"

#include <iostream>
using namespace std;

#include <utility>

void LukeInt::sort(int* x, int n, bool desc, int* a, int method) {
    if(a) {
        for(int i = 0; i < n; i++) a[i]=i;
    }
    if(desc) {
        for(int i = 0; i < (n-1); i++) {
            for(int j = 0; j < (n-1-i); j++) {
                if(x[j]<x[j+1]) {
                    swap(x[j],x[j+1]);
                    if(a) {
                        swap(a[j],a[j+1]);
                    }
                }
            }
        }
    } else {
        for(int i = 0; i < (n-1); i++) {
            for(int j = 0; j < (n-1-i); j++) {
                if(x[j]>x[j+1]) {
                    swap(x[j],x[j+1]);
                    if(a) {
                        swap(a[j],a[j+1]);
                    }
                }
            }
        }
    }

}

int LukeInt::randi(int in_min, int in_max) {
    int d = in_max - in_min + 1;
    int out = rand()%d + in_min;
    return out;
}

void LukeInt::randi(int* x, int n, int in_min, int in_max) {
    for(int i = 0; i < n; i++) {
        x[i] = randi(in_min,in_max);
    }
}

void LukeInt::print(int* x, int n, const char* pre) {
    if(pre) cout << pre;
    for(int i = 0; i < n; i++) {
        cout << x[i];
        if(i!=(n-1)) {
            cout << ", ";
        }
    }
    cout << endl;
}

void LukeInt::copy(int* dst, const int* src, int n) {
    for(int i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}

void LukeInt::sub(int* dst, const int* src, int n, int m, int i0) {
    for(int i = 0; i < m; i++) {
        dst[i] = src[i+i0];
    }
}

void LukeInt::concat(int* top, const int* bottom, int n, int m) {
    for(int i = 0; i < m; i++) {
        top[i+n] = bottom[i];
    }
}

void LukeInt::randperm(int* x, int n) {
    bool* p = new bool[n];
    for(int i = 0; i < n; i++) p[i] = false;
    for(int i = 0; i < n; i++) {
        int r = LukeInt::randi(0,n-1-i);
        int shift = 0;
        for(int j = 0; j <= (r+shift); j++) {
            if(p[j]) {
                shift++;
            }
        }
        x[i] = r+shift;
        p[r+shift] = true;
    }
    delete [] p;
}
