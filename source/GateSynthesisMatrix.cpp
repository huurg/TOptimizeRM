#include "GateSynthesisMatrix.h"

#include <iostream>
using namespace std;

#include <cmath>
#include "LukeMat_GF2.h"
#include "LukeConsoleOut.h"
using namespace LukeConsoleOut;

bool** GateSynthesisMatrix::from_signature(bool*** S, int n, int& mp) {
    int N = (int)pow(2,n);
    bool* a = new bool[N];
    for(int i = 0; i < N; i++) a[i] = 0;

    int I;
    for(int i = 0; i < n; i++) {
        if(S[i][i][i]) {
            I = (int)pow(2,i);
            a[I] = !a[I];
        }
        if(i<(n-1)) {
            for(int j = (i+1); j < n; j++) {
                if(S[i][j][j]) {
                    I = (int)pow(2,i);
                    a[I] = !a[I];
                    I = (int)pow(2,j);
                    a[I] = !a[I];
                    I = (int)pow(2,i)+(int)pow(2,j);
                    a[I] = !a[I];
                }
                if(j<(n-1)) {
                    for(int k = (j+1); k < n; k++) {
                        if(S[i][j][k]) {
                            I = (int)pow(2,i);
                            a[I] = !a[I];
                            I = (int)pow(2,j);
                            a[I] = !a[I];
                            I = (int)pow(2,k);
                            a[I] = !a[I];
                            I = (int)pow(2,i)+(int)pow(2,j);
                            a[I] = !a[I];
                            I = (int)pow(2,i)+(int)pow(2,k);
                            a[I] = !a[I];
                            I = (int)pow(2,j)+(int)pow(2,k);
                            a[I] = !a[I];
                            I = (int)pow(2,i)+(int)pow(2,j)+(int)pow(2,k);
                            a[I] = !a[I];
                        }
                    }
                }
            }
        }
    }
    int m = 0;
    for(int i = 0; i < N; i++) {
        if(a[i]) {
            m++;
        }
    }
    bool** out = LukeMat_GF2::construct(n,m);
    int j = 0;
    for(int I = 0; I < N; I++) {
        if(a[I]) {
            int Ip = I;
            int i = 0;
            while(Ip>0) {
                out[i][j] = (Ip%2);
                Ip /= 2;
                i++;
            }
            j++;
        }
    }

    delete [] a;

    mp = m;
    return out;
}

void GateSynthesisMatrix::cleanup(bool** A, int n, int m, int& mp) {
    for(int j1 = 0; j1 < (m-1); j1++) {
        for(int j2 = (j1+1); j2 < m; j2++) {
            bool same = 1;
            for(int i = 0; same&&(i < n); i++) {
                same *= (A[i][j1]==A[i][j2]);
            }
            if(same) {
                //cout << "j1 = " << j1 << " and j2 = " << j2 << " are the same." << endl;
                for(int i = 0; i < n; i++) {
                    A[i][j1]=0;
                    A[i][j2]=0;
                }
            }
        }
    }

    int j_end = (m-1);
    int non_zero_count = 0;
    for(int j = 0; j < m; j++) {
        int sum = 0;
        for(int i = 0; i < n; i++) {
            sum += A[i][j];
        }
        if(!sum) {
            bool found = false;
            while((!found)&&(j_end>j)) {
                for(int i = 0; (!found)&&(i < n); i++) {
                    found = A[i][j_end];
                }
                if(!found) j_end--;
            }
            if(found) {
                LukeMat_GF2::swapcol(A,n,m,j,j_end);
                non_zero_count++;
            }
        } else {
            non_zero_count++;
        }
    }
    mp = non_zero_count;

}

void GateSynthesisMatrix::LempelX(bool** A, int n, int m, int& omp) {
    int this_m = m;
    bool** x = LukeMat_GF2::construct(n,1);
    bool** nv = LukeMat_GF2::construct(1,m+1);
    bool** xnv = LukeMat_GF2::construct(n,m+1);
    bool** A_xnv = LukeMat_GF2::construct(n,m+1);
    int n_ext = n+((n*(n-1)*(n-2))/6);
    bool** A_ext = LukeMat_GF2::construct(n_ext,m+1);

    bool found = 1;
    int round = 0;
    while(found&&(round<m)) {
        found = 0;
        LOut(); cout << "Round = " << round << endl;
        LukeMat_GF2::copy(A,n,this_m,A_ext);
        for(int j1 = 0; (!found)&&(j1 < (this_m-1)); j1++) {
            for(int j2 = (j1+1); (!found)&&(j2 < this_m); j2++) {
                //cout << "Col pair = (" << j1 << ", " << j2 << ")" << endl;
                for(int i = 0; i < n; i++) {
                    x[i][0] = (A[i][j1] + A[i][j2])%2;
                }
                //LukeMat_GF2::print(x,n,1,"x: ");
                int I = 0;
                for(int a = 0; a < (n-2); a++) {
                    for(int b = (a+1); b < (n-1); b++) {
                        for(int c = (b+1); c < n; c++) {
                            for(int j = 0; j < this_m; j++) {
                                A_ext[n+I][j] = (x[a][0]*A[b][j]*A[c][j] + x[b][0]*A[c][j]*A[a][j] + x[c][0]*A[a][j]*A[b][j])%2;
                            }
                            I++;
                        }
                    }
                }
                //cout << "I = " << I << endl;
                //LukeMat_GF2::print(A_ext,n_ext,this_m,"A_ext: ");
                int d;
                bool** NS = LukeMat_GF2::nullspace(A_ext,n_ext,this_m,d);
                //LukeMat_GF2::print(NS,this_m,d,"NS: ");
                /*bool** gha = LukeMat_GF2::construct(n_ext,d);
                LukeMat_GF2::times(A_ext,NS,n_ext,m,d,gha);
                LukeMat_GF2::print(gha,n_ext,d,"gha: ");
                LukeMat_GF2::destruct(gha,n_ext,d);*/
                found = 0;
                int nsv = -1;
                for(int h = 0; (!found)&&(h < d); h++) {
                    found = (NS[j1][h]+NS[j2][h])%2;
                    if(found) nsv = h;
                }
                if(found) {
                    int weight_nv = 0;
                    for(int h = 0; h < this_m; h++) {
                        nv[0][h] = NS[h][nsv];
                        weight_nv += nv[0][h];
                    }
                    if(weight_nv%2) {
                        nv[0][this_m]=1;
                        for(int h = 0; h < n; h++) A[h][this_m]=0;
                        this_m++;
                    }
                    LukeMat_GF2::times(x,nv,n,1,this_m,xnv);
                    //LukeMat_GF2::print(xnv,n,this_m,"xnv: ");
                    LukeMat_GF2::add(A,xnv,n,this_m,A_xnv);
                    LukeMat_GF2::copy(A_xnv,n,this_m,A);
                    //LukeMat_GF2::print(A,n,this_m, "A: ");
                    int mp;
                    GateSynthesisMatrix::cleanup(A,n,this_m,mp);
                    this_m = mp;
                    //LukeMat_GF2::print(A,n,this_m, "A: ");
                }
                if(d) LukeMat_GF2::destruct(NS,this_m,d);

            }
        }
        round++;
    }

    LukeMat_GF2::destruct(x,n,1);
    LukeMat_GF2::destruct(A_ext,n_ext,m+1);
    LukeMat_GF2::destruct(nv,1,m+1);
    LukeMat_GF2::destruct(xnv,n,m+1);
    LukeMat_GF2::destruct(A_xnv,n,m+1);
    omp = this_m;
}
