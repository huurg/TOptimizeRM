#include "TO_Decoder.h"

#include <iostream>
using namespace std;

#include "GateStringSparse.h"
#include "Signature.h"
#include "GateSigInterface.h"
#include "LukeBool.h"
#include "BMSparse.h"
#include "Interface_BMSGSS.h"
#include "Interface_SigBMS.h"
#include "PhasePolynomial.h"
#include "LukeConsoleOut.h"
#include "LukeMat_GF2.h"
#include "GateSynthesisMatrix.h"
#include "TO_Maps.h"
#include "WeightedPolynomial.h"
using namespace LukeConsoleOut;
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>

string g_output_filename;
string g_csv_filename;
string g_indvar_out;
ostringstream g_h_order_stream;
string g_best_random_h_order_f, g_best_random_h_order_nf;
bool g_error_report = false;
string g_algorithm;
int g_random_circuit_seed = 0;
bool g_lempel_feedback = 1;
int g_Reed_Muller_max = 5;
int g_Hadamard_ancillas = 20;

GateStringSparse ReedMullerSynthesis2(const Signature& inS) {
    int n = inS.get_n();
    int N = (int)pow(2,n);
    GateStringSparse out(n);

    GateStringSparse daft = GateSigInterface::SigToGSS(inS);
    bool* gss_vec = new bool[N];
    for(int i = 0; i < N; i++) gss_vec[i] = daft.E(i);
    bool* result = new bool[N];
    cout << "Reed decoder begin." << endl;
    clock_t tic = clock();
    LukeBool::ReedDecoder(gss_vec,n-4,n,NULL,result);
    clock_t toc = clock();
    cout << "Reed decoder end." << endl;
    double exec_time = ((double)toc-tic)/CLOCKS_PER_SEC;
    cout << "Execution time = " << exec_time << "s" << endl;
    for(int i = 0; i < N; i++) {
        if(result[i]) {
            out.set(i);
        }
    }

    delete [] gss_vec;
    delete [] result;
    return out;
}

GateStringSparse ReedMullerSynthesis(const Signature& inS) {
    int n = inS.get_n();
    GateStringSparse out(n);

    cout << "Circuit to optimize:" << endl;
    inS.print();
    cout << endl;

    GateStringSparse daft = GateSigInterface::SigToGSS(inS);
    //daft.printString("Daft gate string: ");
    daft.print();
    int weight = daft.weight(true);
    cout << "Daft weight = " << weight << endl;
    out = daft;
    cout << endl;

    vector<GateStringSparse*> myRM = GateStringSparse::ReedMullerGenerators(n-4,n);
    vector<GateStringSparse*> optimalStrings;
    int N = myRM.size();
    double N_total = pow(2,myRM.size());
    int* weightdist = new int[(int)pow(2,n)+1];
    for(int i = 0; i < (pow(2,n)+1); i++) weightdist[i] = 0;
    cout << "Total codewords to search = " << N_total << endl;
    cout << "Search begin..." << endl;
    bool thisBS[N],lastBS[N],addBS[N];
    GateStringSparse thisCodeWord(N);
    clock_t tic = clock();
    clock_t start = clock();
    LukeBool::zeros(lastBS,N);
    LukeBool::zeros(thisBS,N);
    bool exit = false;
    double i = 0;
    while(!exit) {
        if(((clock()-tic)/(double)CLOCKS_PER_SEC)>=update_time) {
            double perc_done = i/N_total;
            double done_rate = CLOCKS_PER_SEC*perc_done/(double)(clock()-start);
            double time_left = (1-perc_done)/done_rate;
            cout << "Percentage done = " << perc_done*100 << ";\tcompletion rate = " << 100*done_rate << "%ps;\tTime elapsed = " << (clock()-start)/(double)CLOCKS_PER_SEC << "s;\tTime left = " << time_left << "s;\t";
            cout << "Completion rate (iterations) = " << i*(double)CLOCKS_PER_SEC/(double)(clock()-start) << "ps" << endl;
            tic = clock();
        }
        LukeBool::BitwiseXor(thisBS,lastBS,addBS,N);
        for(int j = 0; j < N; j++) {
            if(addBS[j]) {
                thisCodeWord = (thisCodeWord + (*myRM[j]));
            }
        }
        GateStringSparse newFunc = (daft+thisCodeWord);
        int newWeight = newFunc.weight(true);
        weightdist[newWeight]++;
        if(newWeight<weight) {
            weight = newWeight;
            out = newFunc;
            for(int i = 0; i < (int)optimalStrings.size(); i++) {
                delete optimalStrings[i];
            }
            optimalStrings.clear();
        }
        if(newWeight==weight) {
            GateStringSparse* thisOptimal = new GateStringSparse(N);
            thisOptimal->assign(newFunc);
            optimalStrings.push_back(thisOptimal);
        }
        LukeBool::copy(thisBS,lastBS,N);
        exit = !LukeBool::increment(thisBS,N);
        i++;
    }
    double total_time = (clock()-start)/(double)CLOCKS_PER_SEC;
    cout << "Search done! Total time = " << total_time << "s" << endl << endl;

    cout << "Optimal weight = " << out.weight(true) << endl;
    cout << "Weight reduced by " << (daft.weight(true)-out.weight(true)) << endl;
    cout << endl;
    cout << "Optimal gate strings:" << endl;
    for(int i = 0; i < (int)optimalStrings.size(); i++) {
        cout << "Optimal " << i << ":" << endl;
        GateStringSparse* thisGSS = optimalStrings[i];
        //thisGSS->printString();
        Signature difference(n);
        Signature this_GSS_Sig = GateSigInterface::expandGSSTerm(*thisGSS);
        difference.assign(inS + this_GSS_Sig);
        bool success = difference.isEmpty();
        cout << "Success? " << (success?"Yes":"No") << endl;
        if(!success) difference.print();
    }
    cout << endl;
    cout << "Optimal decompositions:" << endl;
    for(int i = 0; i < (int)optimalStrings.size(); i++) {
        cout << "Optimal " << i << ":" << endl;
        GateStringSparse* thisGSS = optimalStrings[i];
        thisGSS->print();
    }



    cout << endl;

    //cout << "Weight distribution:" << endl;
    //for(int i = 0; i < (pow(2,n)+1); i++) cout << "N(w = " << i << ") = " << weightdist[i] << endl;

    for(int i = 0; i < (int)myRM.size(); i++) {
        delete myRM[i];
    }
    for(int i = 0; i < (int)optimalStrings.size(); i++) {
        delete optimalStrings[i];
    }
    delete [] weightdist;
    weightdist = NULL;

    if((!g_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::REED_MULLER))&&(!g_output_filename.empty())) {
        ofstream my_file(g_output_filename.c_str(),iostream::app);
        my_file << "Input signature" << endl;
        inS.print(my_file);
        my_file << "Output gate string" << endl;
        out.print(my_file);
        result_analysis(inS,out,my_file);
        my_file << "Total execution time = " << total_time << "s" << endl;
        my_file << endl;
        my_file << "Input signature as file" << endl;
        inS.save(my_file);
        my_file << endl;
        my_file << "Output gate synthesis matrix" << endl;
        BMSparse out_BMS = Interface_BMSGSS::GSSToBMS(out);
        out_BMS.printFull(my_file);
        my_file.close();
    }
    if((!g_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::REED_MULLER))&&(!g_csv_filename.empty())&&(g_random_circuit_seed!=0)) {
        ofstream csv_out_file(g_csv_filename.c_str(), iostream::app);
        //CSV format: "2, <#qubits>, <seed#>, <T-count>\n"

        csv_out_file << SYNTHESIS_ALGORITHM_ID::REED_MULLER << "," << n << "," << g_random_circuit_seed << "," << out.weight(true) << endl;

        csv_out_file.close();
    }

    return out;
}

void result_analysis(const Signature& inS, const GateStringSparse& inResult, ostream& inOS) {
    GateStringSparse daft = GateSigInterface::SigToGSS(inS);
    int daft_weight = daft.weight(true);
    int result_weight = inResult.weight(true);
    //Signature restored_sig = GateSigInterface::expandGSSTerm(inResult);
    Signature restored_sig = Interface_SigBMS::BMSToSig(Interface_BMSGSS::GSSToBMS(inResult));
    bool success = (inS + restored_sig).isEmpty();
    inOS << "Successful? " << (success?"Yes":"No") << endl;
    if(!success) {
        Signature difference = inS + restored_sig;
        cout << "Difference:" << endl;
        difference.print();
    }
    inOS << "T-fast? " << (T_fast(inResult)?"Yes":"No") << endl;
    inOS << "Daft T-count = " << daft_weight << endl;
    inOS << "Output T-count = " << result_weight << endl;
    inOS << "Reduced T-count by " << (daft_weight-result_weight) << endl;
    inOS << endl;
}

bool synthesis_success(const Signature& inS, const GateStringSparse& inResult) {
    bool out;

    //Signature restored_sig = GateSigInterface::expandGSSTerm(inResult);
    Signature restored_sig = Interface_SigBMS::BMSToSig(Interface_BMSGSS::GSSToBMS(inResult));
    Signature difference = inS + restored_sig;
    out = difference.isEmpty();

    return out;
}

bool T_fast(const GateStringSparse& inGSS) {
    bool out;

    int k = inGSS.get_n();
    int T = inGSS.weight(true);
    out = (T<=((double)(k*k + 3*k - 14)/2.0));

    return out;
}

GateStringSparse LempelXSynthesis(const Signature& inS) {
    LOut(); cout << "Gate synthesis begin." << endl;
    clock_t start = clock();
    GateStringSparse A_GSS = GateSigInterface::SigToGSS(inS);
    BMSparse A_BMS = Interface_BMSGSS::GSSToBMS(A_GSS);
    int n = A_BMS.get_n();
    int m = A_BMS.get_m();
    bool** A_bool = LukeMat_GF2::construct(n,m);
    A_BMS.toBool(A_bool);
    int t;
    clock_t tic = clock();
    LOut_Pad++;
    if(!g_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_X))
        GateSynthesisMatrix::LempelX(A_bool,n,m,t);
    else if(!g_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_X_2))
        GateSynthesisMatrix::LempelX2(A_bool,n,m,t);
    else
        GateSynthesisMatrix::LempelX(A_bool,n,m,t);
    clock_t toc = clock();
    double exec_time = ((double)toc-(double)tic)/CLOCKS_PER_SEC;
    BMSparse out_BMS(n,t);
    out_BMS.fromBool(A_bool,n,t);
    GateStringSparse out = Interface_BMSGSS::BMSToGSS(out_BMS);
    LukeMat_GF2::destruct(A_bool,n,m);
    LOut() << "LempelX executed in " << exec_time << " seconds." << endl;
    clock_t finish = clock();
    LOut() << "Total time: " << secs(start,finish) << "s" << endl;
    LOut_Pad--;
    LOut() << "Gate synthesis end." << endl;
    if(!g_csv_filename.empty()) {
        ofstream my_csv(g_csv_filename.c_str(),iostream::app);
        int n = inS.get_n();
        my_csv << SYNTHESIS_ALGORITHM_ID::DAFT_GUESS << "," << n << "," << g_random_circuit_seed << "," << A_GSS.weight(true) << endl;
        my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_X << "," << n << "," << g_random_circuit_seed << "," << out.weight(true) << endl;
        my_csv.close();
    }
    return out;
}

GateStringSparse LempelXSynthesis2(const Signature& inS) {
    int n = inS.get_n();

    GateStringSparse out(n);

    GateStringSparse daft = GateSigInterface::SigToGSS(inS);
    int init_weight = daft.weight(true);
    // Step 0
    BMSparse this_B = Interface_BMSGSS::GSSToBMS(daft);
    cout << "Initial B" << endl;
    this_B.printFull();

    cout << "Initial B, randomly permuted" << endl;
    this_B = (this_B*BMSparse::random_perm(this_B.get_m()));

    int synth_round = 0;
    bool exit_outer = false;

    cout << "LempelX synthesis begin." << endl;
    clock_t t_begin = clock();

    int T_count_phase_2_reduction = 0;

    while(!exit_outer) {
        cout << "Synthesis round: " << synth_round << endl;
        bool exit_inner = false;
        int m = this_B.get_m();





        int c = 0;
        int col_pairs_skipped = 0;
        BMSparse unique_xs(n,0);
        for(int i = 0; (!exit_inner)&&(i < (m-1)); i++) {
            for(int j = (i+1); (!exit_inner)&&(j < m); j++, c++) {
                cout << "(i, j) = (" << i << ", " << j << ")" << endl;
                cout << "Col pair " << (c+1) << " of " << 0.5*m*(m-1) << endl;
                // Step 1
                BMSparse this_x = (this_B.col(i) + this_B.col(j));
                cout << "x = c_i + c_j" << endl;
                //this_x.printFull();

                // Check if this_x has been used in previous (i,j) iteration
                bool x_found = false;
                for(int k = 0; (!x_found)&&(k < unique_xs.get_m()); k++) {
                    if(this_x==unique_xs.col(k)) {
                        x_found = true;
                    }
                }

                // Must rewrite this so all column pairs that sum to x are checked against N_i + N_j = 1 (mod 2)
                if(false&&x_found) {
                    col_pairs_skipped++;
                    cout << "This x skipped. Total skipped = " << col_pairs_skipped << endl;
                } else {
                    unique_xs = (unique_xs&&this_x);

                    // Step 2: Extend gate synthesis matrix with additional rows x_a*(r_b^r_c)+...
                    BMSparse this_B_sq = this_B.GSMatrixExtension2(this_x);
                    cout << "B'" << endl;
                    //this_B_sq.printFull();

                    // Step 3: Calculate nullspace
                    BMSparse this_N = this_B_sq.nullspace();
                    cout << "N" << endl;
                    //this_N.printFull();

                    // Step 4: Determine if N_i+N_j = 1 (mod 2)
                    int p = this_N.get_m();
                    int this_k;
                    for(int k = 0; (!exit_inner)&&(k < p); k++) {
                        exit_inner = (this_N.E(i,k)+this_N.E(j,k))%2;
                        if(exit_inner) this_k = k;
                    }

                    if(exit_inner) {
                        // Step 5: Add x to cols in nullspace portion of original matrix
                        BMSparse this_Nk = this_N.col(this_k);
                        cout << "this_Nk" << endl;
                        //this_Nk.T().printFull();

                        if(this_N.E(i,this_k)) { // Forces i to be a column in F and j to be a column in G
                            swap(i,j);
                        }
                        // Construct F
                        BMSparse this_F(n,0);
                        this_F = (this_F && this_B.col(i));
                        for(int k = 0; k < m; k++) {
                            if((k!=i)&&(!this_Nk.E(k,0))) {
                                this_F = (this_F && this_B.col(k));
                            }
                        }
                        cout << "F" << endl;
                        //this_F.printFull();
                        // Construct G
                        BMSparse this_G(n,0);
                        this_G = (this_G && this_B.col(j));
                        for(int k = 0; k < m; k++) {
                            if((k!=j)&&(this_Nk.E(k,0))) {
                                this_G = (this_G && this_B.col(k));
                            }
                        }
                        cout << "G" << endl;
                        //this_G.printFull();

                        BMSparse this_Z(n,(this_G.get_m()%2)?(this_G.get_m()+1):this_G.get_m());
                        this_Z.sub(0,0,this_G);
                        cout << "Z" << endl;
                        //this_Z.printFull();

                        /*
                        // Check sig(Z)
                        Signature this_Z_sig = GateSigInterface::expandGSSTerm(Interface_BMSGSS::BMSToGSS(this_Z));
                        */

                        BMSparse this_Z_sq;
                        this_Z_sq.assign(this_Z);
                        for(int k = 0; k < this_Z_sq.get_m(); k++) {
                            this_Z_sq.sub(0,k,(this_Z.col(k)+this_x));
                        }
                        cout << "Z~" << endl;
                        //this_Z_sq.printFull();

                        /*
                        // Check sig(Z~)
                        Signature this_Z_sq_sig = GateSigInterface::expandGSSTerm(Interface_BMSGSS::BMSToGSS(this_Z_sq.sub(0,1,-1,-1)));

                        Signature Z_diff = (this_Z_sig+this_Z_sq_sig);
                        cout << "Z -> Z~ same signature? " << (Z_diff.isEmpty()?"Yes":"No") << endl;
                        Z_diff.print();
                        */

                        Signature this_B_before_sig = Interface_SigBMS::BMSToSig(this_B);

                        // Step 6: Update gate synthesis matrix
                        this_B = (this_F.sub(0,1,-1,-1) && this_Z_sq.sub(0,1,-1,-1));
                        cout << "new B" << endl;
                        this_B.printFull();

                        Signature this_B_after_sig = Interface_SigBMS::BMSToSig(this_B);
                        Signature this_B_diff = (this_B_before_sig+this_B_after_sig);
                        cout << "B same signature? " << (this_B_diff.isEmpty()?"Yes":"No") << endl;
                        this_B_diff.print();

                        /*
                        Signature this_F_sig = GateSigInterface::expandGSSTerm(Interface_BMSGSS::BMSToGSS(this_F.sub(0,1,-1,-1)));
                        Signature this_full_sig = (this_F_sig+this_Z_sq_sig);
                        Signature this_full_diff = (this_full_sig+this_B_before_sig);
                        cout << "Full same signature? " << (this_full_diff.isEmpty()?"Yes":"No") << endl;
                        this_full_diff.print();

                        cout << "Compare signature: ";
                        switch(this_B.compare_signature(this_B_before)) {
                            case -1:
                                cout << "N/a" << endl;
                                break;
                            case 0:
                                cout << "Same" << endl;
                                break;
                            case 1:
                                cout << "Not same" << endl;
                                break;
                        }
                        */

                        /*
                        cout << "Eliminate columns " << i << " and " << j << endl;

                        if((this_Nk.sum())%2) { // If nullspace vector has odd weight, pad B with all zeros column and Nk with extra 1;
                            cout << "Pad with zeros column" << endl;
                            BMSparse temp_B(n,m+1);
                            temp_B.sub(0,0,this_B);
                            this_B.assign(temp_B);
                            BMSparse temp_Nk(m+1,1);
                            temp_Nk.sub(0,0,this_Nk);
                            temp_Nk.E(m,0,1);
                            this_Nk.assign(temp_Nk);
                        }
                        //cout << "asdf" << endl;
                        cout << "this_Nk" << endl;
                        this_Nk.T().printFull();
                        BMSparse temp_B = this_B + this_x*(this_Nk.T());
                        cout << "temp_B" << endl;
                        temp_B.printFull();
                        //cout << "zxcv" << endl;

                        // Step 6
                        this_B = temp_B.sub(0,0,-1,fmin(i,j)) && temp_B.sub(0,fmin(i,j)+1,-1,fmax(i,j)-fmin(i,j)-1) && temp_B.sub(0,fmax(i,j)+1,-1,-1);
                        cout << "Elimination complete." << endl;
                        */
                    }
                }
            }
        }

        int T_count_before_phase_2 = this_B.get_m();

        // If no columns could be eliminated using x = c_i + c_j, try x = c_i (with corresponding conditions)
        for(int i = 0; (!exit_inner)&&(i < m); i++) {
            cout << "Col " << (i+1) << " of " << m << endl;
            // Step 1: Choose an x = c_i
            BMSparse this_x = this_B.col(i);
            cout << "x = c_i" << endl;
            this_x.T().printFull();
            // Step 2: Extend B to find B~
            BMSparse this_B_sq = this_B.GSMatrixExtension2(this_x);
            cout << "B'" << endl;
            this_B_sq.T().printFull();
            // Step 3: Calculate nullspace(B~)
            BMSparse this_N = this_B_sq.nullspace();
            cout << "N" << endl;
            this_N.T().printFull();
            // Step 4: Search nullspace generators for a nullspace vector with even weight, or the sum of two with odd weight, that contain the column c_i.
            BMSparse this_Nv;

            {
                int N_g = this_N.get_m();
                int N_N = (int)pow(2,N_g);
                for(int p = 0; (!exit_inner)&&(p < N_N); p++) {
                    cout << "p = " << p << endl;
                    BMSparse this_m = BMSparse::integer(p,N_g);
                    cout << "this_m" << endl;
                    this_m.T().printFull();
                    BMSparse this_n = this_N*this_m;
                    cout << "this_n" << endl;
                    this_n.T().printFull();
                    if(this_n.E(i,0)&&!(this_n.sum()%2)) {
                        this_Nv.assign(this_n);
                        exit_inner = true;
                    }
                }
            }

            /*{
                // Partition N into (E D) where E and contain all even and odd weight columns of N, respectively.
                BMSparse this_E(m,0);
                BMSparse this_D(m,0);
                for(int k = 0; k < this_N.get_m(); k++) {
                    if(this_N.col(k).sum()%2) {
                        this_D = (this_D && this_N.col(k));
                    } else {
                        this_E = (this_E && this_N.col(k));
                    }
                }

                cout << "E" << endl;
                this_E.T().printFull();
                cout << "D" << endl;
                this_D.T().printFull();


                BMSparse this_ED = (this_E && this_D);

                // Search through all nullspace vectors with even weight
                int N_p = this_E.get_m();
                int N_q = this_D.get_m();
                for(int p = 0; (!exit_inner)&&(p < pow(2,N_p)); p++) {
                    BMSparse this_m_E; BMSparse this_n_E(m,1);
                    if(N_p) {
                        this_m_E.assign(BMSparse::integer(p,N_p));
                        this_n_E.assign(this_E*this_m_E);
                    }
                    BMSparse this_n_D(m,1);

                    bool* this_m_D_ba = NULL;
                    if(N_q>=2) this_m_D_ba = new bool[N_q];

                    for(int h = 0; (!exit_inner)&&(h <= (int)(N_q/2)); h++) {
                        cout << "h = " << h << endl;
                        if(N_q>=2) for(int y = 0; y < N_q; y++) this_m_D_ba[y] = (y<(2*h));
                        bool contloop = true;
                        int this_perm = 0;
                        while((!exit_inner)&&contloop) {
                            cout << "Permutation number = " << this_perm << endl;
                            this_perm++;
                            if(N_q>=2) LukeBool::print(this_m_D_ba,N_q,"m_D: ");
                            if(N_q>=2) {
                                BMSparse this_m_D(N_q,1);
                                for(int y = 0; y < N_q; y++) {
                                    this_m_D.E(y,0,this_m_D_ba[y]);
                                }
                                this_n_D = (this_D*this_m_D);
                            }

                            cout << "this_n_E" << endl;
                            this_n_E.T().printFull();
                            cout << "this_n_D" << endl;
                            this_n_D.T().printFull();


                            BMSparse this_n = (this_n_E+this_n_D);
                            // TODO: Remove "!(this_n.sum()%2)&&" from below and try to make sure the weight is even anyway.
                            if(this_n.E(i,0)) { // If one is found where n_i = 1, use that nullspace vector and exit
                                this_Nv.assign(this_n);
                                exit_inner = true;
                            }

                            if(N_q>=2) {
                                bool* temp = new bool[N_q];
                                contloop = LukeBool::nextUniquePerm(temp,this_m_D_ba,N_q);
                                LukeBool::copy(temp,this_m_D_ba,N_q);
                                delete [] temp;
                            } else {
                                contloop = false;
                            }
                        }
                    }

                    if(N_q) delete [] this_m_D_ba;


                }
            }*/
            cout << "N_v" << endl;
            this_Nv.T().printFull();

            if(exit_inner) {
                // Step 5: Add x to each column in nullspace portion of this_B
                BMSparse this_B_p = (this_B+this_x*(this_Nv.T()));

                cout << "B + x*N_v'" << endl;
                this_B_p.printFull();

                // Step 6: Update gate synthesis matrix to exclude column c_i
                Signature this_B_before_sig = Interface_SigBMS::BMSToSig(this_B);
                cout << "old B" << endl;
                this_B.printFull();
                this_B = (this_B_p.sub(0,0,-1,i)&&this_B_p.sub(0,i+1,-1,-1));
                cout << "new B" << endl;
                this_B.printFull();

                Signature this_B_after_sig = Interface_SigBMS::BMSToSig(this_B);
                Signature this_B_diff = (this_B_before_sig+this_B_after_sig);
                cout << "B same signature? " << (this_B_diff.isEmpty()?"Yes":"No") << endl;
                this_B_diff.print();
            }

        }
        int T_count_after_phase_2 = this_B.get_m();

        T_count_phase_2_reduction += T_count_before_phase_2-T_count_after_phase_2;

        // If still no columns eliminated, exit outer loop
        if(!exit_inner) exit_outer = true;
        cout << "B" << endl;
        this_B.printFull();
        Signature this_sig = Interface_SigBMS::BMSToSig(this_B);
        Signature this_dif = (inS+this_sig);
        bool same_sig = this_dif.isEmpty();
        cout << "Same sig? " << (same_sig?"Yes":"No") << endl;
        if(!same_sig) {
            cout << "Difference" << endl;
            this_dif.print();
        }
        synth_round++;
    }




    clock_t t_end = clock();
    double exec_time = (double)(t_end-t_begin)/CLOCKS_PER_SEC;
    cout << "LempelX synthesis end. Total time = " << exec_time << "s." << endl;

    out = Interface_BMSGSS::BMSToGSS(this_B);
    int final_weight = out.weight(true);

    cout << "Weight reduced by " << (init_weight-final_weight) << endl;
    cout << "Reduction due to phase 2 = " << T_count_phase_2_reduction << endl;

    ofstream* my_file = NULL;
    if(!g_output_filename.empty()) {
        my_file = new ofstream(g_output_filename.c_str(),iostream::app);
        (*my_file) << "Input signature:" << endl;
        inS.print(*my_file);
        (*my_file) << "Output gate string:" << endl;
        out.print(*my_file);
        result_analysis(inS,out,*my_file);
        (*my_file) << "Total execution time = " << exec_time << "s" << endl << endl;
        (*my_file) << "Signature file:" << endl;
        inS.save(*my_file);
        (*my_file) << endl << "Gate synthesis matrix:" << endl;
        this_B.printFull(*my_file);
        my_file->close();
        delete my_file;
    }
    if((!g_csv_filename.empty())&&(g_random_circuit_seed!=0)) {
        ofstream csv_out_file(g_csv_filename.c_str(), iostream::app);
        //CSV format: "1, <#qubits>, <seed#>, <T-count>\n"

        csv_out_file << SYNTHESIS_ALGORITHM_ID::LEMPEL_X << "," << n << "," << g_random_circuit_seed << "," << final_weight << endl;

        csv_out_file.close();
    }

    return out;
}

PhasePolynomial FullDecoderWrapper(const PhasePolynomial& in, TO_Decoder decoder) {
    clock_t start = clock();
    int n = in.get_n();
    PhasePolynomial out(n);

    // f(x) = in = f'(x) + 2f~(x)
    clock_t tic = clock();
    WeightedPolynomial f = TO_Maps::PhasePolynomial_to_WeightedPolynomial(in); //Expensive line
    clock_t toc = clock();
    //cout << "Time1 = " << secs(tic,toc) << "s" << endl;

    WeightedPolynomial f_prime(f);

    WeightedPolynomial two_f_tilde(f);
    f_prime %= 2;
    two_f_tilde -= f_prime;

    Signature f_prime_sig = TO_Maps::WeightedPolynomial_to_Signature(f);

    //LOut() << "FullDecoderWrapper: Before decoder time = " << secs(start,clock()) << "s" << endl;

    GateStringSparse f_prime_GSS_optimized = decoder(f_prime_sig);

    clock_t mid = clock();

    tic = clock();
    PhasePolynomial g = TO_Maps::GateStringSparse_to_PhasePolynomial(f_prime_GSS_optimized);
    toc = clock();
    //cout << "Time4 = " << secs(tic,toc) << "s" << endl;

    g.print();

    out = g;

    tic = clock();
    out += TO_Maps::WeightedPolynomial_to_PhasePolynomial(two_f_tilde);
    toc = clock();
    //cout << "Time5 = " << secs(tic,toc) << "s" << endl;

    tic = clock();
    WeightedPolynomial two_g_tilde = TO_Maps::PhasePolynomial_to_WeightedPolynomial(g); // Expensive line
    toc = clock();
    //cout << "Time6 = " << secs(tic,toc) << "s" << endl;

    two_g_tilde -= f_prime;

    tic = clock();
    out -= TO_Maps::WeightedPolynomial_to_PhasePolynomial(two_g_tilde);
    toc = clock();
    //cout << "Time5 = " << secs(tic,toc) << "s" << endl;

    out.mod8();

    clock_t finish = clock();
    //LOut() << "FullDecoderWrapper: After decoder time = " << secs(mid,finish) << "s" << endl;
    //LOut() << "FullDecoderWrapper: Total time = " << secs(start,finish) << "s" << endl;
    return out;
}
