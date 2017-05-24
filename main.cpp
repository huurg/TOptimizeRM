#include <iostream>
using namespace std;

#include "Signature.h"
#include "GateStringSparse.h"
#include "GateSigInterface.h"
#include "LukeMaths.h"
#include "LukeBool.h"
#include "BoolMat.h"
#include "BMSparse.h"
#include "Interface_SigBMS.h"
#include "Interface_BMSGSS.h"
#include "LukeInt.h"
#include "Bool_Signature.h"
#include "LukeMat_GF2.h"
#include "GateSynthesisMatrix.h"
#include "SQC_Circuit.h"
#include "LukeConsoleOut.h"
using namespace LukeConsoleOut;
//using namespace LinBox;
//using namespace std;

#include <climits>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <utility>

// TODO (Refactoring): Consolidate all 'interface' namespaces into a single 'converter' namespace.

typedef int (*LempelSelector) (const Signature& inS);

const double update_time = 5.0;

GateStringSparse ReedMullerSynthesis(const Signature& inS);
GateStringSparse ReedMullerSynthesis2(const Signature& inS);
GateStringSparse findNearOptimal(const Signature& inS, int maxRM=5);
bool T_fast(const GateStringSparse& inGSS);
void NullSpaceTest();
void Sig2BMSTest();
void compositionTest();
void hdecomtest();
void debuggingTest();
void n3benchmark(int n, int maxRM=5, int randgen = 0);
void n3benchmark_d(int n, int d, int maxRM=5, int randgen = 0);

//Lempel selectors
int LempelSelector_LeastGreedy(const Signature& inS);
int LempelSelector_Greedy(const Signature& inS);
int LempelSelector_Random(const Signature& inS);

//Main Lempel synthesis functions
GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, LempelSelector lempelSelector, bool feedback, int psrmc=0, int col_order=0, int n_orders=1);
GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, bool feedback, int n_rand = 1, int psrmc=0, int col_order=0, int n_orders=1);
GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, int n_rand = 1, int psrmc=0, int col_order=0, int n_orders=1);
void result_analysis(const Signature& inS, const GateStringSparse& inResult, ostream& inOS = cout);
bool synthesis_success(const Signature& inS, const GateStringSparse& inResult);

//LempelX synthesis functions
GateStringSparse LempelXSynthesis(const Signature& inS);
GateStringSparse LempelXSynthesis2(const Signature& inS);

//Universal T-Count finders
int UniversalTCount(SQC_Circuit* inC, int* out_daft = NULL, double* out_exec_time = NULL, int* n_parts = NULL, int* n_hadancs = NULL);

//Structured circuit generators
Signature CircuitGenerator(const string& inS);
Signature CircuitGenerator_Toffhash(int N_hash);
Signature CircuitGenerator_Toffoli(int N_toff);
Signature CircuitGenerator_RandomComplex(int n, int in_seed = 0);

//Post-synthesis Reed-Muller corrections
GateStringSparse PSRMC_Lightweight_Search(const GateStringSparse& inGSS);
GateStringSparse PSRMC_Complex_Submat(const GateStringSparse& inGSS, int n_RM=6, int col_order = 0, int n_orders = 1);
int trapezoid_width(const BMSparse& inBMS, int height);
//BMSparse decompTest(const BMSparse& inB);

//Tests
void test_PSRMC_1();
void test_Nullspace_t_vs_n();

// Global variables
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
namespace SYNTHESIS_ALGORITHM_TAG {
    const string DAFT_GUESS = "d";
    const string LEMPEL_LEAST_GREEDY = "ll";
    const string LEMPEL_GREEDY = "lg";
    const string LEMPEL_RANDOM = "lr";
    const string LEMPEL_X = "lx";
    const string REED_MULLER = "rm";
    const string ALL_LEMPEL_SELECTORS = "als";
    const string ALL_LEMPEL = "all";
}
enum SYNTHESIS_ALGORITHM_ID {
    DAFT_GUESS,
    REED_MULLER,
    LEMPEL_LEAST_GREEDY_F,
    LEMPEL_GREEDY_F,
    LEMPEL_RANDOM_F,
    LEMPEL_LEAST_GREEDY_NF,
    LEMPEL_GREEDY_NF,
    LEMPEL_RANDOM_NF,
    LEMPEL_X
};

int main(int argc, char* argv[]);

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

void NullSpaceTest() {
    //BMSparse blah = BMSparse::eye(6,6);
    BMSparse blah(6,6);

    BMSparse thisL, thisQ, thisU, thisP;


    blah.E(2,1,1);
    //blah.E(1,0,1);
    blah.E(0,2,1);
    blah.E(0,1,1);
    //blah.E(1,3,1);
    blah.E(2,3,1);
    //blah.E(1,1,1);
    blah.E(3,0,1);
    blah.E(3,1,1);
    blah.E(5,4,1);
    blah.E(5,5,1);

    cout << "Blah" << endl;
    blah.printFull();

    blah.LQUPDecomposition(thisL, thisQ, thisU, thisP);

    cout << "L" << endl;
    thisL.printFull();
    cout << "Q" << endl;
    thisQ.printFull();
    cout << "U" << endl;
    thisU.printFull();
    cout << "P" << endl;
    thisP.printFull();

    BMSparse result = blah.nullspace();
    cout << "nullspace" << endl;
    result.printFull();

    cout << "nullspace is empty? " << result.isEmpty() << endl;

    // TODO: Carry on with RowEchelon
}

void Sig2BMSTest() {
    Signature blah = Signature::SigFromFile("circuits/test3.lsf");
    blah.print();

    BMSparse blahBMS = Interface_SigBMS::SigToBMS(blah);
    cout << "A" << endl;
    blahBMS.printFull();

    BMSparse thisB = blahBMS.elementaryFactorFast();
    cout << "B" << endl;
    thisB.printFull();

    cout << "B'" << endl;
    (thisB.T()).printFull();

    BMSparse result = thisB*thisB.T();
    cout << "A = B*B'" << endl;
    result.printFull();

    BMSparse nullspace_B = thisB.nullspace();
    cout << "nullspace" << endl;
    nullspace_B.printFull();

    BMSparse minfac = blahBMS.minimalFactor();
    cout << "Minimal factor" << endl;
    minfac.printFull();

    GateStringSparse thisGSS = Interface_BMSGSS::BMSToGSS(minfac);
    thisGSS.print();

    /*
    BMSparse blee = BMSparse::eye(6,6);
    blee.E(1,0,1);
    blee.E(2,0,1);
    blee.E(2,1,1);
    blee.E(3,1,1);
    blee.E(5,1,1);
    blee.E(4,3,1);
    blee.E(5,3,1);
    blee.E(5,4,1);
    blee.E(4,5,1);
    blee.E(5,5,0);
    cout << "Factor" << endl;
    blee.printFull();
    BMSparse result2 = blee*blee.T();
    cout << "A = B*B'" << endl;
    result2.printFull();
    */
}

void compositionTest() {
    BMSparse blah = BMSparse::eye(4,4);
    BMSparse blee(1,4);
    BMSparse blan(4,1);
    BMSparse bloo = BMSparse::eye(2,2);

    BMSparse composedLR = (blan&&blah)||blah.sub(2,0,1,-1);
    BMSparse composedUD = (blee||blah)&&composedLR;
    composedUD.sub(1,3,bloo);

    composedLR.printFull();
    composedUD.printFull();

    BMSparse emptyCol(2,0);
    BMSparse emptyRow(0,2);
    BMSparse result = emptyCol*emptyRow;
    result.printFull();
}

void hdecomtest() {
    Signature blah = Signature::SigFromFile("circuits/test3.lsf");
    GateStringSparse blah_GSS = GateSigInterface::SigToGSS(blah);
    cout << "Initial gate string" << endl;
    blah.print();
    Signature dimblah = blah.diminish(1);
    cout << "Diminished gate string" << endl;
    dimblah.print();
    Signature augblah = dimblah.augment(1);
    cout << "Augmented gate string" << endl;
    augblah.print();

    Signature this_alpha, this_A, this_B;

    blah.h_decomposition(1,this_alpha, this_A, this_B);

    BMSparse this_A_BMS = Interface_SigBMS::SigToBMS(this_A);

    BMSparse minfac = this_A_BMS.minimalFactor();

    GateStringSparse this_A_GSS = Interface_BMSGSS::BMSToGSS(minfac);

    GateStringSparse result = GateSigInterface::SigToGSS(this_alpha) + this_A_GSS.mult2xh(1) + GateSigInterface::SigToGSS(this_B);
    cout << "Improved  gate string" << endl;
    result.printString();
}

GateStringSparse findNearOptimal(const Signature& inS, int maxRM) {
    cout << endl;
    cout << "Circuit to optimize:" << endl;
    inS.print();
    cout << endl;

    GateStringSparse out(inS.get_n());

    Signature currentSig(inS.get_n());
    currentSig.assign(inS);
    bool exit = false;

    cout << "Optimization begin..." << endl << endl;
    clock_t start_time = clock();
    for(int h = 1; (!exit)&&(h <= (inS.get_n()-maxRM)); h++) {
        cout << "currentSig iteration h = " << h << endl;
        Signature this_alpha, this_A, this_B;
        currentSig.h_decomposition(1,this_alpha,this_A,this_B);
        BMSparse this_A_BMS = Interface_SigBMS::SigToBMS(this_A);
        BMSparse minfac = this_A_BMS.minimalFactor();
        GateStringSparse this_A_GSS = Interface_BMSGSS::BMSToGSS(minfac);
        GateStringSparse this_alpha_GSS = GateSigInterface::SigToGSS(this_alpha);
        for(int i = 1; i < h; i++) {
            GateStringSparse tempGSS = this_A_GSS.augment(1);
            GateStringSparse tempGSS2 = this_alpha_GSS.augment(1);
            this_A_GSS.assign(tempGSS);
            this_alpha_GSS.assign(tempGSS2);
        }

        {
            GateStringSparse tempout = out + this_alpha_GSS + this_A_GSS.mult2xh(h);
            out = tempout;
        }
        currentSig.assign(this_B);

        {
            Signature tempSig = currentSig.diminish(1);
            currentSig.assign(tempSig);
        }

        exit = currentSig.isEmpty();
    }
    if(!exit) {
        //cout << "currentSig before dim" << endl;
        //currentSig.print();
        /*
        for(int h = 1; h <= (inS.get_n()-maxRM); h++) {
            Signature tempsig = currentSig.diminish(1);
            currentSig.assign(tempsig);
        }
        */
        //cout << "currentSig after dim" << endl;
        //currentSig.print();
        GateStringSparse thisRM = ReedMullerSynthesis(currentSig);
        for(int h = 1; h <= (inS.get_n()-maxRM); h++) {
            GateStringSparse tempGSS = thisRM.augment(1);
            thisRM = tempGSS;
        }
        {
            GateStringSparse tempout = out + thisRM;
            out = tempout;
        }
    }
    cout << endl;
    clock_t total_duration = clock() - start_time;
    double total_duration_s = (double)total_duration/(double)CLOCKS_PER_SEC;
    cout << "...optimization done! Time taken = " << total_duration_s << "s" << endl << endl;

    cout << "Near optimal gate string:" << endl;
    out.print();
    cout << "Near optimal weight = " << out.weight(true) << endl;
    Signature result_sig = GateSigInterface::expandGSSTerm(out);
    cout << "Near optimal signature:" << endl;
    result_sig.print();
    bool successful = (inS + result_sig).isEmpty();
    GateStringSparse daft_string = GateSigInterface::SigToGSS(inS);
    cout << "Successful? " << (successful?"Yes":"No") << endl << endl;
    cout << "Daft gate string. Weight = " << daft_string.weight(true) << "." << endl;
    daft_string.print();
    int improvement = daft_string.weight(true) - out.weight(true);
    cout << "Weight reduced by " << improvement << endl << endl;
    GateStringSparse difference = daft_string + out;
    cout << "Difference:" << endl;
    difference.print();

    return out;
}

bool T_fast(const GateStringSparse& inGSS) {
    bool out;

    int k = inGSS.get_n();
    int T = inGSS.weight(true);
    out = (T<=((double)(k*k + 3*k - 14)/2.0));

    return out;
}

int main(int argc, char* argv[]) {
    g_output_filename.clear();
    g_csv_filename.clear();

    srand(time(NULL));
    /*
    SQC_Circuit blah;
    blah.Load("test.sqc");
    blah.Print();
    blah.toffoli_n_mode = SQC_TOFFOLI_N_MODE_TOFF3;
    blah.ConvertFromToffoli();
    blah.Print();*/
    /*
    SQC_Circuit A,B;
    A.Copy(blah);
    B.Copy(blah);
    A.Clear();
    B.Clear();
    blah.GetPartition(&A,&B);
    *//*
    SQC_Circuit A;
    blah.ConvertHadamard(&A);
    A.Print();*/

    if(argc>1) {
        string this_command = argv[1];
        if(!this_command.compare("optimize")&&(argc>2)) {
            string circuit_filename = argv[2];
            int option_maxRM = 5;
            bool option_fileoutput = false;
            string option_filename;
            bool option_default_algorithm = true;
            string option_algorithm;
            bool option_feedback = true;
            int option_n_rand = 1;
            int option_psrmc = 0;
            int option_col_order = 0;
            int option_n_orders = 1;
            for(int i = 3; i < argc; i++) {
                string this_option = argv[i];
                if((this_option[0]=='-')&&((i+1)<argc)) {
                    string this_value = argv[i+1];
                    char this_option_char = this_option[1];
                    switch(this_option_char) {
                        case 'r':
                            option_maxRM = atoi(this_value.c_str());
                            break;
                        case 'o':
                            option_fileoutput = true;
                            option_filename = this_value;
                            break;
                        case 'a':
                            option_default_algorithm = false;
                            option_algorithm = this_value;
                            g_algorithm = option_algorithm;
                            break;
                        case 'f':
                            if((!this_value.compare("true"))||(!this_value.compare("1"))) {
                                option_feedback = true;
                            } else if((!this_value.compare("false"))||(!this_value.compare("0"))) {
                                option_feedback = false;
                            }
                            break;
                        case 'd':
                            option_n_rand = atoi(this_value.c_str());
                            break;
                        case 'c':
                            g_csv_filename = this_value;
                            break;
                        case 'p':
                            option_psrmc = atoi(this_value.c_str());
                            break;
                        case 'q':
                            option_col_order = atoi(this_value.c_str());
                            break;
                        case 's':
                            option_n_orders = atoi(this_value.c_str());
                            break;
                    }
                    i++;
                }
            }
            g_lempel_feedback = option_feedback;
            g_Reed_Muller_max = option_maxRM;
            g_output_filename = option_filename;
            Signature this_sig;
            char first_char = circuit_filename.c_str()[0];
            if(first_char==':') {
                cout << "Structured circuit." << endl;
                this_sig = CircuitGenerator(circuit_filename.substr(1));
            } else {
                cout << "Circuit from file: " << circuit_filename << endl;
                this_sig = Signature::SigFromFile(circuit_filename.c_str());
            }
            GateStringSparse out(this_sig.get_n());
            cout << "Input signature polynomial:" << endl;
            this_sig.print();
            cout << endl << "Synthesis algorithm: ";
            if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_LEAST_GREEDY)) {
                cout << "Lempel least-greedy" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,LempelSelector_LeastGreedy,option_feedback,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::REED_MULLER)) {
                cout << "Reed-Muller" << endl;
                GateStringSparse tempout = ReedMullerSynthesis2(this_sig);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_GREEDY)) {
                cout << "Lempel greedy" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,LempelSelector_Greedy,option_feedback,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::DAFT_GUESS)) {
                cout << "Daft guess" << endl;
                GateStringSparse tempout = GateSigInterface::SigToGSS(this_sig);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_RANDOM)) {
                cout << "Lempel random" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,LempelSelector_Random,option_feedback,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::ALL_LEMPEL_SELECTORS)) {
                cout << "All Lempel Selectors" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,option_feedback,option_n_rand,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(option_default_algorithm||!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::ALL_LEMPEL)) {
                cout << "All Lempel Selectors w/ && w/o feedback" << endl;
                GateStringSparse tempout = LempelSynthesis(this_sig,option_maxRM,option_n_rand,option_psrmc,option_col_order,option_n_orders);
                out.assign(tempout);
            } else if(!option_algorithm.compare(SYNTHESIS_ALGORITHM_TAG::LEMPEL_X)) {
                cout << "LempelX" << endl;
                GateStringSparse tempout = LempelXSynthesis(this_sig);
                out.assign(tempout);
            }
            result_analysis(this_sig,out);
            cout << "Output phase polynomial:" << endl;
            out.print();

            // PSRMC Test
            /*
            BMSparse out_daft = Interface_BMSGSS::GSSToBMS(GateSigInterface::SigToGSS(this_sig));
            BMSparse out_BMS = Interface_BMSGSS::GSSToBMS(out);
            //BMSparse myL, myE, myU, myK;
            //out_BMS.LEUKTrapezoidal(myL,myE,myU,myK,6);
            cout << "in BMS" << endl;
            out_daft.printFull();
            cout << "4x1x2x3 = " << ((out_daft.row(0).bitwise_AND(out_daft.row(1)).bitwise_AND(out_daft.row(2))).sum()%2) << endl;
            cout << "out BMS" << endl;
            out_BMS.printFull();
            cout << "4x1x2x3 = " << ((out_BMS.row(0).bitwise_AND(out_BMS.row(1)).bitwise_AND(out_BMS.row(2))).sum()%2) << endl;
            /*
            cout << "out U order = 0" << endl;
            myU.printFull();
            out_BMS.LEUKTrapezoidal(myL,myE,myU,myK,6,1);
            cout << "out U order = 1" << endl;
            myU.printFull();
            BMSparse myl, myp, myu;
            out_BMS.LPUTrapezoidal(myl,myp,myu,6);
            myu.printFull();
            */
        } else if(!this_command.compare("test")&&(argc>=3)) {
            string test_name = argv[2];
            if(!test_name.compare("PSRMC_1")) {
                test_PSRMC_1();
            } else if(!test_name.compare("Nullspace_t_vs_n")) {
                test_Nullspace_t_vs_n();
            }
        } else if(!this_command.compare("generate")&&(argc>=3)) {
            cout << "Procedurally generate a circuit." << endl;
            string circuit_name = argv[2];
            Signature this_sig = CircuitGenerator(circuit_name);
            cout << "Output signature:" << endl;
            this_sig.print();
            if(argc>=4) {
                string out_filename = argv[3];
                cout << "Saving to file: " << out_filename << "..." << endl;
                this_sig.save(out_filename.c_str());
                cout << "Saved." << endl;
            }
        } else if(!this_command.compare("list-algorithms")&&(argc>=2)) {
            cout << "ID\tCLA\t\tAlgorithm" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::DAFT_GUESS << "\t" << SYNTHESIS_ALGORITHM_TAG::DAFT_GUESS << "\t\t" << "Daft guess" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::REED_MULLER << "\t" << SYNTHESIS_ALGORITHM_TAG::REED_MULLER << "\t\t" << "Reed Muller decoder" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_LEAST_GREEDY_F << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_LEAST_GREEDY << "\t\t" << "Lempel, least greedy order, with feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_GREEDY_F << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_GREEDY << "\t\t" << "Lempel, greedy order, with feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_RANDOM_F << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_RANDOM << "\t\t" << "Lempel, random order, with feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_LEAST_GREEDY_NF << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_LEAST_GREEDY << " -f false\t" << "Lempel, least greedy order, without feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_GREEDY_NF << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_GREEDY << " -f false\t" << "Lempel, greedy order, without feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_RANDOM_NF << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_RANDOM << " -f false\t" << "Lempel, random order, without feedback" << endl;
            cout << SYNTHESIS_ALGORITHM_ID::LEMPEL_X << "\t" << SYNTHESIS_ALGORITHM_TAG::LEMPEL_X << "\t\t" << "LempelX" << endl;
        } else if(!this_command.compare("UniversalT")&&(argc>=3)) {
            cout << "Find T-Count for Universal Circuit" << endl;
            LOut_Pad++;
            string circuit_filename = argv[2];
            bool option_maslov = 0;
            SQC_ToffoliNMode option_toffmode = SQC_TOFFOLI_N_MODE_TOFF3;
            SQC_HadamardMode option_hadmode = SQC_HADAMARD_MODE_MONTANARO;
            int option_hadmode_max_anc = -1;
            for(int i = 3; i < argc; i++) {
                string this_option = argv[i];
                if((this_option[0]=='-')&&((i+1)<argc)) {
                    string this_value = argv[i+1];
                    char this_option_char = this_option[1];
                    switch(this_option_char) {
                        case 'm':
                            option_maslov = 1;//atoi(this_value.c_str());
                            break;
                        case 't':
                            if(!this_value.compare("PSCJ")) option_toffmode = SQC_TOFFOLI_N_MODE_PSCJ;
                            else if(!this_value.compare("Jones")) option_toffmode = SQC_TOFFOLI_N_MODE_JONES;
                            else if(!this_value.compare("Toff3")) option_toffmode = SQC_TOFFOLI_N_MODE_TOFF3;
                            break;
                        case 'h':
                            if(!this_value.compare("Part")) option_hadmode = SQC_HADAMARD_MODE_PARTITION;
                            else if(!this_value.compare("Mont")) option_hadmode = SQC_HADAMARD_MODE_MONTANARO;
                            break;
                        case 'a':
                            cout << "option_hadmode_max_anc = " << this_value << endl;
                            option_hadmode_max_anc = atoi(this_value.c_str());
                            break;
                    }
                    i++;
                }
            }
            SQC_Circuit this_circuit;
            if(option_maslov) {
                this_circuit.LoadMaslovFile(circuit_filename.c_str());
            } else {
                this_circuit.Load(circuit_filename.c_str());
            }
            this_circuit.Print();
            double exec_time = 0;
            int daft_t_count = 0;
            int n_parts;
            int n_hadancs = 0;
            this_circuit.toffoli_n_mode = option_toffmode;
            this_circuit.hadamard_mode = option_hadmode;
            this_circuit.hadamard_mode_max_ancillas = option_hadmode_max_anc;
            int result = UniversalTCount(&this_circuit, &daft_t_count, &exec_time,&n_parts,&n_hadancs);
            cout << "Final T-count (LempelX) = " << result << endl;
            cout << "Final T-count (Daft) = " << daft_t_count << endl;
            cout << "Final T-count reduced by " << (daft_t_count-result) << endl;
            cout << "Total execution time = " << exec_time << "s" << endl;
            cout << "No. partitions = " << n_parts << endl;
            cout << "No. hadamard ancillas = " << n_hadancs << endl;
            this_circuit.Destruct();
        }
    }
    return 0;
}

int LempelSelector_LeastGreedy(const Signature& inS) {
    // Not yet tested.
    int out = 0;

    int n = inS.get_n();
    int* h_weights = new int[n];
    int* h_order = new int[n];

    for(int i = 0; i < n; i++) {
        h_weights[i] = inS.h_weight(i+1);
    }

    LukeInt::sort(h_weights,n,false,h_order);

    //LukeInt::print(h_weights,n,"h_weights: ");
    //LukeInt::print(h_order,n,"h_order: ");

    bool exit = false;
    for(int i = 0; (!exit)&&(i < n); i++) {
        if(h_weights[i]) {
            out = h_order[i];
            exit = true;
        }
    }

    delete [] h_order;
    delete [] h_weights;

    return (out+1);// Returns qubit number - not array index.
}

int LempelSelector_Greedy(const Signature& inS) {
    // Not yet tested.
    int out = 0;

    int n = inS.get_n();
    int* h_weights = new int[n];
    int* h_order = new int[n];

    for(int i = 0; i < n; i++) {
        h_weights[i] = inS.h_weight(i+1);
    }

    LukeInt::sort(h_weights,n,true,h_order);

    out = h_order[0];

    delete [] h_order;
    delete [] h_weights;

    return (out+1);// Returns qubit number - not array index.
}

int LempelSelector_Random(const Signature& inS) {
    // Not yet tested.
    int out = 0;

    int n = inS.get_n();
    int* h_weights = new int[n];
    int* h_order = new int[n];

    for(int i = 0; i < n; i++) {
        h_weights[i] = inS.h_weight(i+1);
    }

    LukeInt::sort(h_weights,n,true,h_order);

    int n_eff = inS.get_N_eff();

    out = h_order[LukeInt::randi(0,n_eff-1)];

    delete [] h_order;
    delete [] h_weights;

    return (out+1);// Returns qubit number - not array index.
}

GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, LempelSelector lempelSelector, bool feedback, int psrmc, int col_order, int n_orders) {
    cout << "Gate synthesis begin." << endl;
    int n = inS.get_n();
    GateStringSparse out(inS.get_n());

    Signature thisSig;
    thisSig.assign(inS);

    int n_eff = thisSig.get_N_eff();

    cout << "Lempel factoring begin." << endl;
    clock_t start_tic = clock();
    g_h_order_stream.str("");
    while(n_eff>maxRM) {
        int h = lempelSelector(thisSig);
        string stream_contents = g_h_order_stream.str();
        if(!stream_contents.empty()) g_h_order_stream << ", ";
        g_h_order_stream << h;
        clock_t toc = clock();
        double perc_done = 100.0*(double)(n-n_eff)/(n-maxRM);
        double time_remaining = (100.0-perc_done)*perc_done/((double)(toc-start_tic)/CLOCKS_PER_SEC);
        cout << "n_eff = " << n_eff << "; h = " << h << "; ";
        cout << perc_done << "%; " << time_remaining << "s" << endl;
        //cout << "thisSig" << endl;
        //thisSig.print();
        //cout << "outSig" << endl;
        //GateSigInterface::expandGSSTerm(out).print();
        Signature this_alpha, this_Ap, this_B;
        thisSig.h_decomposition(h,this_alpha,this_Ap,this_B);
        //cout << "this_alpha" << endl; this_alpha.print();
        //cout << "this_Ap" << endl; this_Ap.print();
        //cout << "this_B" << endl; this_B.print();
        BMSparse this_Ap_BMS = Interface_SigBMS::SigToBMS(this_Ap);
        BMSparse minfac = this_Ap_BMS.minimalFactor();
        GateStringSparse this_Ap_GSS = Interface_BMSGSS::BMSToGSS(minfac);
        GateStringSparse this_out(n);
        if(feedback) {
            if(!this_Ap.isEmpty()) {
                if(this_Ap_GSS.weight(true)%2) {
                    bool* this_xh = new bool[n];
                    for(int i = 0; i < n; i++) this_xh[i] = (i==(h-1));
                    this_out.set(this_xh);
                    delete [] this_xh;
                }
                this_out = (this_out + this_Ap_GSS.addxh(h) + GateSigInterface::SigToGSS(this_alpha));
            } else {
                this_out = GateSigInterface::SigToGSS(this_alpha);
            }
            thisSig = (GateSigInterface::expandGSSTerm(this_Ap_GSS) + this_B);
        } else {
            this_out = (this_Ap_GSS.mult2xh(h) + GateSigInterface::SigToGSS(this_alpha));
            thisSig = this_B;
        }
        out = (out + this_out);
        n_eff = thisSig.get_N_eff();
    }
    cout << "Lempel factoring end." << endl << endl;

    bool* rm_inc_vec = new bool[n];
    for(int i = 0; i < n; i++) rm_inc_vec[i] = (bool)thisSig.h_weight(i+1);

    Signature rm_sig = thisSig.subsig(rm_inc_vec);

    cout << "Reed-Muller synthesis begin." << endl << endl;
    GateStringSparse rm_GSS = ReedMullerSynthesis(rm_sig);
    cout << "Reed-Muller synthesis end." << endl << endl;

    GateStringSparse restored_rm_GSS = rm_GSS.superGSS(rm_inc_vec,n);

    out = (out + restored_rm_GSS);

    cout << "Gate synthesis end." << endl << endl;

    if(psrmc) {
        cout << "Post-Synthesis Reed-Muller Corrections begin." << endl << endl;
        int pre_psrmc_weight = out.weight(true);
        int post_psrmc_weight;
        switch(psrmc) {
            case 1: // Lightweight search
            {
                cout << "PSRM Algorithm: Lightweight search." << endl << endl;
            }
            break;
            case 2: // Complex gate-synthesis submatrix finding v1
            {
                cout << "PSRM Algorithm: Complex gate-synthesis submatrix. v1" << endl << endl;
                GateStringSparse this_complex_GSS = PSRMC_Complex_Submat(out,maxRM,col_order,n_orders);
                Signature this_complex_sig = GateSigInterface::expandGSSTerm(this_complex_GSS);
                bool* rm_inc_vec2 = new bool[n];
                for(int i = 0; i < n; i++) rm_inc_vec2[i] = (bool)this_complex_sig.h_weight(i+1);

                Signature rm_sig2 = this_complex_sig.subsig(rm_inc_vec2);

                cout << "Reed-Muller synthesis begin." << endl << endl;
                GateStringSparse rm_GSS2 = ReedMullerSynthesis(rm_sig2);
                cout << "Reed-Muller synthesis end." << endl << endl;

                GateStringSparse restored_rm_GSS2 = rm_GSS2.superGSS(rm_inc_vec2,n);

                cout << "GSS before:" << endl;
                this_complex_GSS.print();
                cout << "GSS after:" << endl;
                restored_rm_GSS2.print();
                out = (out + this_complex_GSS + restored_rm_GSS2);
                post_psrmc_weight = out.weight(true);
            }
            break;
            case 3: // Complex sub matrix finding col > rank
            {
                cout << "PSRMC Algorithm: Complex sub-matrix col(U)>rank(U)" << endl << endl;
                BMSparse out_BMS = Interface_BMSGSS::GSSToBMS(out);
                BMSparse thisE, thisU, thisK;
                int this_width = out_BMS.EUKPentagonal(thisE, thisU, thisK, maxRM);

                cout << "E" << endl;
                thisE.printFull();
                cout << "U" << endl;
                thisU.printFull();
                cout << "K" << endl;
                thisK.printFull();

                BMSparse rm_in_BMS = thisU.sub(0,0,maxRM,this_width);
                GateStringSparse rm_in_GSS = Interface_BMSGSS::BMSToGSS(rm_in_BMS);
                Signature rm_in_sig = GateSigInterface::expandGSSTerm(rm_in_GSS);

                cout << "RM Sub-matrix" << endl;
                rm_in_BMS.printFull();
                cout << "RM GSS" << endl;
                rm_in_GSS.print();
                cout << "RM Sig" << endl;
                rm_in_sig.print();

                cout << "Reed-Muller synthesis begin." << endl << endl;
                GateStringSparse rm_out_GSS = ReedMullerSynthesis(rm_in_sig);
                cout << "Reed-Muller synthesis end." << endl << endl;

                BMSparse rm_out_BMS = Interface_BMSGSS::GSSToBMS(rm_out_GSS);
                BMSparse rm_out_BMS_super(n,rm_out_BMS.get_m());
                rm_out_BMS_super.sub(0,0,rm_out_BMS);
                BMSparse rm_out_BMS_full = rm_out_BMS_super&&thisU.sub(0,this_width,-1,-1);
                {
                    BMSparse temp = thisE*rm_out_BMS_full;
                    rm_out_BMS_full.assign(temp);
                }
                GateStringSparse rm_out_GSS_super = Interface_BMSGSS::BMSToGSS(rm_out_BMS_full);

                cout << "GSS before:" << endl;
                out.print();
                cout << "GSS after:" << endl;
                rm_out_BMS_super.print();
                out.assign(rm_out_GSS_super);
                post_psrmc_weight = out.weight(true);
            }
            break;
        }
        cout << "PSRMC reduced weight by " << (pre_psrmc_weight-post_psrmc_weight) << endl;
        cout << "Post-Synthesis Reed-Muller Corrections end." << endl << endl;
    }

    delete [] rm_inc_vec;
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

GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, bool feedback, int n_rand, int psrmc, int col_order, int n_orders) {
    GateStringSparse out(inS.get_n());

    int n_algorithms = 4;
    cout << "Algorithm: Daft" << endl;
    GateStringSparse result_daft = GateSigInterface::SigToGSS(inS);
    cout << "Algorithm: Lempel least-greedy" << endl;
    GateStringSparse result_Lempel_ll = LempelSynthesis(inS, maxRM, LempelSelector_LeastGreedy,feedback);
    cout << "Algorithm: Lempel greedy" << endl;
    GateStringSparse result_Lempel_lg = LempelSynthesis(inS, maxRM,LempelSelector_Greedy,feedback);
    cout << "Algorithm: Lempel random" << endl;
    GateStringSparse result_Lempel_lr(inS.get_n());

    int rand_weight = -1;
    for(int i = 0; i < n_rand; i++) {
        cout << "Lempel random iteration " << (i+1) << " of " << n_rand << endl;
        GateStringSparse this_Lempel_lr = LempelSynthesis(inS,maxRM,LempelSelector_Random,feedback);
        int this_rand_weight = this_Lempel_lr.weight(true);
        if((rand_weight == -1)||(this_rand_weight<rand_weight)) {
            rand_weight = this_rand_weight;
            result_Lempel_lr = this_Lempel_lr;
        }
    }

    int weights[n_algorithms], weight_order[n_algorithms];
    weights[0] = result_daft.weight(true);
    weights[1] = result_Lempel_ll.weight(true);
    weights[2] = result_Lempel_lg.weight(true);
    weights[3] = result_Lempel_lr.weight(true);

    cout << endl << "T-counts:" << endl;
    cout << "Daft = " << weights[0] << endl;
    cout << "Lempel least-greedy = " << weights[1] << endl;
    cout << "Lempel greedy = " << weights[2] << endl;
    cout << "Lempel random = " << weights[3] << endl << endl;

    LukeInt::sort(weights, n_algorithms, false, weight_order);

    cout << "Best solution: ";
    switch(weight_order[0]) {
        case 0:
            cout << "Daft" << endl;
            out = result_daft;
            break;
        case 1:
            cout << "Lempel least-greedy" << endl;
            out = result_Lempel_ll;
            break;
        case 2:
            cout << "Lempel greedy" << endl;
            out = result_Lempel_lg;
            break;
        case 3:
            cout << "Lempel random" << endl;
            out = result_Lempel_lr;
            break;
    }

    return out;
}

GateStringSparse LempelSynthesis(const Signature& inS, int maxRM, int n_rand, int psrmc, int col_order, int n_orders) {
    GateStringSparse out(inS.get_n());

    double time_d;
    double time_ll_f, time_lg_f, time_lr_f;
    double time_ll_nf, time_lg_nf, time_lr_nf;
    clock_t tic, toc;

    int n_algorithms = 7;
    cout << "Algorithm: Daft" << endl;
        tic = clock();
        GateStringSparse result_daft = GateSigInterface::SigToGSS(inS);
        toc = clock();
        time_d = (double)(toc-tic)/CLOCKS_PER_SEC;
    cout << "Algorithm: Lempel least-greedy; Feedback: True" << endl;
        tic = clock();
        GateStringSparse result_Lempel_ll_f = LempelSynthesis(inS, maxRM, LempelSelector_LeastGreedy,true);
        toc = clock();
        time_ll_f = (double)(toc-tic)/CLOCKS_PER_SEC;
    cout << "Algorithm: Lempel greedy; Feedback: True" << endl;
        tic = clock();
        GateStringSparse result_Lempel_lg_f = LempelSynthesis(inS, maxRM,LempelSelector_Greedy,true);
        toc = clock();
        time_lg_f = (double)(toc-tic)/CLOCKS_PER_SEC;
    cout << "Algorithm: Lempel random; Feedback: True" << endl;
        tic = clock();
        GateStringSparse result_Lempel_lr_f(inS.get_n());
        {
            int rand_weight = -1;
            for(int i = 0; i < n_rand; i++) {
                cout << "Lempel random iteration " << (i+1) << " of " << n_rand << endl;
                GateStringSparse this_Lempel_lr = LempelSynthesis(inS,maxRM,LempelSelector_Random,true);
                int this_rand_weight = this_Lempel_lr.weight(true);
                if((rand_weight == -1)||(this_rand_weight<rand_weight)) {
                    rand_weight = this_rand_weight;
                    result_Lempel_lr_f = this_Lempel_lr;
                    g_best_random_h_order_f = g_h_order_stream.str();
                }
            }
        }
        toc = clock();
        time_lr_f = (double)(toc-tic)/CLOCKS_PER_SEC;

    cout << "Algorithm: Lempel least-greedy; Feedback: False" << endl;
        tic = clock();
        GateStringSparse result_Lempel_ll_nf = LempelSynthesis(inS, maxRM, LempelSelector_LeastGreedy,false);
        toc = clock();
        time_ll_nf = (double)(toc-tic)/CLOCKS_PER_SEC;
    cout << "Algorithm: Lempel greedy; Feedback: False" << endl;
        tic = clock();
        GateStringSparse result_Lempel_lg_nf = LempelSynthesis(inS, maxRM,LempelSelector_Greedy,false);
        toc = clock();
        time_lg_nf = (double)(toc-tic)/CLOCKS_PER_SEC;
    cout << "Algorithm: Lempel random; Feedback: False" << endl;
        tic = clock();
        GateStringSparse result_Lempel_lr_nf(inS.get_n());
        {
            int rand_weight = -1;
            for(int i = 0; i < n_rand; i++) {
                cout << "Lempel random iteration " << (i+1) << " of " << n_rand << endl;
                GateStringSparse this_Lempel_lr = LempelSynthesis(inS,maxRM,LempelSelector_Random,false);
                int this_rand_weight = this_Lempel_lr.weight(true);
                if((rand_weight == -1)||(this_rand_weight<rand_weight)) {
                    rand_weight = this_rand_weight;
                    result_Lempel_lr_nf = this_Lempel_lr;
                    g_best_random_h_order_nf = g_h_order_stream.str();
                }
            }
        }
        toc = clock();
        time_lr_nf = (double)(toc-tic)/CLOCKS_PER_SEC;

    int weights[n_algorithms], weight_order[n_algorithms];
    weights[0] = result_daft.weight(true);
    weights[1] = result_Lempel_ll_f.weight(true);
    weights[2] = result_Lempel_lg_f.weight(true);
    weights[3] = result_Lempel_lr_f.weight(true);
    weights[4] = result_Lempel_ll_nf.weight(true);
    weights[5] = result_Lempel_lg_nf.weight(true);
    weights[6] = result_Lempel_lr_nf.weight(true);

    cout << "Execution times:" << endl;
    cout << "Daft: " << time_d << "s" << endl;
    cout << "Lempel ll f: " << time_ll_f << "s" << endl;
    cout << "Lempel lg f: " << time_lg_f << "s" << endl;
    cout << "Lempel lr f: " << time_lr_f << "s" << endl;
    cout << "Lempel ll nf: " << time_ll_nf << "s" << endl;
    cout << "Lempel lg nf: " << time_lg_nf << "s" << endl;
    cout << "Lempel lr nf: " << time_lr_nf << "s" << endl << endl;

    cout << "Success:" << endl;
    cout << "Daft: " << synthesis_success(inS,result_daft) << endl;
    cout << "Lempel ll f: " << synthesis_success(inS,result_Lempel_ll_f) << endl;
    cout << "Lempel lg f: " << synthesis_success(inS,result_Lempel_lg_f) << endl;
    cout << "Lempel lr f: " << synthesis_success(inS,result_Lempel_lr_f) << endl;
    cout << "Lempel ll nf: " << synthesis_success(inS,result_Lempel_ll_nf) << endl;
    cout << "Lempel lg nf: " << synthesis_success(inS,result_Lempel_lg_nf) << endl;
    cout << "Lempel lr nf: " << synthesis_success(inS,result_Lempel_lr_nf) << endl;

    cout << endl << "T-counts:" << endl << endl;
    cout << "Daft = " << weights[0] << endl << endl;
    cout << "Lempel | feedback = True" << endl;
    cout << "Lempel least-greedy = " << weights[1] << endl;
    cout << "Lempel greedy = " << weights[2] << endl;
    cout << "Lempel random (d = " << n_rand << ") = " << weights[3] << endl;
    cout << "\tBest random h-order: " << g_best_random_h_order_f << endl << endl;
    cout << "Lempel | feedback = False" << endl;
    cout << "Lempel least-greedy = " << weights[4] << endl;
    cout << "Lempel greedy = " << weights[5] << endl;
    cout << "Lempel random (d = " << n_rand << ") = " << weights[6] << endl;
    cout << "\tBest random h-order: " << g_best_random_h_order_nf << endl << endl;

    LukeInt::sort(weights, n_algorithms, false, weight_order);

    cout << "Best solution: ";
    switch(weight_order[0]) {
        case 0:
            cout << "Daft" << endl;
            out = result_daft;
            break;
        case 1:
            cout << "Lempel least-greedy w/ feedback" << endl;
            out = result_Lempel_ll_f;
            break;
        case 2:
            cout << "Lempel greedy w/ feedback" << endl;
            out = result_Lempel_lg_f;
            break;
        case 3:
            cout << "Lempel random w/ feedback" << endl;
            out = result_Lempel_lr_f;
            break;
        case 4:
            cout << "Lempel least-greedy w/o feedback" << endl;
            out = result_Lempel_ll_nf;
            break;
        case 5:
            cout << "Lempel greedy w/o feedback" << endl;
            out = result_Lempel_lg_nf;
            break;
        case 6:
            cout << "Lempel random w/o feedback" << endl;
            out = result_Lempel_lr_nf;
            break;
    }

    if(!g_output_filename.empty()) {
        ofstream myfile(g_output_filename.c_str(), ostream::app);
        myfile << "Input signature polynomial:" << endl;
        inS.print(myfile);
        myfile << "Output phase polynomial:" << endl;
        out.print(myfile);
        myfile << "Execution times:" << endl;
        myfile << "Daft: " << time_d << "s" << endl;
        myfile << "Lempel ll f: " << time_ll_f << "s" << endl;
        myfile << "Lempel lg f: " << time_lg_f << "s" << endl;
        myfile << "Lempel lr f: " << time_lr_f << "s" << endl;
        myfile << "Lempel ll nf: " << time_ll_nf << "s" << endl;
        myfile << "Lempel lg nf: " << time_lg_nf << "s" << endl;
        myfile << "Lempel lr nf: " << time_lr_nf << "s" << endl << endl;
        myfile << endl << "Success:" << endl;
        myfile << "Daft: " << synthesis_success(inS,result_daft) << endl;
        myfile << "Lempel ll f: " << synthesis_success(inS,result_Lempel_ll_f) << endl;
        myfile << "Lempel lg f: " << synthesis_success(inS,result_Lempel_lg_f) << endl;
        myfile << "Lempel lr f: " << synthesis_success(inS,result_Lempel_lr_f) << endl;
        myfile << "Lempel ll nf: " << synthesis_success(inS,result_Lempel_ll_nf) << endl;
        myfile << "Lempel lg nf: " << synthesis_success(inS,result_Lempel_lg_nf) << endl;
        myfile << "Lempel lr nf: " << synthesis_success(inS,result_Lempel_lr_nf) << endl;
        myfile << endl << "T-counts:" << endl << endl;
        myfile << "Daft = " << result_daft.weight(true) << endl << endl;
        myfile << "Lempel | feedback = True" << endl;
        myfile << "Lempel least-greedy = " << result_Lempel_ll_f.weight(true) << endl;
        myfile << "Lempel greedy = " << result_Lempel_lg_f.weight(true) << endl;
        myfile << "Lempel random (d = " << n_rand << ") = " << result_Lempel_lr_f.weight(true) << endl;
        myfile << "\tBest random h-order: " << g_best_random_h_order_f << endl << endl;
        myfile << "Lempel | feedback = False" << endl;
        myfile << "Lempel least-greedy = " << result_Lempel_ll_nf.weight(true) << endl;
        myfile << "Lempel greedy = " << result_Lempel_lg_nf.weight(true) << endl;
        myfile << "Lempel random (d = " << n_rand << ") = " << result_Lempel_lr_nf.weight(true) << endl;
        myfile << "\tBest random h-order: " << g_best_random_h_order_nf << endl << endl;
        myfile << "Best algorithm: ";
        switch(weight_order[0]) {
            case 0:
                myfile << "Daft" << endl;
                out = result_daft;
                break;
            case 1:
                myfile << "Lempel least-greedy w/ feedback" << endl;
                out = result_Lempel_ll_f;
                break;
            case 2:
                myfile << "Lempel greedy w/ feedback" << endl;
                out = result_Lempel_lg_f;
                break;
            case 3:
                myfile << "Lempel random w/ feedback" << endl;
                out = result_Lempel_lr_f;
                break;
            case 4:
                myfile << "Lempel least-greedy w/o feedback" << endl;
                out = result_Lempel_ll_nf;
                break;
            case 5:
                myfile << "Lempel greedy w/o feedback" << endl;
                out = result_Lempel_lg_nf;
                break;
            case 6:
                myfile << "Lempel random w/o feedback" << endl;
                out = result_Lempel_lr_nf;
                break;
        }
        result_analysis(inS,out,myfile);
        myfile.close();
    }

    if(!g_csv_filename.empty()) {
        if(g_random_circuit_seed==0) {
            ofstream my_csv(g_csv_filename.c_str(),iostream::app);

            my_csv << g_indvar_out << ",";
            my_csv << result_daft.weight(true) << ",";
            my_csv << result_Lempel_ll_f.weight(true) << ",";
            my_csv << result_Lempel_lg_f.weight(true) << ",";
            my_csv << result_Lempel_lr_f.weight(true) << ",";
            my_csv << result_Lempel_ll_nf.weight(true) << ",";
            my_csv << result_Lempel_lg_nf.weight(true) << ",";
            my_csv << result_Lempel_lr_nf.weight(true) << ",";
            my_csv << n_rand << endl;

            my_csv.close();
        } else {
            ofstream my_csv(g_csv_filename.c_str(),iostream::app);
            int n = inS.get_n();
            my_csv << SYNTHESIS_ALGORITHM_ID::DAFT_GUESS << "," << n << "," << g_random_circuit_seed << "," << result_daft.weight(true) << endl;
            my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_LEAST_GREEDY_F << "," << n << "," << g_random_circuit_seed << "," << result_Lempel_ll_f.weight(true) << endl;
            my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_GREEDY_F << "," << n << "," << g_random_circuit_seed << "," << result_Lempel_lg_f.weight(true) << endl;
            my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_RANDOM_F << "," << n << "," << g_random_circuit_seed << "," << result_Lempel_lr_f.weight(true) << endl;
            my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_LEAST_GREEDY_NF << "," << n << "," << g_random_circuit_seed << "," << result_Lempel_ll_nf.weight(true) << endl;
            my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_GREEDY_NF << "," << n << "," << g_random_circuit_seed << "," << result_Lempel_lg_nf.weight(true) << endl;
            my_csv << SYNTHESIS_ALGORITHM_ID::LEMPEL_RANDOM_NF << "," << n << "," << g_random_circuit_seed << "," << result_Lempel_lr_nf.weight(true) << endl;
            my_csv.close();
        }
    }

    //PSRMC_Lightweight_Search(out);

    return out;
}

Signature CircuitGenerator_Toffhash(int N_hash) {
    int n = 2*N_hash+1;
    Signature out(n);
    for(int t = 0; t < N_hash; t++){
        out.set(1,2*t+2,2*t+3);
    }
    ostringstream temp_ss;
    temp_ss << N_hash;
    g_indvar_out = temp_ss.str();
    return out;
}

Signature CircuitGenerator_Toffoli(int N_toff) {
    int n = 3*N_toff;
    Signature out(n);
    for(int t = 0; t < N_toff; t++){
        out.set(3*t+1,3*t+2,3*t+3);
    }
    ostringstream temp_ss;
    temp_ss << N_toff;
    g_indvar_out = temp_ss.str();
    return out;
}

Signature CircuitGenerator(const string& inS) {
    Signature out;

    char split_char = '-';
    int split_pos = inS.find(split_char);
    string circuit_name = inS.substr(0,split_pos);
    string circuit_args = inS.substr(split_pos+1);

    if(!circuit_name.compare("Toffhash")) {
        int N_hash = atoi(circuit_args.c_str());
        cout << "Toffhash circuit; N_hash = " << N_hash << endl;
        out = CircuitGenerator_Toffhash(N_hash);
    } else if(!circuit_name.compare("Toffoli")) {
        int N_toff = atoi(circuit_args.c_str());
        cout << "Toffoli circuit; N_toff = " << N_toff << endl;
        out = CircuitGenerator_Toffoli(N_toff);
    } else if(!circuit_name.compare("RandomComplex")) {
        int n;
        int this_seed = 0;
        int arg_split_pos = circuit_args.find(split_char);
        if(arg_split_pos!=circuit_args.npos) {
            string seed_string = circuit_args.substr(arg_split_pos+1);
            this_seed = atoi(seed_string.c_str());
        }
        n = atoi(circuit_args.substr(0,arg_split_pos).c_str());

        cout << "Random complex circuit; n = " << n;
        if(this_seed>0) {
            cout << "; seed = " << this_seed;
            g_random_circuit_seed = this_seed;
        }
        cout << endl;
        out = CircuitGenerator_RandomComplex(n,this_seed);
    }

    return out;
}

Signature CircuitGenerator_RandomComplex(int n, int in_seed) {
    Signature out;

    GateStringSparse complex_GSS = GateStringSparse::randomUpTo3Qu(n, in_seed);
    out = GateSigInterface::expandGSSTerm(complex_GSS);

    ostringstream temp_ss;
    temp_ss << n;
    g_indvar_out = temp_ss.str();

    return out;
}

bool synthesis_success(const Signature& inS, const GateStringSparse& inResult) {
    bool out;

    //Signature restored_sig = GateSigInterface::expandGSSTerm(inResult);
    Signature restored_sig = Interface_SigBMS::BMSToSig(Interface_BMSGSS::GSSToBMS(inResult));
    Signature difference = inS + restored_sig;
    out = difference.isEmpty();

    return out;
}

GateStringSparse PSRMC_Lightweight_Search(const GateStringSparse& inGSS) {
    int n = inGSS.get_n();
    GateStringSparse out(n);
    GateStringSparse best_GSS(n), this_GSS(n);

    bool* this_s = new bool[n];
    bool* best_s = new bool[n];
    bool* temp_s = new bool[n];
    for(int i = 0; i < n; i++) this_s[i] = (i < (n-5));
    int this_weight;
    int best_weight = -1;
    bool exit = false;
    int it = 0;
    do {
        this_GSS = inGSS.subGSS_s(this_s);
        this_weight = this_GSS.weight(true);
        if((this_weight>best_weight)||(best_weight==-1)) {
            best_weight = this_weight;
            best_GSS = this_GSS;
            LukeBool::copy(this_s,best_s,n);
        }
        exit = !LukeBool::nextUniquePerm(temp_s,this_s,n);
        LukeBool::copy(temp_s,this_s,n);
        it++;
    } while(!exit);

    if(best_weight>=8) {
        cout << "PSRMC Algorithm can reduce T-count." << endl;
    } else {
        cout << "PSRMC Algorithm can't reduce T-count." << endl;
    }
    cout << "Best weight = " << best_weight << endl;
    LukeBool::print(best_s,n,"Best s: ");
    cout << "Number of iterations = " << it << endl;

    delete [] this_s;
    delete [] best_s;
    delete [] temp_s;

    return out;
}

GateStringSparse PSRMC_Complex_Submat(const GateStringSparse& inGSS, int n_RM, int col_order, int n_orders) {
    int n = inGSS.get_n();
    GateStringSparse out(n);

    BMSparse inBMS = Interface_BMSGSS::GSSToBMS(inGSS);
    BMSparse this_C = inBMS.sub_weight("<=",n_RM,2);

    BMSparse this_M(n,0);

    int* rand_col_order = NULL;
    if(col_order==1) rand_col_order = new int[this_C.get_m()];

    for(int j = 0; j < this_C.get_m(); j++) {
        for(int order = 0; order < n_orders; order++) {
            BMSparse this_M_j = this_C.col(j);
            BMSparse this_d = this_C.col(j);
            if(col_order==1) LukeInt::randperm(rand_col_order, this_C.get_m());
            for(int k = 0; k < this_C.get_m(); k++) { // Result depends on order of k - implement randomisation of order
                int k_col;
                switch(col_order) { // Order of remaining columns
                    case 0: // Numeric order
                        k_col = k;
                        break;
                    case 1: // Random order
                        k_col = rand_col_order[k];
                        break;

                }
                if(k_col!=j) {
                    BMSparse this_c_k = this_C.col(k_col);
                    this_d = this_d.bitwise_OR(this_c_k);
                    if(this_d.sum()<=n_RM) { // Change to column-rank <= 6
                        this_M_j = (this_M_j&&this_c_k);
                    }
                }
            }
            if(this_M_j.get_m()>this_M.get_m()) { // Try different conditions e.g. if(sum(Mj)>sum(M)). Also make sure not same qubits as first RM pass.
                this_M = this_M_j;
            }
        }

    }

    if(col_order==1) delete [] rand_col_order;

    out = Interface_BMSGSS::BMSToGSS(this_M);

    return out;
}

/*
BMSparse decompTest(const BMSparse& inB) {
    BMSparse out;



    return out;
}*/

void test_PSRMC_1() {
    int n = 10;
    int n_RM = 5;
    Signature mySig = GateSigInterface::expandGSSTerm(GateStringSparse::randomUpTo3Qu(n));
    cout << "Input signature" << endl;
    mySig.print();

    GateStringSparse result = LempelSynthesis(mySig,n_RM,LempelSelector_Random,true);
    cout << "Output phase polynomial" << endl;
    result.print();
    result.printString();

    BMSparse result_BMS = Interface_BMSGSS::GSSToBMS(result);
    cout << "Synthesis matrix" << endl;
    result_BMS.printFull();

    {
        BMSparse myL, myE, myU, myK;
        result_BMS.LEUKTrapezoidal(myL,myE,myU,myK,n_RM,1);
        cout << "leUk random" << endl;
        myU.printFull();
        cout << "width = " << trapezoid_width(myU,n_RM) << endl;
        cout << "Nullspace" << endl;
        BMSparse myNull = myU.nullspace();
        myNull.printFull();
        (myU*myNull.col(0)).printFull();
    }
    {
        BMSparse myL, myE, myU, myK;
        result_BMS.LEUKTrapezoidal(myL,myE,myU,myK,n_RM,0);
        cout << "leUk many 1's" << endl;
        myU.printFull();
        cout << "width = " << trapezoid_width(myU,n_RM) << endl;
        cout << "Nullspace" << endl;
        BMSparse myNull = myU.nullspace();
        myNull.printFull();
        (myU*myNull.col(0)).printFull();
    }

}

int trapezoid_width(const BMSparse& inBMS, int height) {
    int out = 0;

    int end_zeros;
    bool exit = false;
    do {
        end_zeros = 0;
        for(int i = height; i < inBMS.get_n(); i++) {
            if(!inBMS.E(i,out)) end_zeros++;
        }
        if(end_zeros>=(inBMS.get_n()-height)) {
            out++;
        } else {
            exit = true;
        }
    } while(!exit);


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

void test_Nullspace_t_vs_n() {
    string out_filename = "./test_Nullspace_t_vs_n.csv";
    int start_n = 6;
    int end_n = 30;
    double max_time = 5*60.0;
    int max_d = 1000;
    double P = 0.5;
    ofstream myfile(out_filename.c_str(),iostream::app);
    cout << "Test: Nullspace_t_vs_n begin." << endl;
    for(int n = start_n; n <= end_n; n++) {
        cout << "n = " << n << endl;
        double total_exec_time = 0.0;
        int d = 0;
        while((total_exec_time<max_time)&&(d<max_d)) {
            cout << "Iteration d = " << d << " begin." << endl;
            BMSparse this_BMS = BMSparse::random(n,n,P);
            clock_t tic = clock();
            BMSparse this_nullspace = this_BMS.nullspace();
            clock_t toc = clock();
            double exec_time = (double)(toc-tic)/CLOCKS_PER_SEC;
            cout << "Executed in " << exec_time << "s" << endl;
            total_exec_time += exec_time;
            d++;
        }
        double mean_exec_time = total_exec_time/d;
        // Line format for csv file: <n>,<t>,<N_samples>
        myfile << n << "," << mean_exec_time << "," << d << endl;
    }
    cout << "Test: Nullspace_t_vs_n end." << endl;
    myfile.close();
}

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

GateStringSparse LempelXSynthesis(const Signature& inS) {
    LOut(); cout << "Gate synthesis begin." << endl;
    GateStringSparse A_GSS = GateSigInterface::SigToGSS(inS);
    BMSparse A_BMS = Interface_BMSGSS::GSSToBMS(A_GSS);
    int n = A_BMS.get_n();
    int m = A_BMS.get_m();
    bool** A_bool = LukeMat_GF2::construct(n,m);
    A_BMS.toBool(A_bool);
    int t;
    clock_t tic = clock();
    LOut_Pad++;
    GateSynthesisMatrix::LempelX(A_bool,n,m,t);
    clock_t toc = clock();
    double exec_time = ((double)toc-(double)tic)/CLOCKS_PER_SEC;
    BMSparse out_BMS(n,t);
    out_BMS.fromBool(A_bool,n,t);
    GateStringSparse out = Interface_BMSGSS::BMSToGSS(out_BMS);
    LukeMat_GF2::destruct(A_bool,n,m);
    LOut(); cout << "LempelX executed in " << exec_time << " seconds." << endl;
    LOut_Pad--;
    LOut(); cout << "Gate synthesis end." << endl;
    return out;
}

int UniversalTCount(SQC_Circuit* inC, int* out_daft, double* out_exec_time, int* n_parts, int* n_hadancs) {
    //cout << "C" << endl;
    int out = 0;
    if(out_daft) *out_daft = 0;
    //inC->Print();
    //cout << endl;
    int tic = clock();
    bool exit = true;
    Signature this_sig;
    int this_part = 1;
    LOut() << "UniversalTCount: Algorithm begin." << endl;
    LOut_Pad++;
    inC->ConvertFromToffoli();
    LOut() << "UniversalTCount: Converted to elementary gate set." << endl;
    //inC->Print(&cout,0,65);
    LOut() << endl;
    //return out;

    do {
        LOut() << "Partition " << this_part << endl;
        //inC->Print();
        LOut_Pad++;
        int init_gates = inC->m;
        LOut() << "Remaining gates = " << init_gates << endl;
        exit = !inC->NextSignature(this_sig);
        int reduces_gates = init_gates-inC->m;
        LOut() << "Converted" << endl;
        //inC->Print();
        this_sig.print();
        GateStringSparse result = LempelXSynthesis(this_sig);
        GateStringSparse result_daft = GateSigInterface::SigToGSS(this_sig);
        int this_t_count = result.weight(true);
        int this_t_count_daft = result_daft.weight(true);
        out += this_t_count;
        *out_daft += this_t_count_daft;
        LOut_Pad--;
        LOut() << "T-count for partition " << this_part << " = " << this_t_count << endl;
        LOut() << "Total T-count so far = " << out << endl;
        LOut() << "T-gate per gate for this partition " << this_part << " = " << (double)this_t_count/(double)reduces_gates << endl;
        LOut() << "Hadamard ancillas for partition " << this_part << " = " << inC->hadamard_mode_max_ancillas << endl;
        LOut_Pad++;
        LOut() << "..." << endl;
        LOut_Pad--;
        /*
        if(!g_algorithm.compare(DAFT_GUESS)) {

        } else if(!g_algorithm.compare(LEMPEL_GREEDY)) {

        } else if(!g_algorithm.compare(LEMPEL_LEAST_GREEDY)) {

        } else if(!g_algorithm.compare(LEMPEL_RANDOM)) {

        } else if(!g_algorithm.compare(LEMPEL_X)) {

        } else if(!g_algorithm.compare(REED_MULLER)) {

        } else if(!g_algorithm.compare(ALL_LEMPEL_SELECTORS)) {

        } else if(!g_algorithm.compare(ALL_LEMPEL)) {

        }
        */
        this_part++;
    } while(!exit);
    if(n_parts) *n_parts = (this_part-1);
    if(n_hadancs) *n_hadancs = inC->hadamard_mode_max_ancillas;
    int toc = clock();
    double exec_time = ((double)toc-(double)tic)/CLOCKS_PER_SEC;
    if(out_exec_time) *out_exec_time = exec_time;
    LOut_Pad--;
    LOut() << "UniversalTCount: Algorithm end." << endl << endl;
    return out;
}
