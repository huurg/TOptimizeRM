#include "SQC_Circuit.h"

void SQC_Circuit::convert_Y_to_XZ() {
    for(int t = 0; t < m; t++) {
        if(operator_list[t][0]==SQC_OPERATOR_Y) {
            SQC_Circuit temp(n);
            int op_X[] = {SQC_OPERATOR_X,operator_list[t][1]};
            int op_Z[] = {SQC_OPERATOR_Z,operator_list[t][1]};
            temp.AddOperator(op_X,2);
            temp.AddOperator(op_Z,2);
            ReplaceOperator(&temp, t,1);
        }
    }
}

void SQC_Circuit::expand_Toffn_to_Toff3(SQC_Circuit& out) const {

}
