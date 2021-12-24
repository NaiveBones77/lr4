#include "mathmodel.h"

void IMathModel::addResult(const vector &addVector, const long double &t){
    for(int i = 0; i < getOrder(); i++){
        resMatrix(N, i) = addVector[i];
    }
    resMatrix(N, getOrder()) = t;
    ++N;
}

void IMathModel::prepareResult(){
    resMatrix.resize((int)((t_fin - t_st)/sampIncrement) + 1, getOrder() + 1);
}

void IMathModel::clearResult(){
    resMatrix.resize(0,0);
    N = 0;
}

void IMathModel::preStart(){

}
