#ifndef MATHMODEL_H
#define MATHMODEL_H

#include "linearalgebra.h"
#include <math.h>

class IMathModel
{
protected:
    long double sampIncrement, t_st, t_fin;//шаг, начало, конец
    matrix resMatrix;
    vector X0;//вектор начальных условий
    int N;//счетичк строк в результирующей матрице
public:
    IMathModel():sampIncrement(1e-1), t_st(0.0), t_fin(1.0), N(0){ resMatrix.resize(0, 0);}
    inline vector getInitialConditions() const {return X0;}
    inline vector setInitialConditions(vector X0) { this->X0 = X0; };
    inline int getOrder()  {return X0.size();}
    inline long double getSampIncrement() const {return this->sampIncrement;}
    inline long double getT_st() const {return this->t_st;}
    inline long double getT_fin() const {return this->t_fin;}
    inline void setT_st(long double t_st) {this->t_st = t_st;}
    inline void setT_fin(long double t_fin) {this->t_fin = t_fin;}
    inline void setSampIncrement(long double sampIncrement) {this->sampIncrement = sampIncrement;}
    inline matrix getResult() const {return this->resMatrix;}
    virtual void getRP(const vector &X, long double t, vector &Y) const = 0; //абстрактная функция правых частей
    virtual void prepareResult();
    virtual void clearResult();
    virtual void preStart();
    virtual void addResult(const vector &addVector, const long double &t);
};

#endif // MATHMODEL_H
