
#ifndef SIMPLEALGORITHMS_H
#define SIMPLEALGORITHMS_H


#include <math.h>
#include "linearalgebra.h"

class simpleAlgorithms
{
public:
    static long double mZero(){
        long double temp{1}, result{1};
        while (result + 1.0 > 1) {
            result = temp;
            temp /= 2.0;
        }
        return result;
    }

    static long double getMax(const long double &a, const long double &b){
        if (a > b){
            return a;
        }else{
            return  b;
        }
    }

    static long double getMin(const long double &a, const long double &b){
        if (a < b){
            return a;
        }else{
            return  b;
        }
    }

    static vector orbitToXYZ(long double big_omega, long double i,
                             long double omega, long double a,
                             long double e, long double nu,
                             long double mu)
    {
        vector r, v;
        r.resize(3); v.resize(3);
        r[0] = (a*(1 - pow(e,2.)))/(1+e*std::cos(nu));            r[1] = 0;                                  r[2] = 0;
        v[0] = std::sqrt(mu/(a*(1 - pow(e,2.))))*e*std::sin(nu);  v[1] = std::sqrt(mu/(a*(1 - pow(e,2.))))*(1+e*std::cos(nu)); v[2] = 0;

        matrix A;
        A.resize(3, 3);
        long double u = omega + nu;
        A(0, 0) =  cos(u)*cos(big_omega) - sin(u)*sin(big_omega)*cos(i); A(0, 1) = -sin(u)*cos(big_omega) - cos(u)*sin(big_omega)*cos(i); A(0, 2) =  sin(i)*sin(big_omega);
        A(1, 0) =  cos(u)*sin(big_omega) + sin(u)*cos(big_omega)*cos(i); A(1, 1) = -sin(u)*sin(big_omega) - cos(u)*cos(big_omega)*cos(i); A(1, 2) = -sin(i)*sin(big_omega);
        A(2, 0) =  sin(u)*sin(i);                                        A(2, 1) =  cos(u)*sin(i);                                        A(2, 2) =  cos(i);

        vector tmp, final;
        tmp.resize(3); final.resize(6);
        tmp = A*r;
        for (int i = 0; i < 3; i++) final[i] = tmp[i];
        tmp = A*v;
        for (int i = 0; i < 3; i++) final[i+3] = tmp[i];
        return final;
    }
    static vector polarToXYZ(long double lat,long double lon, long double h){
        vector r(3);
        r[0] = h*cos(lat)*cos(lon);
        r[1] = h*cos(lat)*sin(lon);
        r[2] = h*sin(lat);
        return r;
    }
};

#endif // SIMPLEALGORITHMS_H
