#pragma once
#include "integrator.h"

dormandPrinceIntgrator::dormandPrinceIntgrator(): IIntegrator(){
    b1.resize(7);
    b2.resize(7);
    c.resize(7);
    K.resize(7);
    a.resize(7,6);

    b1[0] = 35./384;     b1[1] = 0.0; b1[2] = 500./1113;   b1[3] = 125./192; b1[4] = -2187./6784;    b1[5] = 11./84;    b1[6] = 0.0;

    b2[0] = 5179./57600; b2[1] = 0.0; b2[2] = 7571./16695; b2[3] = 393./640; b2[4] = -92097./339200; b2[5] = 187./2100; b2[6] = 1./40;

    c[0] = 0.0;          c[1] = 3.0;  c[2] = 3./10 ;       c[3] = 4./5;      c[4] = 8./9;            c[5] = 1.0 ;       c[6] = 1.0;

    a(0,0) = 0.0;         a(0,1) = 0.0;          a(0,2) = 0.0;           a(0,3) = 0.0;         a(0,4) = 0.0;          a(0,5) = 0.0;
    a(1,0) = 1./5;        a(1,1) = 0.0;          a(1,2) = 0.0;           a(1,3) = 0.0;         a(1,4) = 0.0;          a(1,5) = 0.0;
    a(2,0) = 3./40;       a(2,1) = 9./40;        a(2,2) = 0.0;           a(2,3) = 0.0;         a(2,4) = 0.0;          a(2,5) = 0.0;
    a(3,0) = 44./45;      a(3,1) = -56./15;      a(3,2) = 32./9;         a(3,3) = 0.0;         a(3,4) = 0.0;          a(3,5) = 0.0;
    a(4,0) = 19372./6561; a(4,1) = -25360./2187; a(4,2) = 64448./6561;   a(4,3) = -212./729;   a(4,4) = 0.0;          a(4,5) = 0.0;
    a(5,0) = 9017./3168;  a(5,1) = -355./33;     a(5,2) = 46732./5247;   a(5,3) = 49./176;     a(5,4) = -5103./18656; a(5,5) = 0.0;
    a(6,0) = 35./384;     a(6,1) = 0.0;          a(6,2) = 500./1113;     a(6,3) = 125./192;    a(6,4) = -2187./6784;  a(6,5) = 11./84;

    zero = simpleAlgorithms::mZero();
}

void dormandPrinceIntgrator::run(IMathModel * model)
{
    std::ofstream file;
    file.open("results.txt");
    long double
            t = model->getT_st(),
            t_out = t,
            t1 = model->getT_fin(),
            h,
            h_new = model->getSampIncrement(),
            e = 0;

    vector
            X = model->getInitialConditions(),
            X1( X.size() ),
            X2( X.size() ),
            Xout ( X.size() ),
            Y( X.size() );
    model->clearResult();
    model->prepareResult();
    for(auto &elem : K){
        elem.resize(X.size());
    }

    while ( t < t1 )
    {
        h = h_new;
        for (int j = 0; j < (int)K.size(); ++j) {
            for (int k = 0; k < X.size(); ++k) {
                Y[k] = X[k];
                for (int s = 0; s < j; ++s) {
                    Y[k] += K[s][k] * a(j,s) * h;
                }
            }
            model->getRP(Y, t, K[j]);
        }
        e = 0;
        for (int i = 0; i < X.size(); ++i){
            X1[i] = X2[i] = X[i];
            for (int j = 0; j < b1.size(); ++j) {
                X1[i] += K[j][i] * b1[j] * h;
                X2[i] += K[j][i] * b2[j] * h;
            }
            e += powl(h * (X1[i] - X2[i]) /
                      simpleAlgorithms::getMax(simpleAlgorithms::getMax(abs(X[i]), fabsl(X1[i])), simpleAlgorithms::getMax(1e-5L, 2 * zero / eps)), 2.0);
        }
        e = sqrtl( e / X.size() );
        h_new = h / simpleAlgorithms::getMax( 0.1, simpleAlgorithms::getMin( 5., pow(e / eps, 0.2)/0.9 ) );

        if ( e > eps )
            continue;
        while ( (t_out < t + h) && (t_out <= t1) )
        {
            long double theta = (t_out - t)/h, b[6];

            b[0] = theta  * ( 1 + theta    *(-1337./480.  + theta*(1039./360.    + theta*(-1163./1152.))));
            b[1] = 0;
            b[2] = 100.   * powl(theta, 2) * (1054./9275. + theta*(-4682./27825. + theta*(379./5565.)))/3.;
            b[3] = -5.    * powl(theta, 2) * (27./40.     + theta*(-9./5.        + theta*(83./96.)))/2.;
            b[4] = 18225. * powl(theta, 2) * (-3./250.    + theta*(22./375.      + theta*(-37./600.)))/848.;
            b[5] = -22.   * powl(theta, 2) * (-3./10.     + theta*(29./30.       + theta*(-17./24.)))/7.;

            for ( int k = X.size()-1; k >= 0; k-- )
            {
                long double sum  = 0;
                for ( int j = 5; j >= 0; j-- )
                    sum += b[j] * K[j][k];
                Xout[k] = X[k] + h * sum;
            }

            model->addResult( Xout, t_out );
            for(int i = 0; i < Xout.size(); i++){
                file<<Xout[i]<<"|";
            }
            file<<t_out<<std::endl;
            t_out += model->getSampIncrement();
        }
        X = X1;
        t += h;
    }
}
