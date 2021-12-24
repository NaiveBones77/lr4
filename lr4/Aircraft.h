#pragma once
#include "SNS.h"
#include "INS.h"
#include <vector>
#include "Transition.h"
#include <mutex>
#include <iostream>
#include <fstream>
#include "ASP.h"
#include <list>
#include "integrator.h";	
#include <string>



class Aircraft {
private:
	
	Transition tr;
	std::mutex mutex;
	std::ofstream file1;

	double roll, pitch, yaw;		//крен тангаж рысканье
	double longitude, latitude;		//долгота широта
	double A, theta;				//угол курса и элевации
	double V, Vx, Vz;				
	std::vector <double> startSK = { 0,10000,0 };
	std::vector <double> thetaList = { 0 };
	std::vector <double> AList = { 0 };
	std::vector <ASP> bombs;
	int curIndexBomb = 0;
	//std::list <ASP> bombs;
	//list <ASP> ::iterator it = bombs.begin();

	std::vector <std::vector<double>> coordinates;		//координаты ЛА в стартовой СК
	std::vector <double> coordinatesG;				//координаты ЛА в географической СК
	std::vector <std::vector<double>> PPMs;				//координаты ППМов в стартовой СК
	std::vector <std::vector<double>> PPMsG;				//координаты ППМов в географической СК

	std::vector <double> Xpr = {0};					//Произвоная по рысканью
	int index = 0;								//индекс элемента из ППМ
	double yawmax_pr = 0.14;
	int countPPM = 0;							//количество ППМОВ
	int countOperation = 0;						//количество операций

	double t = 0;
	std::vector<double> distSP = { 3885000, 0 };		//м до сев. полюса от текущей точки

	dormandPrinceIntgrator* dp_integrator;
public:
	INS ins;
	SNS sns;

	Aircraft();
	Aircraft(double longitude, double latitude, double V0, double A0);
		
	void run();
	void run2();

	std::vector<double> OPS(int index);
	void OPS2();
	void startBomb();
	void fillSNS(std::vector<double> vec);
	void fillINS(std::vector<double> vec);
	void run_asp(ASP& asp);

};


