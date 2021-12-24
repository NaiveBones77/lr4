#pragma once
#include <vector>
#include "Transition.h"
#include "mathmodel.h"
#include "simplealgorithms.h"
#include <string>
#include <fstream>;
#include <iostream>;

class ASP: public IMathModel
{

public:
	double T_h, H=10000, A;
	double g = 9.81;
	double t = 0;
	double dt = 0.01;
	double c = 0.9;
	double ro = 0.414;
	double s = 0.06;

	std::string name;


	Transition tr;


	std::vector <double> coordinates = {};
	std::vector <std::vector<double>> coordinatesList = {}, coordinatesG = {};
	std::vector <double> V;
	bool flag = true;
	void run();
	ASP(std::vector <double> coords, std::vector <double> V0);
	ASP();
	void setInitialConditions(std::vector<float>);
	void Ab_f();
	void setX(std::vector <double> X);
	void setV(std::vector <double> V);

	void insertGeogr();

	void getRP(const vector& X, long double t, vector& Y) const override;

	void setName(int index);

	void writeKML(std::vector<std::vector<double>> cG);

	//void setName(int i);
	
};

