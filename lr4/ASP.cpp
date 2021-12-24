#include "ASP.h"
#include <cmath>
#include <string>


void ASP::run()
{
	int ind = 0;
	while (H > 0)
	{
		//coordinates[1] = coordinates[1] + V[1]*dt - g * pow(t, 2) *dt / 2;
		V[1] = V[1] - g * dt;

		coordinates[1] = coordinates[1] + V[1]*dt - g * pow(dt, 2)  / 2;
		H = coordinates[1];
		coordinates[0] = coordinates[0] + V[0] * dt;
		coordinates[2] = coordinates[2] + V[2] * dt;
		/*V[1] = V[1] - g * dt;*/
		coordinatesList.insert(coordinatesList.end(), coordinates);
		coordinatesG.insert(coordinatesG.end(), tr.fromStart2Geogr(coordinates));
		t += dt;
		ind += 1;
	}

}

ASP::ASP(std::vector<double> coords, std::vector<double> V0)
{
	this->coordinates = coords;
	this->V = V0;
}

ASP::ASP() {};

void ASP::setInitialConditions(std::vector<float> vec)
{
	X0.resize(6);
	X0[0] = vec[0];
	X0[1] = vec[1];
	X0[2] = vec[2];
	X0[3] = vec[3];
	X0[4] = vec[4];
	X0[5] = vec[5];
}

void ASP::Ab_f()
{
	A = sqrt(pow(V[0], 2) + pow(V[2], 2)) * (sqrt(2 * H / g));
	T_h = (sqrt(2 * H / g));
}

void ASP::setX(std::vector<double> X)
{
	this->coordinates = X;
}

void ASP::setV(std::vector<double> V)
{
	this->V = V;
}

void ASP::insertGeogr()
{
	for (int i = 0; i < resMatrix.rowsCount(); i++)
	{
		auto X = resMatrix(i,0);
		auto Y = resMatrix(i,1);
		auto Z = resMatrix(i,2);
		std::vector<double> vec{ (double)X, (double)Y, (double)Z };
		coordinatesG.insert(coordinatesG.end(), tr.fromStart2Geogr(vec));
		coordinatesList.insert(coordinatesList.begin(), vec);
	}
}

void ASP::getRP(const vector& X, long double t, vector& Y) const
{
	Y.resize(X.size());
	Y[0] = X[3];
	Y[1] = X[4];
	Y[2] = X[5];
	Y[3] = -c*ro*X[3]*X[3]*s/2;
	Y[4] = -g;
	Y[5] = -c*ro*X[5]*X[5]*s/2;
}

void ASP::setName(int index)
{
	std::string name = "bomb";
	std::string ind = std::to_string(index);
	std::string endN = ".kml";
	this->name = name + ind + endN;
}