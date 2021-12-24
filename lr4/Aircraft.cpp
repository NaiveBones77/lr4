#pragma once
#include "Aircraft.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <thread>
#include "Timer.cpp"
#include <list>
#include "integrator.h";


Aircraft::Aircraft(double longitude, double latitude, double V0, double A0) {

	this->longitude = longitude;
	this->latitude = latitude;
	V = V0;
	A = A0;
	roll = 0;
	pitch = 0;
	yaw = 0;
	Vx = 0;
	Vz = 0;
	tr.setDefault(longitude, latitude);
	PPMs.insert(PPMs.end(), std::vector<double>{20000, 0, 20000});
	PPMs.insert(PPMs.end(), std::vector<double>{20000, 0, 50000});
	//PPMs.insert(PPMs.end(), std::vector<double>{10000, 0, 100000});
	//PPMs.insert(PPMs.end(), std::vector<double>{-10000, 0, 150000});
	PPMsG.insert(PPMsG.end(), tr.fromStart2Geogr(std::vector<double>{20000, 0, 20000}));
	PPMsG.insert(PPMsG.end(), tr.fromStart2Geogr(std::vector<double>{20000, 0, 50000}));
	PPMsG.insert(PPMsG.end(), tr.fromStart2Geogr(std::vector<double>{10000, 0, 100000}));
	PPMsG.insert(PPMsG.end(), tr.fromStart2Geogr(std::vector<double>{-10000, 0, 150000}));
	tr.WritePPM(file1, std::string("PPM1.kml"), PPMsG[0]);
	tr.WritePPM(file1, std::string("PPM2.kml"), PPMsG[1]);
	tr.WritePPM(file1, std::string("PPM3.kml"), PPMsG[2]);
	tr.WritePPM(file1, std::string("PPM4.kml"), PPMsG[3]);
	ASP bomb0(startSK, std::vector <double> {Vx, 0, Vz});
	ASP bomb1(startSK, std::vector <double> {Vx, 0, Vz});
	ASP bomb2(startSK, std::vector <double> {Vx, 0, Vz});
	ASP bomb3(startSK, std::vector <double> {Vx, 0, Vz});
	bombs.push_back(bomb0);
	bombs.push_back(bomb1);
	bombs.push_back(bomb2);
	bombs.push_back(bomb3);

	SNS sns(Val("высота", 10, 20, 65536),
			Val("HDOP", 10.0, 15, 512),
			Val("VDOP", 10.0, 15, 512),
			Val("путевой угол", 5.0, 15, 90),
			Val("текущая широта ", 55, 20, 90),
			Val("текущая широта (точно)", 3, 11, 0.000085830),
			Val("текущая долгота", 35, 20, 90),
			Val("текущая долгота (точно)", 4, 11, 0.000085830),
			Val("задержка выдачи", 13.3, 20, 512),
			Val("текущее время UTC (старшие разряды)", 6.0, 6, 32),
			Val("текущее время UTC (младшие разряды)", 2.0, 20, 512),
			Val("вертикальная скорость", 1.0, 15, 16384)
	);

	INS ins(Val("широта", 28, 20, 90),
			Val("долгота", 55, 20, 90),
			Val("высота", 130, 19, 19975.3728),
			Val("курс истинный", 15.3, 16, 90),
			Val("тангаж", 3.5, 16, 90),
			Val("крен", 6.3245, 16, 90),
			Val("скорость север юг", 400, 19, 1053.5822),
			Val("скорость восток запад", 200, 19, 1053.5822),
			Val("скорость вертикальная инерциальная", 83, 19, 83.2307),
			Val("ускорение продольное", 0, 12, 19.62),
			Val("ускорение поперечное", 0, 12, 19.62),
			Val("ускорение нормальное", 0, 12, 2)
	);

	dp_integrator = new dormandPrinceIntgrator();

	file1.open("123.kml", std::ios::out);
	tr.WriteFile(file1, 1, std::vector<double> {0, 0, 0});
	file1.close();

}

Aircraft::Aircraft(){ }



void Aircraft::run2()
{
	double dt = 0.1;
	
	coordinates.insert(coordinates.end(), startSK);
	if (index < PPMs.size())
	{
		mutex.lock();
		A = A + Xpr[0] * dt;
		Vx = V * cos(A);
		Vz = V * sin(A);
		startSK[0] = startSK[0] + Vx * dt;
		startSK[2] = startSK[2] + Vz * dt;
		coordinatesG = tr.fromStart2Geogr(startSK);
		latitude = coordinatesG[0];
		longitude = coordinatesG[1];

		if (tr.getDistance(std::vector<double> {startSK[0], startSK[2]}, std::vector<double>{PPMs[index][0], PPMs[index][2]}) < 500)
		{
			index += 1;
		}
		countOperation += 1;
		mutex.unlock();
	}
	std::vector<double> SNS_vec = { A, longitude, latitude };
	std::vector<double> INS_vec = { latitude, longitude, A, pitch, roll, Vx, Vz };
	fillSNS(SNS_vec);
	fillINS(INS_vec);

	if (countOperation % 50 == 0)
	{
		file1.open("123.kml", std::ios::app);
		tr.WriteFile(file1, 2, std::vector<double> {latitude, longitude, 10000});
		file1.close();
	}
	if (index == PPMs.size())
	{
		file1.open("123.kml", std::ios::app);
		tr.WriteFile(file1, 3, std::vector<double> {latitude, longitude, 10000});
		file1.close();
		index = 1000;
	}
}

void Aircraft::OPS2()
{
	mutex.lock();
	theta = asin(startSK[1] / tr.getDistance(startSK, PPMs[index]));
	std::vector<double> ort = { 0, 0 };
	ort[0] = (distSP[0] - startSK[0]) / tr.getDistance(distSP, startSK);
	ort[1] = (distSP[1] - startSK[1]) / tr.getDistance(distSP, startSK);

	double delta = tr.getAngleFromScalars(ort, std::vector<double> {PPMs[index][0] - startSK[0], PPMs[index][2] - startSK[2]});
	if (abs(delta - A) <= 0.011)
	{
		Xpr[0] = 0;
	}

	if (delta > A)
	{
		Xpr[0] = yawmax_pr;
		AList.insert(AList.end(), tr.fromRad2Grad(A));
		thetaList.insert(thetaList.end(), tr.fromRad2Grad(theta));
	}

	else if (delta < A)
	{
		Xpr[0]  = -yawmax_pr;
		AList.insert(AList.end(), tr.fromRad2Grad(A));
		thetaList.insert(thetaList.end(), tr.fromRad2Grad(theta));
	}
	
	mutex.unlock();
}

void Aircraft::startBomb()
{
	mutex.lock();
	Timer timer;
	if ((curIndexBomb < bombs.size()) && bombs[curIndexBomb].flag)
	{
		if (countOperation > 900)
		{
			printf("590");
		}
		ASP bomb = bombs[curIndexBomb];
		bomb.setX(startSK);
		bomb.setV(std::vector<double> {V*cos(A), 0, V*sin(A)});


		std::vector<float> initial{ (float)startSK[0], (float)startSK[1], (float)startSK[2], (float)Vx, 0, (float)Vz };
		bomb.Ab_f();
		bomb.setInitialConditions(initial);
		bomb.setSampIncrement(0.1);
		bomb.setT_st(0);
		bomb.setT_fin(bomb.T_h);

		dp_integrator->run(&bomb);

		bomb.insertGeogr();

		double a = bomb.coordinatesList.back()[0];
		if ((abs(tr.getDistance(std::vector<double> {bomb.coordinatesList.back()[0], bomb.coordinatesList.back()[2]},
			std::vector<double> {PPMs[index][0], PPMs[index][2]})) < 1000) && (curIndexBomb==index))
		{			
			
			bomb.setName(index);
			tr.WriteBomb(file1, bomb.name, bomb.coordinatesG);

			bombs[curIndexBomb].flag = false;
			if (curIndexBomb > bombs.size())
			{
				
			}
			else
				curIndexBomb += 1;
		}
	}
	if (curIndexBomb == 2)
	{
		tr.writeAngles(file1, std::string("Angles.txt"), thetaList, AList);
		curIndexBomb = 1000;
	}
		
	mutex.unlock();
	
}

void Aircraft::run_asp(ASP& asp)
{
	asp.run();
}

void Aircraft::fillSNS(std::vector<double> vec)
{
	sns.mutex.lock();
	sns.H.value = 10000;
	sns.trackAngle.value = vec[0];
	sns.curLatitude.value = vec[1];
	sns.curLongitude.value = vec[2];
	sns.mutex.unlock();
}

void Aircraft::fillINS(std::vector<double> vec)
{
	ins.mutex.lock();
	ins.H.value = 10000;
	ins.Latitude.value = vec[0];
	ins.Longitude.value = vec[1];
	ins.CourseTrue.value = vec[2];
	ins.Tungazh.value = vec[3];
	ins.List.value = vec[4];
	ins.VelocityNS.value = vec[5];
	ins.VelocityEW.value = vec[6];
	ins.mutex.unlock();
}