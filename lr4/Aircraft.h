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

	double roll, pitch, yaw;		//���� ������ ��������
	double longitude, latitude;		//������� ������
	double A, theta;				//���� ����� � ��������
	double V, Vx, Vz;				
	std::vector <double> startSK = { 0,10000,0 };
	std::vector <double> thetaList = { 0 };
	std::vector <double> AList = { 0 };
	std::vector <ASP> bombs;
	int curIndexBomb = 0;
	//std::list <ASP> bombs;
	//list <ASP> ::iterator it = bombs.begin();

	std::vector <std::vector<double>> coordinates;		//���������� �� � ��������� ��
	std::vector <double> coordinatesG;				//���������� �� � �������������� ��
	std::vector <std::vector<double>> PPMs;				//���������� ����� � ��������� ��
	std::vector <std::vector<double>> PPMsG;				//���������� ����� � �������������� ��

	std::vector <double> Xpr = {0};					//���������� �� ��������
	int index = 0;								//������ �������� �� ���
	double yawmax_pr = 0.14;
	int countPPM = 0;							//���������� �����
	int countOperation = 0;						//���������� ��������

	double t = 0;
	std::vector<double> distSP = { 3885000, 0 };		//� �� ���. ������ �� ������� �����

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


