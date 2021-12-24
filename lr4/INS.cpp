#pragma once
#include "INS.h"
#include <mutex>
#include <WinSock2.h>
#include <iostream>
#include "protocol.h"
#include <cstdlib>
#include <random>

INS::INS(Val Latitude, Val Longitude, Val H, Val CourseTrue,
	Val Tungazh, Val List, Val VelocityNS, Val VelocityEW, Val VelocityVI,
	Val AccelerationX, Val AccelerationZ, Val AccelerationY) {
	this->Latitude = Latitude;
	this->Longitude = Longitude;
	this->H = H;
	this->CourseTrue = CourseTrue;
	this->Tungazh = Tungazh;
	this->List = List;
	this->VelocityNS = VelocityNS;
	this->VelocityEW = VelocityEW;
	this->VelocityVI = VelocityVI;
	this->AccelerationX = AccelerationX;
	this->AccelerationZ = AccelerationZ;
	this->AccelerationY = AccelerationY;
}

INS::INS() {};

void INS::test() {
	auto start = std::chrono::system_clock::now();
	double time_test = 5;

	mutex.lock();
	delete[] buffer;

	DISCRETE pack;   // Слово состояния с признаками "Исправность ИНС" и "Нет начальных данных"

	pack.value.address = 0270;
	pack.value.preparation = 0;
	pack.value.control = 0;
	pack.value.navigation = 0;
	pack.value.gyrocompas = 0;
	pack.value.restart = 0;
	pack.value.scale = 0;
	pack.value.heating = 1;
	pack.value.temp_control = 1;
	pack.value.init_data = 0;		// "Нет начальных данных"
	pack.value.no_reception = 0;
	pack.value.INS = 1;				// "Исправность ИНС"
	pack.value.acceleration = 0;
	pack.value.readiness = 0;
	pack.value.sm = 0;
	pack.value.p = 0;

	buffer = new char[sizeof(pack)];             
	memset(buffer, 0, sizeof(pack));            
	memcpy(buffer, &pack, sizeof(pack));     

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	while (elapsed_seconds.count() < time_test)
	{
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
	}

	delete[] buffer;
	buffer = new char[sizeof(pack)];            
	memset(buffer, 0, sizeof(pack));            
	memcpy(buffer, &pack, sizeof(pack));        

	mutex.unlock();
	std::cout << "test control end\n";
}

void INS::preparation() {
	mutex.lock();
	std::cout << "prepare start\n";

	DISCRETE pack;
	pack.value.address = 270;
	pack.value.sdi = 1;
	pack.value.preparation = 0;
	pack.value.control = 0;
	pack.value.navigation = 1;
	pack.value.gyrocompas = 0;
	pack.value.emp1 = 0;
	pack.value.restart = 0;
	pack.value.scale = 1;
	pack.value.heating = 1;
	pack.value.temp_control = 0;
	pack.value.init_data = 0;
	pack.value.no_reception = 1;
	pack.value.INS = 1;
	pack.value.acceleration = 0;
	pack.value.readiness = 0;
	pack.value.emp2 = 0;
	pack.value.sm = 0;
	pack.value.p = 1;

	mutex.unlock();

	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	int change = 1;
	double time_pr = 5;
	std::chrono::duration<double> elapsed_seconds = end - start;
	while (elapsed_seconds.count() < (time_pr)) {
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
		if (elapsed_seconds.count() / (time_pr / 6.) - change > 1)
		{
			change += 1;
			mutex.lock();
			pack.value.scale = change;
			mutex.unlock();
		};

	};

	mutex.lock();
	pack.value.scale = 7;
	std::cout << "preparation end\n";
	mutex.unlock();
}

void INS::navigation() {
	if (count_nav == 0)
	{
		std::cout << "navigation start\n";
		delete[] buffer;
		buffer = new char[13 * sizeof(BNR)];
		count_nav++;
	}

	
	std::normal_distribution<double>delta(0.0, 0.01);

	double latitude = Latitude.value + delta(generator);
	double longitude = Longitude.value + delta(generator);
	double h = H.value + delta(generator);
	double coursetrue = CourseTrue.value + delta(generator);
	double tungazh = Tungazh.value + delta(generator);
	double list = List.value + delta(generator);
	double velocityns = VelocityNS.value + delta(generator);
	double velocityew = VelocityEW.value + delta(generator);
	double velocityvi = VelocityVI.value + delta(generator);
	double accelerationx = AccelerationX.value + delta(generator);
	double accelerationz = AccelerationZ.value + delta(generator);
	double accelerationy = AccelerationY.value + delta(generator);

	mutex.lock();
	buff_count = 0;

	BNR pack_1;
	pack_1.value.address = 200;
	pack_1.value.high = conv1.conv_dec_to_bin(Latitude.max_value, 20, Latitude.height, latitude);
	pack_1.value.sign = 1;
	pack_1.value.sm = 3;
	pack_1.value.p = 1;

	fill_buff<BNR>(pack_1);

	BNR pack_2;
	pack_2.value.address = 201;
	pack_2.value.high = conv1.conv_dec_to_bin(Longitude.max_value, 20, Longitude.height, longitude);
	pack_2.value.sign = 1;
	pack_2.value.sm = 3;
	pack_2.value.p = 1;

	fill_buff<BNR>(pack_2);

	BNR pack_3;
	pack_3.value.address = 241;
	pack_3.value.high = conv1.conv_dec_to_bin(H.max_value, 20, H.height, h);
	pack_3.value.sign = 1;
	pack_3.value.sm = 3;
	pack_3.value.p = 1;

	fill_buff<BNR>(pack_3);

	BNR pack_4;
	pack_4.value.address = 204;
	pack_4.value.high = conv1.conv_dec_to_bin(CourseTrue.max_value, 20, CourseTrue.height, coursetrue);
	pack_4.value.sign = 1;
	pack_4.value.sm = 3;
	pack_4.value.p = 1;

	fill_buff<BNR>(pack_4);

	BNR pack_5;
	pack_5.value.address = 212;
	pack_5.value.high = conv1.conv_dec_to_bin(Tungazh.max_value, 20, Tungazh.height, tungazh);
	pack_5.value.sign = 1;
	pack_5.value.sm = 3;
	pack_5.value.p = 1;

	fill_buff<BNR>(pack_5);

	BNR pack_6;
	pack_6.value.address = 213;
	pack_6.value.high = conv1.conv_dec_to_bin(List.max_value, 20, List.height, list);
	pack_6.value.sign = 1;
	pack_6.value.sm = 3;
	pack_6.value.p = 1;

	fill_buff<BNR>(pack_6);

	BNR pack_7;
	pack_7.value.address = 246;
	pack_7.value.high = conv1.conv_dec_to_bin(VelocityNS.max_value, 20, VelocityNS.height, velocityns);
	pack_7.value.sign = 1;
	pack_7.value.sm = 3;
	pack_7.value.p = 1;

	fill_buff<BNR>(pack_7);

	BNR pack_8;
	pack_8.value.address = 247;
	pack_8.value.high = conv1.conv_dec_to_bin(VelocityEW.max_value, 20, VelocityEW.height, velocityew);
	pack_8.value.sign = 1;
	pack_8.value.sm = 3;
	pack_8.value.p = 1;

	fill_buff<BNR>(pack_8);

	BNR pack_9;
	pack_9.value.address = 245;
	pack_9.value.high = conv1.conv_dec_to_bin(VelocityVI.max_value, 20, VelocityVI.height, velocityvi);
	pack_9.value.sign = 1;
	pack_9.value.sm = 3;
	pack_9.value.p = 1;

	fill_buff<BNR>(pack_9);

	BNR pack_10;
	pack_10.value.address = 217;
	pack_10.value.high = conv1.conv_dec_to_bin(AccelerationX.max_value, 20, AccelerationX.height, accelerationx);
	pack_10.value.sign = 1;
	pack_10.value.sm = 3;
	pack_10.value.p = 1;

	fill_buff<BNR>(pack_10);

	BNR pack_11;
	pack_11.value.address = 218;
	pack_11.value.high = conv1.conv_dec_to_bin(AccelerationY.max_value, 20, AccelerationY.height, accelerationy);
	pack_11.value.sign = 1;
	pack_11.value.sm = 3;
	pack_11.value.p = 1;

	fill_buff<BNR>(pack_11);

	BNR pack_12;
	pack_12.value.address = 219;
	pack_12.value.high = conv1.conv_dec_to_bin(AccelerationZ.max_value, 20, AccelerationZ.height, accelerationz);
	pack_12.value.sign = 1;
	pack_12.value.sm = 3;
	pack_12.value.p = 1;

	fill_buff<BNR>(pack_12);

	DISCRETE pack_13;
	pack_13.value.address = 184;
	pack_13.value.sdi = 1;
	pack_13.value.preparation = 0;
	pack_13.value.control = 0;
	pack_13.value.navigation = 1;
	pack_13.value.gyrocompas = 0;
	pack_13.value.emp1 = 0;
	pack_13.value.restart = 0;
	pack_13.value.scale = 0;
	pack_13.value.heating = 1;
	pack_13.value.temp_control = 0;
	pack_13.value.init_data = 0;
	pack_13.value.no_reception = 1;
	pack_13.value.INS = 1;
	pack_13.value.acceleration = 0;
	pack_13.value.readiness = 0;
	pack_13.value.emp2 = 0;
	pack_13.value.sm = 0;
	pack_13.value.p = 1;

	fill_buff<DISCRETE>(pack_13);

	mutex.unlock();
}

void INS::start() {
	//if (not isStart)
	//{
	//	test();
	//	preparation();
	//	isStart = true;
	//}
	isStart = true;
	navigation();
}

int INS::bindPort(SOCKET s, sockaddr_in destAddr) {
	_s = s;
	_destAddr = destAddr;
	return 1;
}

void INS::sendPack()
{
	mutex.lock();
	if (!isStart)
	{
		sendto(_s, &buffer[0], sizeof(buffer), 0,           
			(sockaddr*)&_destAddr, sizeof(_destAddr));		
															
	}
	else
	{
		int count_word = 0;
		for (int i = 0; i < 13; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				buff2send[j] = buffer[count_word * 4 + j];
			}

			sendto(_s, &buff2send[0], sizeof(buff2send), 0,
				(sockaddr*)&_destAddr, sizeof(_destAddr));
			count_word++;
		}

	}

	mutex.unlock();
	std::cout << "send " << std::endl;
}

template <class T>
void INS::fill_buff(T pack)
{
	char* buff;
	buff = new char[sizeof(pack)];					 
	memset(buff, 0, sizeof(pack));                   
	memcpy(buff, &pack, sizeof(pack));				  

	for (int i = 0; i < 4; ++i)
	{
		buffer[4 * buff_count + i] = buff[i];
	}
	buff_count++;
	delete[] buff;
}

