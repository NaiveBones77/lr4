#pragma once
#include <mutex>
#include <WinSock2.h>
#include <iostream>
#include "Conversion.h"
#include <random>
#include "protocol.h"


class INS {
public:

	Val Latitude, Longitude, H, CourseTrue, Tungazh, List, VelocityNS, VelocityEW,
		VelocityVI, AccelerationX, AccelerationZ, AccelerationY;
	std::mutex mutex;
	char* buffer;
	char buff2send[4];
	bool isStart = false;
	SOCKET _s;
	sockaddr_in _destAddr;
	int buff_count;
	int count_nav;
	Conversion conv1;
	std::mt19937 generator;



	INS(Val Latitude, Val Longitude, Val H, Val CourseTrue,
		Val Tungazh, Val List, Val VelocityNS, Val VelocityEW, Val VelocityVI,
		Val AccelerationX, Val AccelerationZ, Val AccelerationY);

	INS();

	void test();

	void preparation();

	void navigation();

	void start();

	void sendPack();

	int bindPort(SOCKET s, sockaddr_in destAddr);

	template <class T>
	void fill_buff(T pack);
};