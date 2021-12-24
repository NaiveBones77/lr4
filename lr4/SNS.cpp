#pragma once
#include "SNS.h"
#include <mutex>
#include <WinSock2.h>
#include <iostream>]

SNS::SNS(Val H, Val HDOP, Val VDOP, Val trackAngle, Val curLatitude, Val curLatitude_pr,
    Val curLongitude, Val curLongitude_pr, Val delay, Val UTC_high, Val UTC_low,
    Val V_h)
{
	this->H = H;
	this->HDOP = HDOP;
	this->VDOP = VDOP;
	this->trackAngle = trackAngle;
	this->curLatitude = curLatitude;
	this->curLatitude_pr = curLatitude_pr;
	this->curLongitude = curLongitude;
	this->curLongitude_pr = curLongitude_pr;
	this->delay = delay;
	this->UTC_high = UTC_high;
	this->UTC_low = UTC_low;
	this->V_h = V_h;
}

SNS::SNS() {};


void SNS::test() {
	std::cout << "test control start\n";
	mutex.lock();
	auto start = std::chrono::system_clock::now();

	FEATURES pack;
	pack.value.address = 187;
	pack.value.request_data = 0;
	pack.value.type_srns = 3;
	pack.value.almanach_gps = 1;
	pack.value.almanach_glonass = 1;
	pack.value.time_ind = 0;
	pack.value.emp1 = 0;
	pack.value.mode = 0;
	pack.value.diff_mode = 0;
	pack.value.submode_srns = 0;
	pack.value.emp2 = 0;
	pack.value.rejection = 1;
	pack.value.signal_threshold = 0;
	pack.value.cs = 0;
	pack.value.emp3 = 0;
	pack.value.sm = 0;
	pack.value.p = 0;

	buffer = new char[sizeof(pack)];
	memset(buffer, 0, sizeof(pack));
	memcpy(buffer, &pack, sizeof(pack));

	mutex.unlock();


	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	while (elapsed_seconds.count() < 12) {
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
	};

	mutex.lock();
	pack.value.mode = 1;
	pack.value.submode_srns = 1;
	pack.value.rejection = 0;
	pack.value.p = 1;


	buffer = new char[sizeof(pack)];
	memset(buffer, 0, sizeof(pack));
	memcpy(buffer, &pack, sizeof(pack));

	mutex.unlock();
	std::cout << "test control end\n";
}

void SNS::preparation() {
	std::cout << "prepare start\n";

	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	while (elapsed_seconds.count() < 5) {
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
	};

	std::cout << "prepare end\n";
}

void SNS::navigation() {
    if (count_nav == 0)
    {
        std::cout << "navigation start\n";
        delete[] buffer;
        buffer = new char[14 * sizeof(BNR)];
        count_nav++;
    }

    std::normal_distribution<double> delta(0.0, 0.01);

    double h = H.value + delta(generator);
    double hDOP = HDOP.value + delta(generator);
    double vDOP = VDOP.value + delta(generator);
    double trackangle = trackAngle.value + delta(generator);
    double curlatitude = curLatitude.value + delta(generator);
    double curlatitude_pr = curLatitude_pr.value + delta(generator);
    double curlongitude = curLongitude.value + delta(generator);
    double curlongitude_pr = curLongitude_pr.value + delta(generator);
    double delay_ = delay.value + delta(generator);
    double utc_high = UTC_high.value + delta(generator);
    double utc_low = UTC_low.value + delta(generator);
    double v_h = V_h.value + delta(generator);

    mutex.lock();

    buff_count = 0;


    BNR pack_1;
    pack_1.value.address = 62;
    pack_1.value.high = conv1.conv_dec_to_bin(H.max_value, 20, H.height, h);
    pack_1.value.sign = 1;
    pack_1.value.sm = 3;
    pack_1.value.p = 1;

    fill_buff<BNR>(pack_1);

    BNR pack_2;
    pack_2.value.address = 65;
    pack_2.value.high = conv1.conv_dec_to_bin(HDOP.max_value, 20, HDOP.height, hDOP);
    pack_2.value.sign = 1;
    pack_2.value.sm = 3;
    pack_2.value.p = 1;

    fill_buff<BNR>(pack_2);

    BNR pack_3;
    pack_3.value.address = 66;
    pack_3.value.high = conv1.conv_dec_to_bin(VDOP.max_value, 20, VDOP.height, vDOP);
    pack_3.value.sign = 1;
    pack_3.value.sm = 3;
    pack_3.value.p = 1;

    fill_buff<BNR>(pack_3);

    BNR pack_4;
    pack_4.value.address = 67;
    pack_4.value.high = conv1.conv_dec_to_bin(trackAngle.max_value, 20, trackAngle.height, trackangle);
    pack_4.value.sign = 1;
    pack_4.value.sm = 3;
    pack_4.value.p = 1;

    fill_buff<BNR>(pack_4);

    BNR pack_5;
    pack_5.value.address = 72;
    pack_5.value.high = conv1.conv_dec_to_bin(curLatitude.max_value, 20, curLatitude.height, curlatitude);
    pack_5.value.sign = 1;
    pack_5.value.sm = 3;
    pack_5.value.p = 1;

    fill_buff<BNR>(pack_5);

    BNR pack_6;
    pack_6.value.address = 80;
    pack_6.value.high = conv1.conv_dec_to_bin(curLatitude_pr.max_value, 20, curLatitude_pr.height, curlatitude_pr);
    pack_6.value.sign = 1;
    pack_6.value.sm = 3;
    pack_6.value.p = 1;

    fill_buff<BNR>(pack_6);

    BNR pack_7;
    pack_7.value.address = 73;
    pack_7.value.high = conv1.conv_dec_to_bin(curLongitude.max_value, 20, curLongitude.height, curlongitude);
    pack_7.value.sign = 1;
    pack_7.value.sm = 3;
    pack_7.value.p = 1;

    fill_buff<BNR>(pack_7);

    BNR pack_8;
    pack_8.value.address = 81;
    pack_8.value.high = conv1.conv_dec_to_bin(curLongitude_pr.max_value, 20, curLongitude_pr.height, curlongitude_pr);
    pack_8.value.sign = 1;
    pack_8.value.sm = 3;
    pack_8.value.p = 1;

    fill_buff<BNR>(pack_8);

    BNR pack_9;
    pack_9.value.address = 75;
    pack_9.value.high = conv1.conv_dec_to_bin(delay.max_value, 20, delay.height, delay_);
    pack_9.value.sign = 1;
    pack_9.value.sm = 3;
    pack_9.value.p = 1;

    fill_buff<BNR>(pack_9);

    BNR pack_10;
    pack_10.value.address = 104;
    pack_10.value.high = conv1.conv_dec_to_bin(UTC_high.max_value, 20, UTC_high.height, utc_high);
    pack_10.value.sign = 1;
    pack_10.value.sm = 3;
    pack_10.value.p = 1;

    fill_buff<BNR>(pack_10);

    BNR pack_11;
    pack_11.value.address = 96;
    pack_11.value.high = conv1.conv_dec_to_bin(UTC_low.max_value, 20, UTC_low.height, utc_low);
    pack_11.value.sign = 1;
    pack_11.value.sm = 3;
    pack_11.value.p = 1;

    fill_buff<BNR>(pack_11);

    BNR pack_12;
    pack_12.value.address = 117;
    pack_12.value.high = conv1.conv_dec_to_bin(V_h.max_value, 20, V_h.height, v_h);
    pack_12.value.sign = 1;
    pack_12.value.sm = 3;
    pack_12.value.p = 1;

    fill_buff<BNR>(pack_12);

    DATE_ pack_13;
    pack_13.value.address = 176;
    pack_13.value.emp1 = 0;
    pack_13.value.y_u = Year - Year / 10 * 10;
    pack_13.value.y_d = Year / 10;
    pack_13.value.m_u = Mounth - Mounth / 10 * 10;
    pack_13.value.m_d = Mounth / 10;
    pack_13.value.d_u = Day - Day / 10 * 10;
    pack_13.value.d_d = Day / 10;
    pack_13.value.emp2 = 1;
    pack_13.value.sm = 0;
    pack_13.value.p = 1;

    fill_buff<DATE_>(pack_13);

    FEATURES pack_14;
    pack_14.value.address = 187;
    pack_14.value.request_data = 0;
    pack_14.value.type_srns = 3;
    pack_14.value.almanach_gps = 1;
    pack_14.value.almanach_glonass = 1;
    pack_14.value.time_ind = 0;
    pack_14.value.emp1 = 0;
    pack_14.value.mode = 1;
    pack_14.value.diff_mode = 0;
    pack_14.value.submode_srns = 1;
    pack_14.value.emp2 = 0;
    pack_14.value.rejection = 0;
    pack_14.value.signal_threshold = 0;
    pack_14.value.cs = 0;
    pack_14.value.emp3 = 0;
    pack_14.value.sm = 0;
    pack_14.value.p = 1;

    fill_buff<FEATURES>(pack_14);

    mutex.unlock();
}

void SNS::start() {
	//if (not isStart)
	//{
	//	test();
	//	preparation();
	//	isStart = true;
	//}
    isStart = true;
	navigation();
}

void SNS::sendPack()
{
	mutex.lock();
	if (!isStart)
	{
		sendto(_s, &buffer[0], sizeof(buffer), 0,           
			(sockaddr*)&_destAddr, sizeof(_destAddr));		
															
	}
	else
	{
		char buff[4];
		int count_word = 0;
		for (int i = 0; i < 14; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				buff[j] = buffer[count_word * 4 + j];
			}

			sendto(_s, &buff[0], sizeof(buff), 0,
				(sockaddr*)&_destAddr, sizeof(_destAddr));
			count_word++;
		}

	}

	mutex.unlock();
	std::cout << "send " << std::endl;
}

int SNS::bindPort(SOCKET s, sockaddr_in destAddr) {
	_s = s;
	_destAddr = destAddr;
	return 1;
}

template <class T>
void SNS::fill_buff(T pack)
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
