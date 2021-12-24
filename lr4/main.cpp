#pragma once
#include "Conversion.h"
#include "INS.h"
#include "protocol.h"
#include "SNS.h"
#include "Timer.cpp"

#include <WinSock2.h>
#include "Aircraft.h"
#include "integrator.h"
#pragma comment(lib,"Ws2_32.lib")

#define PORT 12346
#define SERVERADDR "127.0.0.1"
#pragma warning(disable: 4996)

SOCKET _s;
sockaddr_in _destAddr;

void test_run(Aircraft& a)
{
    a.run2();
}

void test_ops(Aircraft& a)
{
    a.OPS2();
}

void test_bomb(Aircraft& a)
{
    a.startBomb();
}


int main()
{
	setlocale(LC_ALL, "rus");



	INS ins(Val("������", 28, 20, 90),
            Val("�������", 55, 20, 90),
            Val("������", 130, 19, 19975.3728),
            Val("���� ��������", 15.3, 16, 90),
            Val("������", 3.5, 16, 90),
            Val("����", 6.3245, 16, 90),
            Val("�������� ����� ��", 400, 19, 1053.5822),
            Val("�������� ������ �����", 200, 19, 1053.5822),
            Val("�������� ������������ ������������", 83, 19, 83.2307),
            Val("��������� ����������", 0, 12, 19.62),
            Val("��������� ����������", 0, 12, 19.62),
            Val("��������� ����������", 0, 12, 2)
    );

    SNS sns(Val("������", 10, 20, 65536),
            Val("HDOP", 10.0, 15, 512),
            Val("VDOP", 10.0, 15, 512),
            Val("������� ����", 5.0, 15, 90),
            Val("������� ������ ", 55, 20, 90),
            Val("������� ������ (�����)", 3, 11, 0.000085830),
            Val("������� �������", 35, 20, 90),
            Val("������� ������� (�����)", 4, 11, 0.000085830),
            Val("�������� ������", 13.3, 20, 512),
            Val("������� ����� UTC (������� �������)", 6.0, 6, 32),
            Val("������� ����� UTC (������� �������)", 2.0, 20, 512),
            Val("������������ ��������", 1.0, 15, 16384)
    );

    char buff[10 * 1014];

    if (WSAStartup(0x202, (WSADATA*)&buff[0]))
    {
        printf("WSAStartup error: %d\n",
            WSAGetLastError());
        return -1;
    }

    _s = socket(AF_INET, SOCK_DGRAM, 0);
    if (_s == INVALID_SOCKET)
    {
        printf("socket() error: %d\n", WSAGetLastError());
        WSACleanup();
        return -1;
    }

    HOSTENT* hst;
    _destAddr;
    _destAddr.sin_family = AF_INET;
    _destAddr.sin_port = htons(PORT);

    if (inet_addr(SERVERADDR))
        _destAddr.sin_addr.s_addr = inet_addr(SERVERADDR);

    else
    {
        if (hst = gethostbyname(SERVERADDR))
            _destAddr.sin_addr.s_addr = ((unsigned long**)hst->h_addr_list)[0][0];

        else
        {
            printf("Unknown host: %d\n", WSAGetLastError());
            closesocket(_s);
            WSACleanup();
            return -1;
        }
    }
	
	
	
    Timer timer;

    Aircraft a1(37.41255708413501, 55.97313079458042, 300, 0);

    a1.ins.bindPort(_s, _destAddr);
    a1.sns.bindPort(_s, _destAddr);

    //a1.run();
    

    timer.add(std::chrono::microseconds(1000), [&]() {test_run(a1); });
    timer.add(std::chrono::microseconds(1000), [&]() {test_ops(a1); });
    timer.add(std::chrono::microseconds(1000), [&]() {test_bomb(a1); });

    int i = 0;
    while (true) { std::this_thread::sleep_for(std::chrono::seconds(20));/* i++;*/ };
}

