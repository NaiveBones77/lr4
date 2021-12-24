#pragma once
#include "Conversion.h"

int Conversion::conv_dec_to_bin(double max_value, int max_digit, int digit, double value)
{
	int* arr = massbin(max_value, max_digit, digit, value);
	int sum = 0;
	for (int i = 0; i < max_digit; ++i) {
		sum = sum + (arr[i] * pow(2, max_digit - (1 + i)));
	}
	delete[] arr;
	return sum;
}

int* Conversion::massbin(double max_value, int max_digit, int digit, double value) {
	int* arr = new int[max_digit]();
	int last_i = 0;
	float sum = 0;
	for (int i = 0; i <= digit - 1; i++)
	{
		sum = sum + (max_value / pow(2, i));
		if (sum <= value) {
			arr[i] = 1;
		}
		else {
			arr[i] = 0;
			sum = sum - (max_value / pow(2, i));
		}
		last_i = i;
	}
	if (last_i < max_digit - 1)
		for (int i = last_i; i < max_digit; ++i)
			arr[i] = 0;
	return arr;
}

double Conversion::conv_bin_to_dec(double max_value, int max_digit, int digit, int dec)
{
	double* bin = new double[max_digit]();
	int sum = 0;
	int i = 0;
	while (dec > 0)
	{
		bin[max_digit - (i + 1)] = (dec % 2);
		dec /= 2;
		i++;
	}
	i = 0;
	while (i <= digit - 1)
	{
		sum = sum + bin[i] * (max_value / pow(2, i));
		i++;
	}
	delete[] bin;
	return sum;
}

