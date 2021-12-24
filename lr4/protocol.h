#pragma once
#include <map>
#include <iostream>

#pragma pack(push,1)
union BNR {
	struct data {
		unsigned address : 8;
		unsigned high : 20;
		unsigned sign : 1;
		unsigned sm : 2;
		unsigned p : 1;
	} value;
};
#pragma pack(pop)

#pragma pack(push,1)
union BCD {
	struct data {
		unsigned address : 8;
		unsigned emp1 : 3;
		unsigned sec : 6;
		unsigned min : 6;
		unsigned hour : 5;
		unsigned emp2 : 1;
		unsigned sm : 2;
		unsigned p : 1;
	} value;
};
#pragma pack(pop)

#pragma pack(push,1)
union DISCRETE {
	struct data {
		unsigned address : 8;
		unsigned sdi : 2;
		unsigned preparation : 1;
		unsigned control : 1;
		unsigned navigation : 1;
		unsigned gyrocompas : 1;
		unsigned emp1 : 1;
		unsigned restart : 1;
		unsigned scale : 3;
		unsigned heating : 1;
		unsigned temp_control : 1;
		unsigned init_data : 1;
		unsigned no_reception : 1;
		unsigned INS : 1;
		unsigned acceleration : 1;
		unsigned readiness : 1;
		unsigned emp2 : 3;
		unsigned sm : 2;
		unsigned p : 1;
	} value;
};
#pragma pack(pop)

#pragma pack(push,1)
union DATE_ {
	struct data {
		unsigned address : 8;
		unsigned emp1 : 2;
		unsigned y_u : 4;
		unsigned y_d : 4;
		unsigned m_u : 4;
		unsigned m_d : 1;
		unsigned d_u : 4;
		unsigned d_d : 1;
		unsigned emp2 : 1;
		unsigned sm : 2;
		unsigned p : 1;
	} value;
};
#pragma pack(pop)

#pragma pack(push,1)
union FEATURES {
	struct data {
		unsigned address : 8;
		unsigned request_data : 1;
		unsigned type_srns : 3;
		unsigned almanach_gps : 1;
		unsigned almanach_glonass : 1;
		unsigned mode : 2;
		unsigned submode_srns : 1;
		unsigned time_ind : 1;
		unsigned emp1 : 2;
		unsigned diff_mode : 1;
		unsigned emp2 : 1;
		unsigned rejection : 1;
		unsigned signal_threshold : 1;
		unsigned cs : 2;
		unsigned emp3 : 3;
		unsigned sm : 2;
		unsigned p : 1;
	} value;
};
#pragma pack(pop)

struct Val {
	std::string name;
	double value;
	int height;
	double max_value;
	Val(std::string name, double value, int height, double max_value)
	{
		this->name = name;
		this->value = value;
		this->height = height;
		this->max_value = max_value;
	}
	Val() {}
};