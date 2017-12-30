#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <iostream>

#include "plan.hpp"

oclFFTPlan1D::oclFFTPlan1D(uint32_t x) {
	FFT_BATCH = 1;
	FFT_LEN = (int32_t)(x);
	FFT_ChirpZ = -1;
	Calc_Counts();
}

void oclFFTPlan1D::Calc_Counts() {
	bool flag;
	int32_t K_W = FFT_LEN;
	for (int idx = 0; idx < 6; idx++) {
		do {
			flag = true;
			uint32_t test = K_W % FFT_TYPES[idx];
			if ((test == 0) && (K_W > 0)) {
				K_W /= FFT_TYPES[idx];
				NCounts[idx] += 1;
			} else {
				if (K_W <= 0) {
					std::wcout << "Error!: K_W not a valid integer!" << std::endl;
				} else {
					Next_K_W[idx] = K_W;
				}
				flag = false;
			}
		} while (flag);
	}
	FFT_ChirpZ = K_W;
}

void oclFFTPlan1D::Set_KW(uint32_t K_W) {
	FFT_LEN = K_W;
	Calc_Counts();
}

void oclFFTPlan1D::Set_Batch(uint32_t x) {
	FFT_BATCH = x;
}

int32_t oclFFTPlan1D::Get_FFT_Type(int32_t idx) {
	if ((idx > 5) || (idx < 0)) {
		return 0;
	}
	return FFT_TYPES[idx];
}

uint32_t oclFFTPlan1D::Get_Count(int32_t idx) {
	if ((idx > 5) || (idx < 0)) {
		return 0;
	}
	return NCounts[idx];
}

uint32_t oclFFTPlan1D::Get_Next_K_W(int32_t idx) {
	if ((idx > 5) || (idx < 0)) {
		return 0;
	}
	return Next_K_W[idx];
}
