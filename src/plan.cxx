#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <iostream>

#include "plan.h"

#define LOCAL_GRP_SZ 256

oclFFTPlan1D::oclFFTPlan1D(uint32_t x) {
	FFT_BATCH = 1;
	FFT_LEN = (int32_t)(x);
	FFT_ChirpZ = -1;
	baked = false;
	executing = false;
	direction = oclFFT_FORWARD;
	Calc_Counts();
}

oclFFTStatus oclFFTPlan1D::Calc_Counts() {
	if (!executing) {
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
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::Set_KW(uint32_t K_W) {
	if (!executing) {
		FFT_LEN = K_W;
		baked = false;
		Calc_Counts();
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::Set_Batch(uint32_t x) {
	if (!executing) {
		FFT_BATCH = x;
		baked = false;
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
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

oclFFTStatus oclFFTPlan1D::Set_Forward() {
	if (!executing) {
		direction = oclFFT_FORWARD;
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::Set_Inverse() {
	if (!executing) {
		direction = oclFFT_INVERSE;
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::Gen_Kernel() {
	if (!executing) {
		if (FFT_LEN > 1) {
			PlanKernel.SetFFTLength(FFT_BATCH * FFT_LEN / FFT_ChirpZ);
			PlanKernel.SetKernelFFTLength(FFT_LEN / FFT_ChirpZ);
			PlanKernel.SetInitKW(FFT_LEN);
			PlanKernel.SetIOFSScaleFactor(FFT_LEN);
			PlanKernel.SetLocalGroupSize(LOCAL_GRP_SZ);
			PlanKernel.SetGlobalGroupSize(LOCAL_GRP_SZ);
			if (direction == oclFFT_INVERSE) {
				PlanKernel.SetInverseDirection(true);
			} else {
				PlanKernel.SetInverseDirection(false);
			}
		} else {
			if (FFT_LEN < 0) {
				return oclFFT_FAILED;
			} else {
				return oclFFT_SUCCESS;
			}
		}
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}