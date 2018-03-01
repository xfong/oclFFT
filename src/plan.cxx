#include "plan.h"

#define LOCAL_GRP_SZ 256

oclFFTPlan1D::oclFFTPlan1D(uint32_t x) {
	FFT_BATCH = 1;
	FFT_LEN = (int32_t)(x);
	FFT_ChirpZ = -1;
	baked = false;
	executing = false;
	direction = oclFFT_FORWARD;
	fft_precision = oclFFT_FLOAT;
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

oclFFTStatus oclFFTPlan1D::Set_Float() {
	if (!executing) {
		fft_precision = oclFFT_FLOAT;
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::Set_Double() {
	if (!executing) {
		fft_precision = oclFFT_DOUBLE;
		return oclFFT_SUCCESS;
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::DiscoverChirpZ(int32_t *outVal) {
	*outVal = -1;
	if (!executing) {
		if (FFT_ChirpZ < 1) {
			return oclFFT_FAILED;
		} else if (FFT_ChirpZ == 1) {
			*outVal = 1;
			return oclFFT_SUCCESS;
		} else {
			bool flag = true;
			int32_t ChirpLength = FFT_ChirpZ + 1;
			do {
				oclFFTPlan1D tmpPlan(2*ChirpLength);
				int32_t	testValue = tmpPlan.Get_Chirp_Len();
				if (testValue == 1) {
					*outVal = 2 * ChirpLength;
					flag = false;
				} else {
					ChirpLength++;
				}
			} while (flag);
			return oclFFT_SUCCESS;
		}
	} else {
		return oclFFT_PLAN_EXECUTING;
	}
}

oclFFTStatus oclFFTPlan1D::Gen_Main_FFT_Kernel() {
	if (!executing) {
		if (FFT_LEN > 1) {
			PlanKernel.SetFFTLength(FFT_BATCH * FFT_LEN);
			PlanKernel.SetBatchCount(FFT_BATCH * FFT_ChirpZ);
			PlanKernel.SetIOFSScaleFactor(FFT_LEN);
			PlanKernel.SetLocalGroupSize(LOCAL_GRP_SZ);
			PlanKernel.SetGlobalGroupSize(LOCAL_GRP_SZ);
			if (direction == oclFFT_INVERSE) {
				PlanKernel.SetInverseDirection(true);
			} else {
				PlanKernel.SetInverseDirection(false);
			}
			main_fft_kernel = PlanKernel.ReturnCommonMacros();
			main_fft_kernel += PlanKernel.GetSourceHeader();
			main_fft_kernel += PlanKernel.ReturnMacro();
			bool first_kern = true;
			PlanKernel.SetInitKW(FFT_LEN);
			int32_t init_num = 0;
			for (uint8_t idx = 0; idx < 6; idx++) {
				if (NCounts[idx] > 0) {
					PlanKernel.SetKernelFFTLength(FFT_TYPES[idx]);
					if (first_kern) {
						init_num = 0;
						PlanKernel.SetInitSCount(init_num);
						PlanKernel.SetTwiddleSCount(init_num);
					} else {
						init_num = 1;
						PlanKernel.SetInitSCount(init_num);
						PlanKernel.SetTwiddleSCount(init_num - 1);
					}
					PlanKernel.SetOuterCount(NCounts[idx] + init_num);
					PlanKernel.SetInnerCount(FFT_BATCH * FFT_LEN / FFT_TYPES[idx]);
					PlanKernel.SetMatrixInverseSCount(NCounts[idx] - 1 + init_num);
					PlanKernel.SetFinalMatrixInverseSCount(NCounts[idx] + init_num);
					kernel_names[idx] = PlanKernel.genKernelName(8);
					first_kern = false;
					main_fft_kernel.append(PlanKernel.print_kernel_name(kernel_names[idx]));
					if (fft_precision == oclFFT_FLOAT) {
						main_fft_kernel.append(PlanKernel.print_kernel_float2_inputs(true));
					} else {
						main_fft_kernel.append(PlanKernel.print_kernel_float2_inputs(false));
					}
					main_fft_kernel.append(" {\n");
					main_fft_kernel.append(PlanKernel.print_kernel_initialization());
					main_fft_kernel.append(PlanKernel.print_kernel_outer_loop());
					main_fft_kernel.append("\n}\n\n");
				} else {
				}
				PlanKernel.SetInitKW(Next_K_W[idx]);
			}
			return oclFFT_SUCCESS;
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
