#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#include "kernel_header.h"
#include "CL/cl.h"

enum oclFFTStatus : uint32_t {
	oclFFT_SUCCESS					=			0,
	oclFFT_FAILED					=			1,
	oclFFT_PLAN_EXECUTING			=			2
};

enum oclFFT_DIRECTION : int32_t {
	oclFFT_FORWARD		=		-1,
	oclFFT_INVERSE		=		1
};

class oclFFTPlan1D {
	private:
		int32_t FFT_LEN, FFT_BATCH, FFT_ChirpZ;
		uint32_t FFT_TYPES[6] = { 8, 7, 5, 4, 3, 2 };
		uint32_t Next_K_W[6] = { 0, 0, 0, 0, 0, 0 };
		uint32_t NCounts[6] = { 0, 0, 0, 0, 0, 0 };
		bool	baked;
		bool	executing;
		oclFFT_DIRECTION	direction;
		std::string 	fft8_kernel;
		std::string 	fft7_kernel;
		std::string 	fft5_kernel;
		std::string 	fft4_kernel;
		std::string 	fft3_kernel;
		std::string 	fft2_kernel;
		oclFFTHeader	PlanKernel;

		oclFFTStatus Calc_Counts();

	public:
		oclFFTPlan1D(uint32_t);
		oclFFTStatus Set_KW (uint32_t);
		oclFFTStatus Set_Batch (uint32_t);
		oclFFTStatus Gen_Kernel ();
		oclFFTStatus Set_Forward ();
		oclFFTStatus Set_Inverse ();
		uint32_t Get_Len() { return FFT_LEN; };
		uint32_t Get_Batch() { return FFT_BATCH; };
		uint32_t Get_Chirp_Len() { return FFT_ChirpZ; };
		uint32_t Get_Count(int32_t);
		uint32_t Get_Next_K_W(int32_t);
		int32_t Get_FFT_Type(int32_t);
		oclFFT_DIRECTION Get_Direction() { return direction; };
};

