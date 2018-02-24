#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#include "kernel_header.h"
#include "CL/cl.h"

enum oclFFTStatus : uint8_t {
	oclFFT_SUCCESS					=			0,
	oclFFT_FAILED					=			1,
	oclFFT_PLAN_EXECUTING			=			2
};

enum oclFFT_DIRECTION : int8_t {
	oclFFT_FORWARD		=		-1,
	oclFFT_INVERSE		=		1
};

enum oclFFT_PRECISION : uint8_t {
	oclFFT_FLOAT		=		1,
	oclFFT_DOUBLE		=		2
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
		oclFFT_PRECISION	fft_precision;
		std::string 	main_fft_kernel;
		std::string		chirpz_kernel;
		std::string		kernel_names[6] = { "", "", "", "", "", "" };
		oclFFTHeader	PlanKernel;

		oclFFTStatus Calc_Counts();

	public:
		oclFFTPlan1D(uint32_t);
		oclFFTStatus Set_KW (uint32_t);
		oclFFTStatus Set_Batch (uint32_t);
		oclFFTStatus Gen_Main_FFT_Kernel ();
		oclFFTStatus Set_Forward ();
		oclFFTStatus Set_Inverse ();
		oclFFTStatus Set_Float ();
		oclFFTStatus Set_Double ();
		uint32_t Get_Len() { return FFT_LEN; };
		uint32_t Get_Batch() { return FFT_BATCH; };
		uint32_t Get_Chirp_Len() { return FFT_ChirpZ; };
		uint32_t Get_Count(int32_t);
		uint32_t Get_Next_K_W(int32_t);
		int32_t Get_FFT_Type(int32_t);
		oclFFT_DIRECTION Get_Direction() { return direction; };
		oclFFT_PRECISION Get_Precision() { return fft_precision; };
		std::string Get_Main_FFT_Kernel() { return main_fft_kernel; };
		std::string Get_ChirpZ_Kernel() { return chirpz_kernel; };
};

