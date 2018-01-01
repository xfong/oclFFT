#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

class oclFFTPlan1D {
	private:
		int32_t FFT_LEN, FFT_BATCH, FFT_ChirpZ;
		uint32_t FFT_TYPES[6] = { 8, 7, 5, 4, 3, 2 };
		uint32_t Next_K_W[6] = { 0, 0, 0, 0, 0, 0 };
		uint32_t NCounts[6] = { 0, 0, 0, 0, 0, 0 };
		void Calc_Counts();
	public:
		oclFFTPlan1D(uint32_t);
		void Set_KW (uint32_t);
		void Set_Batch (uint32_t);
		uint32_t Get_Len() { return FFT_LEN; };
		uint32_t Get_Batch() { return FFT_BATCH; };
		uint32_t Get_Chirp_Len() { return FFT_ChirpZ; };
		uint32_t Get_Count(int32_t);
		uint32_t Get_Next_K_W(int32_t);
		int32_t Get_FFT_Type(int32_t);
};

