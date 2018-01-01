#include <cstring>
#include <cstdint>
#include <iostream>

#include "plan.h"

using namespace std;

int main() {
	int32_t testInput = 68;
	oclFFTPlan1D plan1 (testInput);
	wcout << "When performing FFT of length " << plan1.Get_Len() << ", the library will perform:" << endl;
	for (int idx = 0; idx < 6; idx++) {
		wcout << "FFT Length: " << plan1.Get_FFT_Type(idx) << "; Stages: " << plan1.Get_Count(idx) << "; Next K_W: " << plan1.Get_Next_K_W(idx) << endl;
	}
	wcout << "Chirp Z length: " << plan1.Get_Chirp_Len() << endl;
	oclFFTPlan1D plan2 (1);
	uint32_t test = plan1.Get_Chirp_Len();
	while (test != 1) {
		plan2.Set_KW(2 * (test + 1));
		test = plan2.Get_Chirp_Len();
		wcout << "Chirp Z new length: " << plan2.Get_Len() << "; Chirp Z of new plan: " << test << endl;
	}
}
