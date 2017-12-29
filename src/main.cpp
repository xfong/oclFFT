#include <cstring>
#include <cstdint>
#include <iostream>

#include "plan.hpp"

using namespace std;

int main() {
	int32_t testInput = 8*8*8*7*5*5*3*4*4*19;
	oclFFTPlan plan1 (testInput);
	wcout << "When performing FFT of length " << plan1.Get_Len() << ", the library will perform:" << endl;
	for (int idx = 0; idx < 6; idx++) {
		wcout << "FFT Length: " << plan1.Get_FFT_Type(idx) << "; Stages: " << plan1.Get_Count(idx) << "; Next K_W: " << plan1.Get_Next_K_W(idx) << endl;
	}
	wcout << "Chirp Z length: " << plan1.Get_Chirp_Len() << endl;
}