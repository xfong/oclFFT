#include <iostream>

#include "plan.h"

using namespace std;

int main(int argc, char* argv[]) {
	int32_t testInput;
	if (argc < 2) {
		testInput = 1;
		wcerr << "No inputs given: proper call is \"main <val>\"" << endl;
		return 1;
	}
	string input_arg = argv[1];
	testInput = stol(input_arg, NULL, 10);
	oclFFTPlan1D plan1 (testInput);
	wcout << "When performing FFT of length " << plan1.Get_Len() << ", the library will perform:" << endl;
	for (int idx = 0; idx < 6; idx++) {
		wcout << "FFT Length: " << plan1.Get_FFT_Type(idx) << "; Stages: " << plan1.Get_Count(idx) << "; Next K_W: " << plan1.Get_Next_K_W(idx) << endl;
	}
	wcout << "Chirp Z length: " << plan1.Get_Chirp_Len() << endl;
	oclFFTPlan1D plan2 (1);
	uint32_t test = plan1.Get_Chirp_Len() + 1;
	while (test != 1) {
		plan2.Set_KW(2 * test);
		uint32_t newLen = plan2.Get_Chirp_Len();
		wcout << "Chirp Z new length: " << plan2.Get_Len() << "; Chirp Z of new plan: " << newLen << endl;
		if (newLen != 1) {
			test++;
		} else {
			test = 1;
		}
	}
	int32_t ChirpPlanLen;
	plan1.DiscoverChirpZ(&ChirpPlanLen);
	wcout << "Chirp Z needs FFT length: " << ChirpPlanLen << endl;
	wcout << "\nPrinting kernel string...\n" << endl;
	plan1.Gen_Main_FFT_Kernel();
	string	tmpString = plan1.Get_Main_FFT_Kernel();
	char * cstr_ = new char [tmpString.length()+1];
	strcpy(cstr_, tmpString.c_str());
	wcout << cstr_ << endl;
	delete [] cstr_;
	return 0;
}
