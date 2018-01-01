#include <cstring>
#include <string>
#include <cstdint>
#include <iostream>

#include "header.h"

#define LOCAL_GRP_SZ 256

int main() {
	oclFFTHeader printer;
	std::string tmp;
	uint32_t local_fft_len = 7;
	uint32_t full_fft_len = 11 * local_fft_len;
	printer.SetFFTLength(full_fft_len);
	printer.SetKernelFFTLength(local_fft_len);
	printer.SetInitSCount(0);
	printer.SetBatchCount(0);
	printer.SetOuterCount(1);
	printer.SetInnerCount(3);
	printer.SetInitKW(full_fft_len);
	printer.SetInverseDirection(true);
	printer.SetPrecision(true);
	printer.SetIOFSScaleFactor(full_fft_len);
	printer.SetLocalGroupSize(LOCAL_GRP_SZ);
	printer.SetGlobalGroupSize(LOCAL_GRP_SZ);
	printer.SetMatrixInverseSCount(1+1);
	printer.SetTwiddleSCount(0);
	printer.SetFinalMatrixInverseSCount(1);
	tmp = printer.ReturnCommonMacros();
	tmp += printer.GetSourceHeader();
	tmp += printer.ReturnMacro();
	tmp.append(printer.print_kernel_name("phase1"));
	tmp.append(printer.print_kernel_float2_inputs());
	tmp.append(" {\n");
	tmp.append(printer.print_kernel_initialization());
	tmp.append(printer.print_kernel_outer_loop());
	tmp.append("\n}\n");
	char * cstr_ = new char [tmp.length()+1];
	std::strcpy(cstr_, tmp.c_str());
	std::wcout << cstr_ << std::endl;
	delete [] cstr_;
	return 0;
}
