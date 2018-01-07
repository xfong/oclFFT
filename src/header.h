#include <cstring>
#include <string>
#include <cstdint>
#include <iostream>
#include <ctime>

class oclFFTHeader {
	private:
		std::string preprocess_define = "#define";
		std::string preprocess_ifdef = "#ifdef";
		std::string preprocess_ifndef = "#ifndef";
		std::string preprocess_endif = "#endif";

		std::string kernel_global_input_float_type = "__global float";
		std::string kernel_global_input_double_type = "__global double";

		std::string def_2s;
		std::string def_3s;
		std::string def_4s;
		std::string def_5s;
		std::string def_7s;
		std::string def_8s;

		std::string def_2macro;
		std::string def_3rev;
		std::string def_4rev;
		std::string def_5rev;
		std::string def_7rev;
		std::string def_8rev;
		std::string def_3for;
		std::string def_4for;
		std::string def_5for;
		std::string def_7for;
		std::string def_8for;
		std::string def_twiddle_macro;
		std::string def_fft_sum_function;
		std::string def_fft_diff_function;
		std::string def_two_sum_function;
		std::string def_two_diff_function;

		uint32_t batch_count = 1;
		int sCount_init;
		uint32_t outer_loop_count;
		uint32_t inner_loop_count;
		uint32_t tps;
		uint32_t full_fft_length;
		uint32_t kernel_init_kw = 1;
		bool full_fft_length_set = false;
		uint32_t kernel_fft_length;
		bool kernel_fft_length_set = false;
		bool nameGenInit = false;
		int mat_inverse_sCount;
		int final_mat_inverse_sCount;
		int twiddle_sCount;
		bool inverse_fft;
		uint32_t iofs_scaleFac;
		uint32_t global_grp_size;
		uint32_t local_grp_size;
		uint32_t fft_type;
		uint32_t total_fft_count;
		bool isfloat;
		bool isterminal;
		uint32_t fft_input_stride;

		std::string print_kernel_inner_loop();
		std::string print_fft2_register_initialization();
		std::string print_fft_register_initialization();
		std::string print_kernel_register_initialization(bool isfloat, uint32_t vLen, std::string reg_prefix, uint32_t count);
		std::string print_kernel_array_fetch(std::string arrName, std::string offsVar, uint32_t stride, uint32_t inLen);
		std::string print_kernel_internal_twiddle_mult();
		std::string print_shift_input_to_registers();
		std::string print_scale_input_registers();
		std::string print_fft_macro();
		std::string print_matrix_inverse_logic();
		std::string print_local_buffer_copy();
		std::string print_final_matrix_inverse_logic();
		void SetTPS(uint32_t x);
		void SetTwiddleMacro();
		void SetTwoDiffFunction();
		void SetTwoSumFunction();
		void SetFFTSumFunction();
		void SetFFTDiffFunction();
		void GenFFT2Macros();
		void GenFFT3Macros();
		void GenFFT4Macros();
		void GenFFT5Macros();
		void GenFFT7Macros();
		void GenFFT8Macros();
		void initKernelNameGenerator();
		
	public:
		oclFFTHeader();
		std::string genKernelName(uint32_t nLen);
		std::string print_kernel_input_vector_type(bool isfloat, uint32_t vlen, std::string vName);
		std::string print_kernel_name(std::string kName);
		std::string print_kernel_float2_inputs();
		std::string print_kernel_initialization();
		std::string print_kernel_outer_loop();
		std::string GetSourceHeader();
		std::string ReturnCommonMacros();
		std::string ReturnMacro();

		void SetIOFSScaleFactor(uint32_t x) { iofs_scaleFac = x; };
		void SetGlobalGroupSize(uint32_t x) { global_grp_size = x; };
		void SetLocalGroupSize(uint32_t x) { local_grp_size = x; };
		void SetInitSCount(int x) { sCount_init = x; };
		void SetTwiddleSCount(int x) { twiddle_sCount = x; };
		void SetOuterCount(uint32_t x) { outer_loop_count = x; };
		void SetInnerCount(uint32_t x) { inner_loop_count = x; };
		void SetBatchCount(uint32_t x) { batch_count = x; };
		void SetFFTInputStride(uint32_t x) { fft_input_stride = x; };
		void SetMatrixInverseSCount(int x) { mat_inverse_sCount = x; };
		void SetFinalMatrixInverseSCount(int x) { final_mat_inverse_sCount = x; };
		void SetInverseDirection(bool x) { inverse_fft = x; };
		void SetFFTType(uint32_t x) { fft_type = x; };
		void SetTotalFFTCount(uint32_t x) { total_fft_count = x; };
		void SetPrecision(bool x);
		void SetInitKW(uint32_t x) { kernel_init_kw = x; };

		void SetFFTLength(uint32_t x);
		void SetKernelFFTLength(uint32_t x);
};
