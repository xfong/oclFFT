#include <cstring>
#include <string>
#include <cstdint>
#include <iostream>
#include <ctime>

#include "kernel_header.h"

oclFFTHeader::oclFFTHeader() {
	def_2s = preprocess_ifndef + " M_SCALE_1_2F\n\t#define M_SCALE_1_2F 0.500000000000000000000000f\n" + preprocess_endif + "\n\n";

	def_3s = def_2s + preprocess_ifndef + " M_SCALE_1_3F\n\t#define M_SCALE_1_3F 0.333333333333333333333333f\n" + preprocess_endif + "\n";
	def_3s += "\n" + preprocess_ifndef + " M_SQRT_3F\n\t#define M_SQRT_3F 1.73205080756887729352f\n" + preprocess_endif + "\n\n";

	def_4s = preprocess_ifndef + " M_SCALE_1_4F\n\t#define M_SCALE_1_4F 0.250000000000000000000000f\n" + preprocess_endif + "\n\n";

	def_5s = preprocess_ifndef + " M_SCALE_1_5F\n\t#define M_SCALE_1_5F 0.20000000000000000000f\n" + preprocess_endif + "\n";
	def_5s += "\n" + preprocess_ifndef + " M_COS_2PI_5F\n\t#define M_COS_2PI_5F 0.309016994374947451262869f\n" + preprocess_endif + "\n";
	def_5s += "\n" + preprocess_ifndef + " M_SIN_2PI_5F\n\t#define M_SIN_2PI_5F 0.951056516295153531181938f\n" + preprocess_endif + "\n";
	def_5s += "\n" + preprocess_ifndef + " M_COS_4PI_5F\n\t#define M_COS_4PI_5F -0.809016994374947451262869f\n" + preprocess_endif + "\n";
	def_5s += "\n" + preprocess_ifndef + " M_SIN_4PI_5F\n\t#define M_SIN_4PI_5F 0.587785252292473248125759f\n" + preprocess_endif + "\n\n";

	def_7s = preprocess_ifndef + " M_SCALE_1_7F\n\t#define M_SCALE_1_7F 0.142857142857142849212693f\n" + preprocess_endif + "\n";
	def_7s += "\n" + preprocess_ifndef + " M_COS_2PI_7F\n\t#define M_COS_2PI_7F 0.623489801858733594386308f\n" + preprocess_endif + "\n";
	def_7s += "\n" + preprocess_ifndef + " M_SIN_2PI_7F\n\t#define M_SIN_2PI_7F 0.781831482468029803634124f\n" + preprocess_endif + "\n";
	def_7s += "\n" + preprocess_ifndef + " M_COS_4PI_7F\n\t#define M_COS_4PI_7F -0.222520933956314337365257f\n" + preprocess_endif + "\n";
	def_7s += "\n" + preprocess_ifndef + " M_SIN_4PI_7F\n\t#define M_SIN_4PI_7F 0.974927912181823730364272f\n" + preprocess_endif + "\n";
	def_7s += "\n" + preprocess_ifndef + " M_COS_6PI_7F\n\t#define M_COS_6PI_7F -0.900968867902419034976447f\n" + preprocess_endif + "\n";
	def_7s += "\n" + preprocess_ifndef + " M_SIN_6PI_7F\n\t#define M_SIN_6PI_7F 0.433883739117558231423999f\n" + preprocess_endif + "\n\n";

	def_8s = preprocess_ifndef + " M_SCALE_1_8F\n\t#define M_SCALE_1_8F 0.125000000000000000000000f\n" + preprocess_endif + "\n\n";

	SetTwiddleMacro();
	SetTwoDiffFunction();
	SetTwoSumFunction();
	SetFFTSumFunction();
	SetFFTDiffFunction();
	GenFFT2Macros();
	GenFFT3Macros();
	GenFFT4Macros();
	GenFFT5Macros();
	GenFFT7Macros();
	GenFFT8Macros();
}

void oclFFTHeader::SetTwiddleMacro() {
	def_twiddle_macro = "#define twiddle_factor(k, angle, in) {\\\n\t";
	if (isfloat) {
		def_twiddle_macro += "float2 tw, v;\\\n\ttw.x = cos((float)k*angle);\\\n\ttw.y = sin((float)k*angle);\\\n\t";
	} else {
		def_twiddle_macro += "double2 tw, v;\\\n\ttw.x = cos((double)k*angle);\\\n\ttw.y = sin((double)k*angle);\\\n\t";
	}
	def_twiddle_macro += "v.x = fma(tw.x, in.x, 0.0f);\\\n\tv.x = fma(-1.0f * tw.y, in.y, v.x);\\\n\tv.y = fma(tw.x, in.y, 0.0f);\\\n\tv.y = fma(tw.y, in.x, v.y);\\\n\tin = v;\\\n}\n\n";
}

void oclFFTHeader::SetPrecision(bool x) {
	if (isfloat != x) {
		isfloat = x;
		SetTwiddleMacro();
		SetFFTDiffFunction();
		SetFFTSumFunction();
		SetTwoDiffFunction();
		SetTwoSumFunction();
		GenFFT2Macros();
		GenFFT3Macros();
		GenFFT4Macros();
		GenFFT5Macros();
		GenFFT7Macros();
		GenFFT8Macros();
	}
}

void oclFFTHeader::SetFFTSumFunction() {
	std::string precision_string;
	if (isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	def_fft_sum_function = "inline " + precision_string + "2 fft_sum(" + precision_string + "2 a, " + precision_string + "2 b) {\n\t";
	def_fft_sum_function += precision_string + "2 s;\n\t";
	def_fft_sum_function += "if (fabs(a.x) > fabs(b.x)) {\n\t\t";
	def_fft_sum_function += "s.x = fma(1.0f, a.x, b.x);\n\t";
	def_fft_sum_function += "} else {\n\t\t";
	def_fft_sum_function += "s.x = fma(1.0f, b.x, a.x);\n\t";
	def_fft_sum_function += "}\n\t";
	def_fft_sum_function += "if (fabs(a.y) > fabs(b.y)) {\n\t\t";
	def_fft_sum_function += "s.y = fma(1.0f, a.y, b.y);\n\t";
	def_fft_sum_function += "} else {\n\t\t";
	def_fft_sum_function += "s.y = fma(1.0f, b.y, a.y);\n\t";
	def_fft_sum_function += "}\n\t";
	def_fft_sum_function += "return s;\n";
	def_fft_sum_function += "}\n\n";
}

void oclFFTHeader::SetTwoSumFunction() {
	std::string precision_string;
	if (isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	def_two_sum_function = "inline " + precision_string + " two_sum(" + precision_string + " a, " + precision_string + " b) {\n\t";
	def_two_sum_function += precision_string + " s;\n";
	def_two_sum_function += "if (fabs(a) > fabs(b)) {\n\t\t";
	def_two_sum_function += "s = fma(1.0f, a, b);\n\t";
	def_two_sum_function += "} else {\n\t\t";
	def_two_sum_function += "s = fma(1.0f, b, a);\n\t";
	def_two_sum_function += "}\n\t";
	def_two_sum_function += "return s;\n";
	def_two_sum_function += "}\n\n";
}

void oclFFTHeader::SetFFTDiffFunction() {
	std::string precision_string;
	if (isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	def_fft_diff_function = "inline " + precision_string + "2 fft_diff(" + precision_string + "2 a, " + precision_string + "2 b) {\n\t";
	def_fft_diff_function += precision_string + "2 s;\n\t";
	def_fft_diff_function += "if (fabs(a.x) > fabs(b.x)) {\n\t\t";
	def_fft_diff_function += "s.x = fma(1.0f, a.x, -1.0f * b.x);\n\t";
	def_fft_diff_function += "} else {\n\t\t";
	def_fft_diff_function += "s.x = fma(-1.0f, b.x, a.x);\n\t";
	def_fft_diff_function += "}\n\t";
	def_fft_diff_function += "if (fabs(a.y) > fabs(b.y)) {\n\t\t";
	def_fft_diff_function += "s.y = fma(1.0f, a.y, -1.0f * b.y);\n\t";
	def_fft_diff_function += "} else {\n\t\t";
	def_fft_diff_function += "s.y = fma(-1.0f, b.y, a.y);\n\t";
	def_fft_diff_function += "}\n\t";
	def_fft_diff_function += "return s;\n";
	def_fft_diff_function += "}\n\n";
}

void oclFFTHeader::SetTwoDiffFunction() {
	std::string precision_string;
	if (isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	def_two_diff_function = "inline " + precision_string + " two_diff(" + precision_string + " a, " + precision_string + " b) {\n\t";
	def_two_diff_function += precision_string + " s;\n\t";
	def_two_diff_function += "if (fabs(a) > fabs(b)) {\n\t\t";
	def_two_diff_function += "s = fma(1.0f, a, -1.0f * b);\n\t";
	def_two_diff_function += "} else {\n\t\t";
	def_two_diff_function += "s = fma(-1.0f, b, a);\n\t";
	def_two_diff_function += "}\n\t";
	def_two_diff_function += "return s;\n";
	def_two_diff_function += "}\n\n";
}

std::string oclFFTHeader::print_kernel_input_vector_type(bool isfloat, uint32_t vlen, std::string vName) {
	std::string tmp;
	if (isfloat) {
		tmp = kernel_global_input_float_type;
	} else {
		tmp = kernel_global_input_double_type;
	}
	switch (vlen) {
		case 1:
			tmp.append("* ");
			break;
		case 2:
			tmp.append("2* ");
			break;
		case 3:
			tmp.append("3* ");
			break;
		case 4:
			tmp.append("4* ");
			break;
		default:
			std::wcout << "Cannot determine length of vlen!" << std::endl;
			return "-1";
	}
	tmp.append(vName);
	return tmp;
}

std::string oclFFTHeader::print_kernel_name(std::string kName) {
	std::string tmp ("__kernel void\n");
	tmp.append(kName);
	return tmp;
}

std::string oclFFTHeader::print_kernel_register_initialization(bool isfloat, uint32_t vLen, std::string reg_prefix, uint32_t count) {
	std::string tmp = "";
	for (uint32_t idx = 0; idx < count; idx++) {
		std::string curr_def;
		if (isfloat) {
			curr_def = "\t\t\tfloat";
		} else {
			curr_def = "\t\t\tdouble";
		}
		switch (vLen) {
			case 1:
				break;
			case 2:
			case 3:
			case 4:
				curr_def.append(std::to_string(vLen));
				break;
			default:
				std::wcout << "Cannot determine length of vlen!" << std::endl;
				return "-1";
		}
		tmp.append(curr_def);
		tmp.append(" " + reg_prefix + std::to_string(idx) + ";\n");
	}
	return tmp;
}

std::string oclFFTHeader::print_kernel_array_fetch(std::string arrName, std::string offsVar, uint32_t stride, uint32_t inLen) {
	std::string tmp;
	if (inLen > 0) {
		for (uint32_t idx = 0; idx < inLen; idx++) {
			if (idx == 0) {
				tmp.append("in" + std::to_string(idx) + " = dataIn[" + offsVar + "];\n");
			} else {
				tmp.append("in" + std::to_string(idx) + " = dataIn[" + offsVar + " + " + std::to_string(idx*stride) + "];\n");
			}
		}
	}
	return tmp;
}

std::string oclFFTHeader::print_kernel_float2_inputs() {
	std::string tmp = "(";
	tmp.append(print_kernel_input_vector_type(true, 2, "dataIn"));
	tmp.append(", ");
	tmp.append(print_kernel_input_vector_type(true, 2, "buffer"));
	tmp.append(")");
	return tmp;
}

std::string oclFFTHeader::print_kernel_initialization() {
	std::string tmp;
	tmp = "\t// Initialize indices for kernel\n";
	tmp += "\tint lid = get_local_id(0);\n";
	tmp += "\tint itr = get_group_id(0)*" + std::to_string(local_grp_size) + " + lid;\n";
	tmp += "\tint sCount = " + std::to_string(sCount_init) + ";\n";
	if (inverse_fft) {
		tmp += "\tfloat angle_ = 2.0f * M_PI_F / (float)(" + std::to_string(kernel_init_kw) + ");\n";
	} else {
		tmp += "\tfloat angle_ = -2.0f * M_PI_F / (float)(" + std::to_string(kernel_init_kw) + ");\n";
	}
	return tmp;
}

std::string oclFFTHeader::print_kernel_outer_loop() {
	std::string tmp;
	tmp = "\n\t// Outer loop for each work item to emulate many work items\n";
	tmp += "\twhile (sCount < " + std::to_string(outer_loop_count) + ") {";
	tmp += "\n\t\tint idx = itr;\n";
	tmp += "\t\tint lti = idx % " + std::to_string(tps) + ";\n";
	tmp += "\t\tint iofs = idx - lti;\n";
	tmp += "\t\tiofs /= " + std::to_string(tps) + ";\n";
	tmp += "\t\tiofs *= " + std::to_string(iofs_scaleFac) + ";\n";
	tmp += print_kernel_inner_loop();
	tmp += "\n\t\tsCount++;\n";
	tmp.append(print_final_matrix_inverse_logic());
	tmp += "\t}";
	return tmp;
}

std::string oclFFTHeader::print_kernel_inner_loop() {
	std::string tmp;
	tmp = "\n\t\t// Inner loop for each work item to emulate many work items\n";
	tmp += "\t\twhile (idx < " + std::to_string(inner_loop_count) + ") {";
	switch (kernel_fft_length) {
		case 2:
			tmp.append(print_fft2_register_initialization());
			break;
		default:
			tmp.append(print_fft_register_initialization());
			break;
	};
	tmp += print_shift_input_to_registers();
	tmp += "\n";
	if (inverse_fft) {
		tmp += "\t\t\t// Scaling input values for inverse FFT\n";
		tmp.append(print_scale_input_registers());
		tmp += "\n";
	}
	tmp += print_kernel_internal_twiddle_mult();
	tmp += "\n" + print_fft_macro();
	tmp += print_matrix_inverse_logic();
	tmp += "\n\t\t\tidx += " + std::to_string(global_grp_size) + ";\n";
	tmp += "\t\t\tangle_ = fma(" + std::to_string(kernel_fft_length) + ".0f, angle_, 0.0f);\n\t\t}";
	return tmp;
}

std::string oclFFTHeader::print_fft2_register_initialization() {
	std::string tmp = "\n\t\t\t// Initializing local registers for length-2 FFT";
	tmp.append("\n" + print_kernel_register_initialization(isfloat, 2, "in", 2));
	tmp.append("\n");
	tmp.append("\n" + print_kernel_register_initialization(isfloat, 2, "v", 1));
	return tmp;
}

std::string oclFFTHeader::print_fft_register_initialization() {
	std::string tmp = "\n\t\t\t// Initializing local registers for FFT";
	tmp.append("\n" + print_kernel_register_initialization(true, 2, "in", kernel_fft_length));
	if (kernel_fft_length == 3) {
		tmp.append("\n" + print_kernel_register_initialization(true, 2, "v", kernel_fft_length+1));
	}
	if ((kernel_fft_length == 4) || (kernel_fft_length == 8)) {
		tmp.append("\n" + print_kernel_register_initialization(true, 2, "v", kernel_fft_length));
	}
	if ((kernel_fft_length == 5) || (kernel_fft_length == 7)) {
		tmp.append("\n" + print_kernel_register_initialization(true, 2, "v", kernel_fft_length+2));
	}
	if ((kernel_fft_length != 3) && (kernel_fft_length != 4) && (kernel_fft_length != 5) && (kernel_fft_length != 7) && (kernel_fft_length != 8)) {
		tmp.append("Unsupported FFT Length detected!");
	}
	tmp.append("\n");
	return tmp;
}

std::string oclFFTHeader::print_kernel_internal_twiddle_mult() {
	std::string tmp = "\t\t\t// Performing twiddle factor multiplication before FFT\n";
	tmp += "\t\t\tif (sCount > ";
	tmp.append(std::to_string(twiddle_sCount) + ") {\n");
	tmp.append("\t\t\t\tfloat angle = fma(angle_, (float)(lti), 0.0f);\n");
	for (uint32_t idx = 1; idx < kernel_fft_length; idx++) {
		tmp.append("\t\t\t\ttwiddle_factor(" + std::to_string(idx) + ", angle, in" + std::to_string(idx) + ");\n");
	}
	tmp.append("\t\t\t}\n");
	return tmp;
}

std::string oclFFTHeader::print_shift_input_to_registers() {
	std::string tmp = "";
	std::string tmp1 = "";
	for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
		tmp.append("\t\t\tint Idout" + std::to_string(idx) + " = idx + iofs");
		if (idx > 0) {
			tmp.append(" + " + std::to_string(idx*fft_input_stride));
		}
		tmp.append(";\n");
		tmp1.append("\t\t\tin" + std::to_string(idx) + " = dataIn[Idout" + std::to_string(idx) + "];\n");
	}
	return tmp + "\n" + tmp1;
}

std::string oclFFTHeader::print_scale_input_registers() {
	std::string tmp = "";
	std::string tmp1 = "M_SCALE_1_" + std::to_string(kernel_fft_length) + "F";
	for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
		tmp.append("\t\t\tin" + std::to_string(idx) + " = fma(" + tmp1 + ", in" + std::to_string(idx) + ", 0.0f);\n");
	}
	return tmp;
}

std::string oclFFTHeader::print_fft_macro() {
	std::string tmp = "\t\t\t// Call FFT macro\n\t\t\t";
	if (kernel_fft_length == 2) {
		tmp += "FFT2(in0, in1);\n";
	} else {
		if (inverse_fft) {
			tmp += "I";
		}
		tmp.append("FFT" + std::to_string(kernel_fft_length) + "(");
		for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
			tmp.append("in" + std::to_string(idx));
			if (idx < kernel_fft_length - 1) {
				tmp.append(", ");
			}
		}
		tmp.append(");\n");
	}
	return tmp;
}

std::string oclFFTHeader::print_matrix_inverse_logic() {
	std::string tmp = "\n\t\t\t// Printing logic to perform matrix inverse\n";
	tmp.append("\t\t\tif (sCount == " + std::to_string(mat_inverse_sCount) + ") {\n");
	for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
		tmp.append("\t\t\t\tdataIn[Idout" + std::to_string(idx) + "] = in" + std::to_string(idx) + ";\n");
	}
	tmp.append("\t\t\t} else {\n" + print_local_buffer_copy());
	return tmp + "\n\t\t\t}\n";
}

std::string oclFFTHeader::print_local_buffer_copy() {
	std::string tmp = "\t\t\t\t";
	if (isfloat) {
		tmp.append("__local float2 block_cache[");
	} else {
		tmp.append("__local double2 block_cache[");
	}
	tmp += std::to_string(local_grp_size*kernel_fft_length*2) + "];\n";
	for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
		std::string tmp1 = "\t\t\t\t";
		if (idx == 0) {
			tmp1.append("block_cache[" + std::to_string(kernel_fft_length) + "*lid] = in0;\n");
		} else {
			tmp1.append("block_cache[" + std::to_string(kernel_fft_length) + "*lid + " + std::to_string(idx) + "] = in" + std::to_string(idx) + ";\n");
		}
		tmp.append(tmp1);
	}
	tmp += "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
	for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
		std::string tmp1 = "\t\t\t\tin" + std::to_string(idx) + " = block_cache[lid";
		if (idx == 0) {
			tmp1.append("];\n");
		} else {
			tmp1.append(" + " + std::to_string(idx*local_grp_size) + "];\n");
		}
		tmp.append(tmp1);
	}
	tmp += "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
	for (uint32_t idx = 0; idx < kernel_fft_length; idx++) {
		if (idx == 0) {
			tmp.append("\t\t\t\tbuffer[Idout0] = in0;\n");
		} else {
			tmp.append("\t\t\t\tbuffer[Idout0 + " + std::to_string(idx*local_grp_size) + "] = in" + std::to_string(idx) + ";\n");
		}
	}
	tmp.append("\t\t\t\tbarrier(CLK_GLOBAL_MEM_FENCE);");
	return tmp;
}

std::string oclFFTHeader::print_final_matrix_inverse_logic() {
	std::string tmp = "\t\t// Final stage matrix transpose logic\n\t\tif (sCount < ";
	tmp.append(std::to_string(final_mat_inverse_sCount) + ") {\n\t\t\tidx = itr;\n");
	tmp.append("\t\t\twhile (idx < " + std::to_string(batch_count*full_fft_length));
	tmp.append(") {\n\t\t\t\t");
	if (isfloat) {
		tmp.append("float2 buf_val = buffer[idx];\n");
	} else {
		tmp.append("double2 buf_val = buffer[idx];\n");
	}
	tmp.append("\t\t\t\tdata[idx] = buf_val;\n");
	tmp.append("\t\t\t\tidx += " + std::to_string(local_grp_size) + ";\n");
	tmp.append("\t\t\t}\n\t\t}\n");
	return tmp;
}

void oclFFTHeader::SetFFTLength(uint32_t x) {
	full_fft_length = x;
	full_fft_length_set = true;
	if (full_fft_length_set && kernel_fft_length_set) {
		SetTPS(full_fft_length / kernel_fft_length);
	}
}

void oclFFTHeader::SetKernelFFTLength(uint32_t x) {
	kernel_fft_length = x;
	kernel_fft_length_set = true;
	if (full_fft_length_set && kernel_fft_length_set) {
		SetTPS(full_fft_length / kernel_fft_length);
	}
}

std::string oclFFTHeader::GetSourceHeader() {
	switch (kernel_fft_length) {
		case 8:
			return def_8s;
		case 7:
			return def_7s;
		case 5:
			return def_5s;
		case 4:
			return def_4s;
		case 3:
			return def_3s;
		case 2:
			return def_2s;
		default:
			return "";
	}
}

void oclFFTHeader::SetTPS(uint32_t x) {
	tps = x;
	SetFFTInputStride(tps);
}

std::string oclFFTHeader::ReturnCommonMacros() {
	std::string tmp = "";
	tmp += def_twiddle_macro + def_fft_diff_function + def_fft_sum_function;
	return tmp;
}

std::string oclFFTHeader::ReturnMacro() {
	std::string macro_string = "";
	switch (kernel_fft_length) {
		case 8:
			if (inverse_fft) {
				macro_string = def_8rev;
			} else {
				macro_string = def_8for;
			}
			macro_string += "\n" + def_two_diff_function + def_two_sum_function;
			break;
		case 7:
			if (inverse_fft) {
				macro_string = def_7rev;
			} else {
				macro_string = def_7for;
			}
			macro_string += "\n";
			break;
		case 5:
			if (inverse_fft) {
				macro_string = def_5rev;
			} else {
				macro_string = def_5for;
			}
			macro_string += "\n";
			break;
		case 4:
			if (inverse_fft) {
				macro_string = def_4rev;
			} else {
				macro_string = def_4for;
			}
			macro_string += "\n";
			break;
		case 3:
			if (inverse_fft) {
				macro_string = def_3rev;
			} else {
				macro_string = def_3for;
			}
			macro_string += "\n";
			break;
		case 2:
			macro_string = def_2macro + "\n";
			break;
		default:
			break;
	}
	return macro_string;
}

void oclFFTHeader::GenFFT2Macros() {
	def_2macro = "#define FFT2(in0, in1) {\\\n\tv0 = in0;\\\n\tin0 = fft_sum(v0, in1);\\\n\tin1 = fft_diff(v0, in1);\\\n}\n";
}

void oclFFTHeader::GenFFT3Macros() {
	def_3for = "#define FFT3(in0, in1, in2) {\\\n\tv0 = in0;\\\n\tv1 = fft_sum(in1, in2);\\\n\tin0 = fft_sum(in0, v1);\\\n\t";
	def_3for += "v1 = fma(M_SCALE_1_2F, v1, 0.0f);\\\n\tv2 = fft_diff(in2, in1);\\\n\tv3 = fft_diff(in1, in2);\\\n\t";
	def_3for += "in1 = fft_diff(v0, v1);\\\n\tin2 = in1;\\\n\tv3.x = v3.y;\\\n\tv3.y = v2.x;\\\n\tv3 = fma(M_SCALE_1_2F, v3, 0.0f);\\";
	def_3for += "\n\tv3 = fma(M_SQRT_3F, v3, 0.0f);\\\n\tv2 = in1;\\\n\tin1 = fft_sum(v2, v3);\\\n\tin2 = fft_diff(v2, v3);\\\n}\n";

	def_3rev = "#define IFFT3(in0, in1, in2) {\\\n\tv0 = in0;\\\n\tv1 = fft_sum(in1, in2);\\\n\tin0 = fft_sum(in0, v1);\\\n\tv1 = fma(M_SCALE_1_2F, v1, 0.0f);\\";
	def_3rev += "\n\tv2 = fft_diff(in2, in1);\\\n\tv3 = fft_diff(in1, in2);\\\n\tin1 = fft_diff(v0, v1);\\\n\tin2 = in1;\\\n\tv3.x = v3.y;\\";
	def_3rev += "\n\tv3.y = v2.x;\\\n\tv3 = fma(M_SCALE_1_2F, v3, 0.0f);\\\n\tv3 = fma(M_SQRT_3F, v3, 0.0f);\\\n\tv2 = in1;\\\n\tin1 = fft_diff(v2, v3);\\\n\tin2 = fft_sum(v2, v3);\\\n}\n";
}

void oclFFTHeader::GenFFT4Macros() {
	def_4for = "#define FFT4(in0, in1, in2, in3) {\\\n\tv0 = fft_sum(in0, in2);\\\n\tv1 = fft_diff(in3, in1);\\\n\t";
	def_4for += "v2 = fft_diff(in1, in3);\\\n\tv3.x = v2.y;\\\n\tv3.y = v1.x;\\\n\tv1 = fft_sum(in1, in3);\\\n\tv2 = fft_diff(in0, in2);\\\n\t";
	def_4for += "in0 = fft_sum(v0, v1);\\\n\tin2 = fft_diff(v0, v1);\\\n\tin1 = fft_sum(v2, v3);\\\n\tin3 = fft_diff(v2, v3);\\\n}\n";

	def_4rev = "#define IFFT4(in0, in1, in2, in3) {\\\n\tv0 = fft_sum(in0, in2);\\\n\tv1 = fft_diff(in3, in1);\\\n\t";
	def_4rev += "v2 = fft_diff(in1, in3);\\\n\tv3.x = v1.y;\\\n\tv3.y = v2.x;\\\n\tv1 = fft_sum(in1, in3);\\\n\tv2 = fft_diff(in0, in2);\\\n\t";
	def_4rev += "in0 = fft_sum(v0, v1);\\\n\tin2 = fft_diff(v0, v1);\\\n\tin1 = fft_sum(v2, v3);\\\n\tin3 = fft_diff(v2, v3);\\\n}\n";
}

void oclFFTHeader::GenFFT5Macros() {
	def_5for = "#define FFT5(in0, in1, in2, in3, in4) {\\\n\tv0 = in0;\\\n\tv1 = fft_sum(in1, in4);\\\n\t";
	def_5for += "v2 = fft_sum(in2, in3);\\\n\tv3 = fft_sum(v1, v2);\\\n\tin0 = fft_sum(v0, v3);\\\n\tv3 = fft_diff(in1, in4);\\\n\tv4 = fft_diff(in2, in3);\\\n\t";
	def_5for += "in1 = fma(M_COS_2PI_5F, v1, 0.0f);\\\n\tin2 = fma(M_COS_4PI_5F, v1, 0.0f);\\\n\tin3 = fma(M_COS_2PI_5F, v2, 0.0f);\\\n\tin4 = fma(M_COS_4PI_5F, v2, 0.0f);\\\n\tv1 = fft_sum(in1, in4);\\\n\t";
	def_5for += "v2 = fft_sum(in2, in3);\\\n\tin1 = fma(M_SIN_2PI_5F, v3, 0.0f);\\\n\tin2 = fma(M_SIN_4PI_5F, v3, 0.0f);\\\n\tin3 = fma(M_SIN_2PI_5F, v4, 0.0f);\\\n\tin4 = fma(M_SIN_4PI_5F, v4, 0.0f);\\\n\t";
	def_5for += "v3 = fft_sum(in1, in4);\\\n\tv4 = fft_diff(in2, in3);\\\n\tv5.x = v4.y;\\\n\tv5.y = -1.0f * v4.x;\\\n\tv6.x = v3.y;\\\n\t";
	def_5for += "v6.y = fma(-1.0f, v3.x, 0.0f);\\\n\tin1 = fft_sum(v0, v1);\\\n\tin2 = fft_sum(v0, v2);\\\n\tv1 = in1;\\\n\tv2 = in2;\\\n\t";
	def_5for += "in1 = fft_sum(v1, v6);\\\n\tin2 = fft_sum(v2, v5);\\\n\tin3 = fft_diff(v2, v5);\\\n\tin4 = fft_diff(v1, v6);\\\n}\n";

	def_5rev = "#define IFFT5(in0, in1, in2, in3, in4) {\\\n\tv0 = in0;\\\n\tv1 = fft_sum(in1, in4);\\\n\tv2 = fft_sum(in2, in3);\\\n\t";
	def_5rev += "v3 = fft_sum(v1, v2);\\\n\tin0 = fft_sum(v0, v3);\\\n\tv3 = fft_diff(in1, in4);\\\n\tv4 = fft_diff(in2, in3);\\\n\t";
	def_5rev += "in1 = fma(M_COS_2PI_5F, v1, 0.0f);\\\n\tin2 = fma(M_COS_4PI_5F, v1, 0.0f);\\\n\tin3 = fma(M_COS_2PI_5, v2, 0.0f);\\\n\tin4 = fma(M_COS_4PI_5F, v2, 0.0f);\\\n\t";
	def_5rev += "v1 = fft_sum(in1, in4);\\\n\tv2 = fft_sum(in2, in3);\\\n\tin1 = fma(M_SIN_2PI_5F, v3, 0.0f);\\\n\tin2 = fma(M_SIN_4PI_5F, v3, 0.0f);\\\n\t";
	def_5rev += "in3 = fma(M_SIN_2PI_5F, v4, 0.0f);\\\n\tin4 = fma(M_SIN_4PI_5F, v4, 0.0f);\\\n\tv3 = fft_sum(in1, in4);\\\n\tv4 = fft_diff(in2, in3);\\\n\t";
	def_5rev += "v5.x = fma(-1.0f, v4.y, 0.0f);\\\n\tv5.y = v4.x;\\\n\tv6.x = fma(-1.0f, v3.y, 0.0f);\\\n\tv6.y = v3.x;\\\n\tin1 = fft_sum(v0, v1);\\\n\t";
	def_5rev += "in2 = fft_sum(v0, v2);\\\n\tv1 = in1;\\\n\tv2 = in2;\\\n\tin1 = fft_sum(v1, v6);\\\n\tin2 = fft_sum(v2, v5);\\\n\t";
	def_5rev += "in3 = fft_diff(v2, v5);\\\n\tin4 = fft_diff(v1, v6);\\\n}\n";
}

void oclFFTHeader::GenFFT7Macros() {
	def_7for = "#define FFT7(in0, in1, in2, in3, in4, in5, in6) {\\\n\tv0 = in0;\\\n\tv1 = fft_sum(in1, in6);\\\n\tv2 = fft_sum(in2, in5);\\\n\t";
	def_7for += "v3 = fft_sum(in3, in4);\\\n\tv7 = M_COS_2PI_7F*v1;\\\n\tv8 = M_COS_2PI_7F*v2;\\\n\tv9 = M_COS_2PI_7F*v3;\\\n\tv4 = fft_sum(v0, v1);\\\n\t";
	def_7for += "v5 = fft_sum(v2, v3);\\\n\tin0 = fft_sum(v4, v5);\\\n\tv4 = fft_diff(in6, in1);\\\n\tv5 = fft_diff(in5, in2);\\\n\t";
	def_7for += "v6 = fft_diff(in4, in3);\\\n\tin1 = fft_sum(v0, v7);\\\n\tin2 = fft_sum(v0, v9);\\\n\tin3 = fft_sum(v0, v8);\\\n\tv7 = M_COS_4PI_7F*v1;\\\n\t";
	def_7for += "v8 = M_COS_4PI_7F*v2;\\\n\tv9 = M_COS_4PI_7F*v3;\\\n\tin4 = fft_sum(in1, v8);\\\n\tin5 = fft_sum(in2, v7);\\\n\tin6 = fft_sum(in3, v9);\\\n\t";
	def_7for += "v7 = M_COS_6PI_7F*v1;\\\n\tv8 = M_COS_6PI_7F*v2;\\\n\tv9 = M_COS_6PI_7F*v3;\\\n\tin1 = fft_sum(in4, v9);\\\n\tin2 = fft_sum(in5, v8);\\\n\t";
	def_7for += "in3 = fft_sum(in6, v7);\\\n\tv1.y = v4.x;\\\n\tv1.x = -1.0f * v4.y;\\\n\tv2.y = v5.x;\\\n\tv2.x = -1.0f * v5.y;\\\n\t";
	def_7for += "v3.y = v6.x;\\\n\tv3.x = -1.0f * v6.y;\\\n\tv7 = M_SIN_2PI_7F*v1;\\\n\tv8 = M_SIN_2PI_7F*v2;\\\n\tv9 = M_SIN_2PI_7F*v3;\\\n\t";
	def_7for += "v4 = fft_sum(in1, v7);\\\n\tv5 = fft_diff(in2, v9);\\\n\tv6 = fft_diff(in3, v8);\\\n\tin4 = fft_sum(in3, v8);\\\n\t";
	def_7for += "in5 = fft_sum(in2, v9);\\\n\tin6 = fft_diff(in1, v7);\\\n\tv7 = M_SIN_4PI_7F*v1;\\\n\tv8 = M_SIN_4PI_7F*v2;\\\n\tv9 = M_SIN_4PI_7F*v3;\\\n\t";
	def_7for += "in1 = fft_sum(v4, v8);\\\n\tin2 = fft_sum(v5, v7);\\\n\tin3 = fft_sum(v6, v9);\\\n\tv4 = fft_diff(in4, v9);\\\n\t";
	def_7for += "v5 = fft_diff(in5, v7);\\\n\tv6 = fft_diff(in6, v8);\\\n\tv7 = M_SIN_6PI_7F*v1;\\\n\tv8 = M_SIN_6PI_7F*v2;\\\n\tv9 = M_SIN_6PI_7F*v3;\\\n\tv1 = in1;\\\n\t";
	def_7for += "v2 = in2;\\\n\tv3 = in3;\\\n\tin1 = fft_sum(v1, v9);\\\n\tin2 = fft_diff(v2, v8);\\\n\tin3 = fft_sum(v3, v7);\\\n\tin4 = fft_diff(v4, v7);\\\n\t";
	def_7for += "in5 = fft_sum(v5, v8);\\\n\tin6 = fft_diff(v6, v9);\\\n}\n";

	def_7rev = "#define IFFT7(in0, in1, in2, in3, in4, in5, in6) {\\\n\tv0 = in0;\\\n\tv1 = fft_sum(in1, in6);\\\n\tv2 = fft_sum(in2, in5);\\\n\tv3 = fft_sum(in3, in4);\\\n\t";
	def_7rev += "v7 = M_COS_2PI_7F*v1;\\\n\tv8 = M_COS_2PI_7F*v2;\\\n\tv9 = M_COS_2PI_7F*v3;\\\n\tv4 = fft_sum(v0, v1);\\\n\tv5 = fft_sum(v2, v3);\\\n\t";
	def_7rev += "in0 = fft_sum(v4, v5);\\\n\tv4 = fft_diff(in1, in6);\\\n\tv5 = fft_diff(in2, in5);\\\n\tv6 = fft_diff(in3, in4);\\\n\t";
	def_7rev += "in1 = fft_sum(v0, v7);\\\n\tin2 = fft_sum(v0, v9);\\\n\tin3 = fft_sum(v0, v8);\\\n\tv7 = M_COS_4PI_7F*v1;\\\n\t";
	def_7rev += "v8 = M_COS_4PI_7F*v2;\\\n\tv9 = M_COS_4PI_7F*v3;\\\n\tin4 = fft_sum(in1, v8);\\\n\tin5 = fft_sum(in2, v7);\\\n\t";
	def_7rev += "in6 = fft_sum(in3, v9);\\\n\tv7 = M_COS_6PI_7F*v1;\\\n\tv8 = M_COS_6PI_7F*v2;\\\n\tv9 = M_COS_6PI_7F*v3;\\\n\t";
	def_7rev += "in1 = fft_sum(in4, v9);\\\n\tin2 = fft_sum(in5, v8);\\\n\tin3 = fft_sum(in6, v7);\\\n\tv1.y = v4.x;\\\n\t";
	def_7rev += "v1.x = -1.0f * v4.y;\\\n\tv2.y = v5.x;\\\n\tv2.x = -1.0f * v5.y;\\\n\tv3.y = v6.x;\\\n\t";
	def_7rev += "v3.x = -1.0f * v6.y;\\\n\tv7 = M_SIN_2PI_7F*v1;\\\n\tv8 = M_SIN_2PI_7F*v2;\\\n\tv9 = M_SIN_2PI_7F*v3;\\\n\t";
	def_7rev += "v4 = fft_sum(in1, v7);\\\n\tv5 = fft_diff(in2, v9);\\\n\tv6 = fft_diff(in3, v8);\\\n\tin4 = fft_sum(in3, v8);\\\n\t";
	def_7rev += "in5 = fft_sum(in2, v9);\\\n\tin6 = fft_diff(in1, v7);\\\n\tv7 = M_SIN_4PI_7F*v1;\\\n\tv8 = M_SIN_4PI_7F*v2;\\\n\t";
	def_7rev += "v9 = M_SIN_4PI_7F*v3;\\\n\tin1 = fft_sum(v4, v8);\\\n\tin2 = fft_sum(v5, v7);\\\n\tin3 = fft_sum(v6, v9);\\\n\t";
	def_7rev += "v4 = fft_diff(in4, v9);\\\n\tv5 = fft_diff(in5, v7);\\\n\tv6 = fft_diff(in6, v8);\\\n\tv7 = M_SIN_6PI_7F*v1;\\\n\tv8 = M_SIN_6PI_7F*v2;\\\n\t";
	def_7rev += "v9 = M_SIN_6PI_7F*v3;\\\n\tv1 = in1;\\\n\tv2 = in2;\\\n\tv3 = in3;\\\n\tin1 = fft_sum(v1, v9);\\\n\t";
	def_7rev += "in2 = fft_diff(v2, v8);\\\n\tin3 = fft_sum(v3, v7);\\\n\tin4 = fft_diff(v4, v7);\\\n\t";
	def_7rev += "in5 = fft_sum(v5, v8);\\\n\tin6 = fft_diff(v6, v9);\\\n}\n";
}

void oclFFTHeader::GenFFT8Macros() {
	def_8for = "#define FFT8(in0, in1, in2, in3, in4, in5, in6, in7) {\\\n\tv0 = fft_sum(in0, in4);\\\n\tv1 = fft_diff(in6, in2);\\\n\t";
	def_8for += "v2 = fft_diff(in2, in6);\\\n\tv3.x = v2.y;\\\n\tv3.y = v1.x;\\\n\tv1 = fft_sum(in2, in6);\\\n\t";
	def_8for += "v2 = fft_diff(in0, in4);\\\n\tin0 = fft_sum(v0, v1);\\\n\tin4 = fft_diff(v0, v1);\\\n\tin2 = fft_sum(v2, v3);\\\n\t";
	def_8for += "in6 = fft_diff(v2, v3);\\\n\tv0 = fft_sum(in1, in5);\\\n\tv1 = fft_diff(in7, in3);\\\n\tv2 = fft_diff(in3, in7);\\\n\t";
	def_8for += "v3.x = v2.y;\\\n\tv3.y = v1.x;\\\n\tv1 = fft_sum(in3, in7);\\\n\tv2 = fft_diff(in1, in5);\\\n\t";
	def_8for += "in1 = fft_sum(v0, v1);\\\n\tin5 = fft_diff(v0, v1);\\\n\tin3 = fft_sum(v2, v3);\\\n\tin7 = fft_diff(v2, v3);\\\n\t";
	def_8for += "v0 = in0;\\\n\tv1 = in1;\\\n\tv2 = in2;\\\n\tv3.x = two_sum(in3.x, in3.y);\\\n\t";
	def_8for += "v3.y = two_diff(in3.y, in3.x);\\\n\tv3.x *= M_SQRT1_2F;\\\n\tv3.y *= M_SQRT1_2F;\\\n\tv4 = in4;\\\n\t";
	def_8for += "v5.x = in5.y;\\\n\tv5.y = -1.0f * in5.x;\\\n\tv6 = in6;\\\n\tv7.x = two_diff(in7.y, in7.x);\\\n\t";
	def_8for += "v7.y = -1.0f * two_sum(in7.x, in7.y);\\\n\tv7.x *= M_SQRT1_2F;\\\n\tv7.y *= M_SQRT1_2F;\\\n\tin0 = fft_sum(v0, v1);\\\n\t";
	def_8for += "in1 = fft_sum(v2, v3);\\\n\tin2 = fft_sum(v4, v5);\\\n\tin3 = fft_sum(v6, v7);\\\n\tin4 = fft_diff(v0, v1);\\\n\t";
	def_8for += "in5 = fft_diff(v2, v3);\\\n\tin6 = fft_diff(v4, v5);\\\n\tin7 = fft_diff(v6, v7);\\\n}\n";

	def_8rev = "#define IFFT8(in0, in1, in2, in3, in4, in5, in6, in7) {\\\n\tv0 = fft_sum(in0, in4);\\\n\tv1 = fft_diff(in6, in2);\\\n\tv2 = fft_diff(in2, in6);\\\n\t";
	def_8rev += "v3.x = v1.y;\\\n\tv3.y = v2.x;\\\n\tv1 = fft_sum(in2, in6);\\\n\tv2 = fft_diff(in0, in4);\\\n\t";
	def_8rev += "in0 = fft_sum(v0, v1);\\\n\tin4 = fft_diff(v0, v1);\\\n\tin2 = fft_sum(v2, v3);\\\n\tin6 = fft_diff(v2, v3);\\\n\t";
	def_8rev += "v0 = fft_sum(in1, in5);\\\n\tv1 = fft_diff(in7, in3);\\\n\tv2 = fft_diff(in3, in7);\\\n\tv3.x = v1.y;\\\n\t";
	def_8rev += "v3.y = v2.x;\\\n\tv1 = fft_sum(in3, in7);\\\n\tv2 = fft_diff(in1, in5);\\\n\tin1 = fft_sum(v0, v1);\\\n\t";
	def_8rev += "in5 = fft_diff(v0, v1);\\\n\tin3 = fft_sum(v2, v3);\\\n\tin7 = fft_diff(v2, v3);\\\n\tv0 = in0;\\\n\t";
	def_8rev += "v1 = in1;\\\n\tv2 = in2;\\\n\tv3.x = two_diff(in3.x, in3.y);\\\n\tv3.y = two_sum(in3.x, in3.y);\\\n\t";
	def_8rev += "v3.x *= M_SQRT1_2F;\\\n\tv3.y *= M_SQRT1_2F;\\\n\tv4 = in4;\\\n\tv5.x = -1.0f * in5.y;\\\n\t";
	def_8rev += "v5.y = in5.x;\\\n\tv6 = in6;\\\n\tv7.x = -1.0f * two_sum(in7.x, in7.y);\\\n\tv7.y = two_diff(in7.x, in7.y);\\\n\t";
	def_8rev += "v7.x *= M_SQRT1_2F;\\\n\tv7.y *= M_SQRT1_2F;\\\n\tin0 = fft_sum(v0, v1);\\\n\tin1 = fft_sum(v2, v3);\\\n\t";
	def_8rev += "in2 = fft_sum(v4, v5);\\\n\tin3 = fft_sum(v6, v7);\\\n\tin4 = fft_diff(v0, v1);\\\n\tin5 = fft_diff(v2, v3);\\\n\t";
	def_8rev += "in6 = fft_diff(v4, v5);\\\n\tin7 = fft_diff(v6, v7);\\\n}\n";
}

void oclFFTHeader::initKernelNameGenerator() {
	srand(time(0));
}

std::string oclFFTHeader::genKernelName(uint32_t nLen) {
	static const char alphaNum[] = "aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ";
	int stringLength = sizeof(alphaNum) - 1;
	std::string tmp = "";
	if (!nameGenInit) {
		initKernelNameGenerator();
	}
	for (uint32_t idx = 0; idx < nLen; idx++) {
		tmp += alphaNum[rand() % stringLength];
	}
	return tmp;
}
