package oclFFT

const (
	preprocess_define = string("#define")
	preprocess_ifdef = string("#ifdef");
	preprocess_ifndef = string("#ifndef");
	preprocess_endif = string("#endif");

	kernel_global_input_half_type = string("__global half");
	kernel_global_input_float_type = string("__global float");
	kernel_global_input_double_type = string("__global double");
)

type kernelWriter struct {
	def_2s						string
	def_3s						string
	def_4s						string
	def_5s						string
	def_7s						string
	def_8s						string

	def_2macro					string
	def_3rev					string
	def_4rev					string
	def_5rev					string
	def_7rev					string
	def_8rev					string
	def_3for					string
	def_4for					string
	def_5for					string
	def_7for					string
	def_8for					string

	def_twiddle_macro			string
	def_fft_sum_function		string
	def_fft_diff_function		string
	def_two_sum_function		string
	def_two_diff_function		string
	
	batch_count					int
	sCount_init					int
	outer_loop_count			int
	inner_loop_count			int
	tps							int
	full_fft_length				int
	kernel_init_kw				int
	full_fft_length_set			bool
	kernel_fft_length			int
	kernel_fft_length_set		bool
	nameGenInit					bool
	mat_inverse_sCount			int
	final_mat_inverse_sCount	int
	twiddle_sCount				int
	inverse_fft					bool
	iofs_scaleFac				int
	global_grp_size				int
	local_grp_size				int
	fft_type					int
	total_fft_count				int
	isfloat						bool
	isterminal					bool
	fft_input_stride			int
}

func NewKernelWriter() *kernelWriter {
	tmp := new(kernelWriter)

	tmp.def_2s = preprocess_ifndef + " M_SCALE_1_2F\n\t#define M_SCALE_1_2F 0.500000000000000000000000f\n" + preprocess_endif + "\n\n"

	tmp.def_3s = tmp.def_2s + preprocess_ifndef + " M_SCALE_1_3F\n\t#define M_SCALE_1_3F 0.333333333333333333333333f\n" + preprocess_endif + "\n"
	tmp.def_3s += "\n" + preprocess_ifndef + " M_SQRT_3F\n\t#define M_SQRT_3F 1.73205080756887729352f\n" + preprocess_endif + "\n\n"

	tmp.def_4s = preprocess_ifndef + " M_SCALE_1_4F\n\t#define M_SCALE_1_4F 0.250000000000000000000000f\n" + preprocess_endif + "\n\n"

	tmp.def_5s = preprocess_ifndef + " M_SCALE_1_5F\n\t#define M_SCALE_1_5F 0.20000000000000000000f\n" + preprocess_endif + "\n"
	tmp.def_5s += "\n" + preprocess_ifndef + " M_COS_2PI_5F\n\t#define M_COS_2PI_5F 0.309016994374947451262869f\n" + preprocess_endif + "\n"
	tmp.def_5s += "\n" + preprocess_ifndef + " M_SIN_2PI_5F\n\t#define M_SIN_2PI_5F 0.951056516295153531181938f\n" + preprocess_endif + "\n"
	tmp.def_5s += "\n" + preprocess_ifndef + " M_COS_4PI_5F\n\t#define M_COS_4PI_5F -0.809016994374947451262869f\n" + preprocess_endif + "\n"
	tmp.def_5s += "\n" + preprocess_ifndef + " M_SIN_4PI_5F\n\t#define M_SIN_4PI_5F 0.587785252292473248125759f\n" + preprocess_endif + "\n\n"

	tmp.def_7s = preprocess_ifndef + " M_SCALE_1_7F\n\t#define M_SCALE_1_7F 0.142857142857142849212693f\n" + preprocess_endif + "\n"
	tmp.def_7s += "\n" + preprocess_ifndef + " M_COS_2PI_7F\n\t#define M_COS_2PI_7F 0.623489801858733594386308f\n" + preprocess_endif + "\n"
	tmp.def_7s += "\n" + preprocess_ifndef + " M_SIN_2PI_7F\n\t#define M_SIN_2PI_7F 0.781831482468029803634124f\n" + preprocess_endif + "\n"
	tmp.def_7s += "\n" + preprocess_ifndef + " M_COS_4PI_7F\n\t#define M_COS_4PI_7F -0.222520933956314337365257f\n" + preprocess_endif + "\n"
	tmp.def_7s += "\n" + preprocess_ifndef + " M_SIN_4PI_7F\n\t#define M_SIN_4PI_7F 0.974927912181823730364272f\n" + preprocess_endif + "\n"
	tmp.def_7s += "\n" + preprocess_ifndef + " M_COS_6PI_7F\n\t#define M_COS_6PI_7F -0.900968867902419034976447f\n" + preprocess_endif + "\n"
	tmp.def_7s += "\n" + preprocess_ifndef + " M_SIN_6PI_7F\n\t#define M_SIN_6PI_7F 0.433883739117558231423999f\n" + preprocess_endif + "\n\n"

	tmp.def_8s = preprocess_ifndef + " M_SCALE_1_8F\n\t#define M_SCALE_1_8F 0.125000000000000000000000f\n" + preprocess_endif + "\n\n"

	tmp.batch_count = 1
	tmp.kernel_init_kw = 1
	tmp.full_fft_length_set = false
	tmp.kernel_fft_length_set = false
	tmp.nameGenInit = false

	tmp.setTwiddleMacro()
	tmp.setTwoDiffFunction()
	tmp.setTwoSumFunction()
	tmp.setFFTDiffFunction()
	tmp.setFFTSumFunction()

	tmp.genFFT2Macros()
	tmp.genFFT3Macros()
	tmp.genFFT4Macros()
	tmp.genFFT5Macros()
	tmp.genFFT7Macros()
	tmp.genFFT8Macros()

	return tmp
}

func (s *kernelWriter) setTwiddleMacro() {
	s.def_twiddle_macro = "#define twiddle_factor(k, angle, in) {\\\n\t";
	if (s.isfloat) {
		s.def_twiddle_macro += "float2 tw, v;\\\n\ttw.x = cos((float)k*angle);\\\n\ttw.y = sin((float)k*angle);\\\n\t";
	} else {
		s.def_twiddle_macro += "double2 tw, v;\\\n\ttw.x = cos((double)k*angle);\\\n\ttw.y = sin((double)k*angle);\\\n\t";
	}
	s.def_twiddle_macro += "v.x = fma(tw.x, in.x, 0.0f);\\\n\tv.x = fma(-1.0f * tw.y, in.y, v.x);\\\n\tv.y = fma(tw.x, in.y, 0.0f);\\\n\tv.y = fma(tw.y, in.x, v.y);\\\n\tin = v;\\\n}\n\n";
}

func (s *kernelWriter) setTwoDiffFunction() {
	precision_string := string("")
	if (s.isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	s.def_two_diff_function = "inline " + precision_string + " two_diff(" + precision_string + " a, " + precision_string + " b) {\n\t";
	s.def_two_diff_function += precision_string + " s;\n\t";
	s.def_two_diff_function += "if (fabs(a) > fabs(b)) {\n\t\t";
	s.def_two_diff_function += "s = fma(1.0f, a, -1.0f * b);\n\t";
	s.def_two_diff_function += "} else {\n\t\t";
	s.def_two_diff_function += "s = fma(-1.0f, b, a);\n\t";
	s.def_two_diff_function += "}\n\t";
	s.def_two_diff_function += "return s;\n";
	s.def_two_diff_function += "}\n\n";
}

func (s *kernelWriter) setTwoSumFunction() {
	precision_string := string("");
	if (s.isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	s.def_two_sum_function = "inline " + precision_string + " two_sum(" + precision_string + " a, " + precision_string + " b) {\n\t";
	s.def_two_sum_function += precision_string + " s;\n";
	s.def_two_sum_function += "if (fabs(a) > fabs(b)) {\n\t\t";
	s.def_two_sum_function += "s = fma(1.0f, a, b);\n\t";
	s.def_two_sum_function += "} else {\n\t\t";
	s.def_two_sum_function += "s = fma(1.0f, b, a);\n\t";
	s.def_two_sum_function += "}\n\t";
	s.def_two_sum_function += "return s;\n";
	s.def_two_sum_function += "}\n\n";
}

func (s *kernelWriter) setFFTDiffFunction() {
	precision_string := string("")
	if (s.isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	s.def_fft_diff_function = "inline " + precision_string + "2 fft_diff(" + precision_string + "2 a, " + precision_string + "2 b) {\n\t";
	s.def_fft_diff_function += precision_string + "2 s;\n\t";
	s.def_fft_diff_function += "if (fabs(a.x) > fabs(b.x)) {\n\t\t";
	s.def_fft_diff_function += "s.x = fma(1.0f, a.x, -1.0f * b.x);\n\t";
	s.def_fft_diff_function += "} else {\n\t\t";
	s.def_fft_diff_function += "s.x = fma(-1.0f, b.x, a.x);\n\t";
	s.def_fft_diff_function += "}\n\t";
	s.def_fft_diff_function += "if (fabs(a.y) > fabs(b.y)) {\n\t\t";
	s.def_fft_diff_function += "s.y = fma(1.0f, a.y, -1.0f * b.y);\n\t";
	s.def_fft_diff_function += "} else {\n\t\t";
	s.def_fft_diff_function += "s.y = fma(-1.0f, b.y, a.y);\n\t";
	s.def_fft_diff_function += "}\n\t";
	s.def_fft_diff_function += "return s;\n";
	s.def_fft_diff_function += "}\n\n";
}

func (s *kernelWriter) setFFTSumFunction() {
	precision_string := string("")
	if (s.isfloat) {
		precision_string = "float";
	} else {
		precision_string = "double";
	}
	s.def_fft_sum_function = "inline " + precision_string + "2 fft_sum(" + precision_string + "2 a, " + precision_string + "2 b) {\n\t";
	s.def_fft_sum_function += precision_string + "2 s;\n\t";
	s.def_fft_sum_function += "if (fabs(a.x) > fabs(b.x)) {\n\t\t";
	s.def_fft_sum_function += "s.x = fma(1.0f, a.x, b.x);\n\t";
	s.def_fft_sum_function += "} else {\n\t\t";
	s.def_fft_sum_function += "s.x = fma(1.0f, b.x, a.x);\n\t";
	s.def_fft_sum_function += "}\n\t";
	s.def_fft_sum_function += "if (fabs(a.y) > fabs(b.y)) {\n\t\t";
	s.def_fft_sum_function += "s.y = fma(1.0f, a.y, b.y);\n\t";
	s.def_fft_sum_function += "} else {\n\t\t";
	s.def_fft_sum_function += "s.y = fma(1.0f, b.y, a.y);\n\t";
	s.def_fft_sum_function += "}\n\t";
	s.def_fft_sum_function += "return s;\n";
	s.def_fft_sum_function += "}\n\n";
}

func (s *kernelWriter) genFFT2Macros() {
}

func (s *kernelWriter) genFFT3Macros() {
}

func (s *kernelWriter) genFFT4Macros() {
}

func (s *kernelWriter) genFFT5Macros() {
}

func (s *kernelWriter) genFFT7Macros() {
}

func (s *kernelWriter) genFFT8Macros() {
}

func (s *kernelWriter) SetPrecision(x bool) {
	if (s.isfloat != x) {
		s.isfloat = x;
		s.setTwiddleMacro();
		s.setFFTDiffFunction();
		s.setFFTSumFunction();
		s.setTwoDiffFunction();
		s.setTwoSumFunction();
		s.genFFT2Macros();
		s.genFFT3Macros();
		s.genFFT4Macros();
		s.genFFT5Macros();
		s.genFFT7Macros();
		s.genFFT8Macros();
	}
}

func (s *kernelWriter) SetIOFSScaleFactor(x int) {
	s.iofs_scaleFac = x
}

func (s *kernelWriter) SetGlobalGroupSize(x int) {
	s.global_grp_size = x
}

func (s *kernelWriter) SetLocalGroupSize(x int) {
	s.local_grp_size = x
}

func (s *kernelWriter) SetInitSCount(x int) {
	s.sCount_init = x
}

func (s *kernelWriter) SetTwiddleSCount(x int) {
	s.twiddle_sCount = x
}

func (s *kernelWriter) SetOuterCount(x int) {
	s.outer_loop_count = x
}

func (s *kernelWriter) SetInnerCount(x int) {
	s.inner_loop_count = x
}

func (s *kernelWriter) SetBatchCount(x int) {
	s.batch_count = x
}

func (s *kernelWriter) SetFFTInputStride(x int) {
	s.fft_input_stride = x
}

func (s *kernelWriter) SetMatrixInverseSCount(x int) {
	s.mat_inverse_sCount = x
}

func (s *kernelWriter) SetFinalMatrixInverseSCount(x int) {
	s.final_mat_inverse_sCount = x
}

func (s *kernelWriter) SetInverseDirection(x bool) {
	s.inverse_fft = x
}

func (s *kernelWriter) SetFFTType(x int) {
	s.fft_type = x
}

func (s *kernelWriter) SetTotalFFTCount(x int) {
	s.total_fft_count = x
}

func (s *kernelWriter) SetInitKW(x int) {
	s.kernel_init_kw = x
}

func (s *kernelWriter) SetTPS(x int) {
	s.tps = x;
	s.SetFFTInputStride(s.tps);
}

func (s *kernelWriter) SetFFTLength (x int) {
	s.full_fft_length = x;
	s.full_fft_length_set = true;
	if (s.full_fft_length_set && s.kernel_fft_length_set) {
		s.SetTPS(s.full_fft_length / s.kernel_fft_length);
	}
}

func (s *kernelWriter) SetKernelFFTLength (x int) {
	s.kernel_fft_length = x;
	s.kernel_fft_length_set = true;
	if (s.full_fft_length_set && s.kernel_fft_length_set) {
		s.SetTPS(s.full_fft_length / s.kernel_fft_length);
	}
}
