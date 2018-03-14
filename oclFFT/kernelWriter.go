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
	
}

func NewKernelWriter() kernelWriter {
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
	
	return *tmp
}