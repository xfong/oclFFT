#include <sstream>
#include "plan.h"

#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x100000)
 
int main(int argc, char *argv[])
{
	cl_device_id *device_id_arr_ptr;
	cl_device_id device_id = NULL;
	cl_context context = NULL;
	cl_command_queue command_queue = NULL;
	cl_mem memobj = NULL;
	// cl_program program = NULL;
	// cl_kernel kernel = NULL;
	cl_platform_id *platform_id_arr_ptr;
	cl_platform_id platform_id = NULL;
	cl_uint ret_num_devices;
	cl_uint ret_num_platforms;
	cl_int ret;
	int32_t fft_length = -1;
	uint32_t platform_num = 0;
	uint32_t gpu_num = 0;
 
	if (argc == 6) {
		for (int idx = 1; idx < argc;) {
			if (strcmp("-platform", argv[idx])) {
				std::istringstream ss(argv[idx+1]);
				if (!(ss >> platform_num)) {
					std::cerr << "Invalid platform number: " << argv[idx+1] << '\n';
					return -1;
				}
				idx += 2;
				ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
				platform_id_arr_ptr = new cl_platform_id[ ret_num_platforms ];
				ret = clGetPlatformIDs(ret_num_platforms, platform_id_arr_ptr, NULL);
				platform_id = platform_id_arr_ptr[platform_num];
			} else if (strcmp("-gpu", argv[idx])) {
				std::istringstream ss(argv[idx+1]);
				if (!(ss >> gpu_num)) {
					std::cerr << "Invalid gpu number: " << argv[idx+1] << '\n';
					return -1;
				}
				idx += 2;
			} else {
				std::istringstream ss(argv[idx]);
				if (!(ss >> fft_length)) {
					std::cerr << "Invalid fft_length: " << argv[idx] << '\n';
					return -1;
				}
				idx += 1;
			}
		}
	} else {
		std::wcout << "Fail to parse arguments!" << std::endl;
		return -1;
	}
	ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
	device_id_arr_ptr = new cl_device_id[ ret_num_devices ];
	ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, ret_num_devices, device_id_arr_ptr, NULL);
	device_id = device_id_arr_ptr[gpu_num];
	 
	// /* Create OpenCL context */
	context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
	 
	// /* Create Command Queue */
	command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
	 
	// /* Create Memory Buffer */
	memobj = clCreateBuffer(context, CL_MEM_READ_WRITE,MEM_SIZE * sizeof(char), NULL, &ret);
	 
	// /* Create Kernel Program from the source */
	// program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
	// (const size_t *)&source_size, &ret);
	 
	// /* Build Kernel Program */
	// ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	 
	// /* Create OpenCL Kernel */
	// kernel = clCreateKernel(program, "hello", &ret);
	 
	// /* Set OpenCL Kernel Parameters */
	// ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobj);
	 
	// /* Execute OpenCL Kernel */
	// ret = clEnqueueTask(command_queue, kernel, 0, NULL,NULL);
	 
	// /* Copy results from the memory buffer */
	// ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0,
	// MEM_SIZE * sizeof(char),string, 0, NULL, NULL);
	 
	// /* Display Result */
	// puts(string);
	 
	// /* Finalization */
	ret = clFlush(command_queue);
	ret = clFinish(command_queue);
	// ret = clReleaseKernel(kernel);
	// ret = clReleaseProgram(program);
	ret = clReleaseMemObject(memobj);
	ret = clReleaseCommandQueue(command_queue);
	ret = clReleaseContext(context);
	 
	// free(source_str);
	 
	return 0;
}
