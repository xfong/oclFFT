// C2C FFT kernel loop
// EnqueueNDRangeKernel should be given 1-D global (FFT5_GLOBAL_SZ) and local (FFT5_LOCAL_SZ) sizes
// Need to define the following during compile time:
// FFT5_LOCAL_SZ
// FFT5_GLOBAL_SZ
// NFFT5: total number of FFT5 to perform
// NsFFT5: total number of stages of FFT5 to perform
// K_W: actual length of FFT that needs to be performed
// LFFT5: length of FFT kernel performs (= 5 for FFT5)
// TPS: number of work-items working on a single batch (K_W / LFFT5)
__kernel void
ifft5_c2c_long_interleaved_oop(__global float2* dataIn, __global float2* buffer) {
	int lid = get_local_id(0);
	int itr = get_group_id(0)*FFT5_LOCAL_SZ + lid;
	int sCount = 0;
	float angle_ = 2.0f * M_PI_F / (float)(K_W);

	while (sCount < Ns) {
		int idx = itr;
		int lti = idx % TPS;
		int iofs = idx - lti;
		iofs /= TPS;
		iofs *= K_W;
		while (idx < N) {
			// Load input into registers
			float2 in0, in1, in2, in3, in4;
			int Idout0 = idx + iofs;
			int Idout1 = Idout0+K_W / LFFT5;
			int Idout2 = Idout0+2*K_W / LFFT5;
			int Idout3 = Idout0+3*K_W / LFFT5;
			int Idout4 = Idout0+4*K_W / LFFT5;

			in0 = dataIn[Idout0];
			in1 = dataIn[Idout1];
			in2 = dataIn[Idout2];
			in3 = dataIn[Idout3];
			in4 = dataIn[Idout4];

			// Scale values
			in0 *= 0.200000000000f;
			in1 *= 0.200000000000f;
			in2 *= 0.200000000000f;
			in3 *= 0.200000000000f;
			in4 *= 0.200000000000f;

			// Perform length-5 inverse FFT calculations
			if (sCount > 0) {
				float angle = angle_ * (float)(lti);
				twiddle_factor(1, angle, in1);
				twiddle_factor(2, angle, in2);
			}

			IFFT5(in0, in1, in2, in3, in4);

			// Store results in place if we do not need to transpose
			// Otherwise, we have to store in sequence in the buffer, and copy back later
			if (sCount == NsFFT5 - 1) {
				dataIn[Idout0] = in0;
				dataIn[Idout1] = in1;
				dataIn[Idout2] = in2;
				dataIn[Idout3] = in3;
				dataIn[Idout4] = in4;
			} else {
				__local float2 block_cache[FFT5_LOCAL_SZ*5];
				block_cache[5*lid] = in0;
				block_cache[5*lid+1] = in1;
				block_cache[5*lid+2] = in2;
				block_cache[5*lid+3] = in3;
				block_cache[5*lid+4] = in4;
				barrier(CLK_LOCAL_MEM_FENCE);
				in0 = block_cache[lid];
				in1 = block_cache[lid + FFT5_LOCAL_SZ];
				in2 = block_cache[lid + 2*FFT5_LOCAL_SZ];
				in3 = block_cache[lid + 3*FFT5_LOCAL_SZ];
				in4 = block_cache[lid + 4*FFT5_LOCAL_SZ];
				buffer[Idout0] = in0;
				buffer[Idout0 + FFT5_LOCAL_SZ] = in1;
				buffer[Idout0 + 2*FFT5_LOCAL_SZ] = in2;
				buffer[Idout0 + 3*FFT5_LOCAL_SZ] = in3;
				buffer[Idout0 + 4*FFT5_LOCAL_SZ] = in4;
				barrier(CLK_GLOBAL_MEM_FENCE);
			}

			// Increment idx to go to next logical work item
			idx += FFT5_GLOBAL_SZ;
			angle_ *= 5.0f;
		}
		// If we needed to do tranpose, the FFT result is in buffer and not in dataIn
		// because we do not want to overwrite the input signals. So we need to copy
		// from buffer back to dataIn before the next FFT
		if (sCount < NsFFT5 - 1) {
			idx = itr;
			while (idx < 5*NFFT5) {
				float2 buf_val = buffer[idx];
				dataIn[idx] = buffer[idx];
				idx += FFT5_GLOBAL_SZ;
			}
		}
		sCount++;
	}
}
