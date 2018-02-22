// C2C length 4 FFT kernel loop
// EnqueueNDRangeKernel should be given 1-D global (FFT4_GLOBAL_SZ) and local (FFT4_LOCAL_SZ) sizes
// Need to define the following during compile time:
// FFT4_LOCAL_SZ
// FFT4_GLOBAL_SZ
// NFFT4: total number of FFT4 to perform
// NsFFT4: total number of stages of FFT4 to perform
// K_W: actual length of FFT that needs to be performed
// LFFT4: length of FFT kernel performs (= 4 for FFT4)
// TPS: number of work-items working on a single batch (K_W / LFFT4)
__kernel void
ifft4_c2c_long_interleaved_inp(__global float2* dataIn, __global float2* buffer) {
	int lid = get_local_id(0);
	int itr = get_group_id(0)*FFT4_LOCAL_SZ + lid;
	int sCount = 0;
	float angle_ = 2.0f * M_PI_F / (float)(K_W);

	while (sCount < NsFFT4) {
		int idx = itr;
		int lti = idx % TPS;
		int iofs = idx - lti;
		iofs /= TPS;
		iofs *= K_W;
		while (idx < NFFT4) {
			// Load input into registers
			float2 in0, in1, in2, in3;
			int Idout0 = idx + iofs;
			int Idout1 = Idout0+K_W / LFFT4;
			int Idout2 = Idout0+2*K_W / LFFT4;
			int Idout3 = Idout0+3*K_W / LFFT4;

			in0 = dataIn[Idout0];
			in1 = dataIn[Idout1];
			in2 = dataIn[Idout2];
			in3 = dataIn[Idout3];

			// Scale values
			in0 *= 0.2500000f;
			in1 *= 0.2500000f;
			in2 *= 0.2500000f;
			in3 *= 0.2500000f;

			// Perform length-4 inverse FFT calculations
			if (sCount > 0) {
				float angle = angle_ * (float)(lti);
				twiddle_factor(1, angle, in1);
				twiddle_factor(2, angle, in2);
				twiddle_factor(3, angle, in3);
			}

			IFFT4(in0, in1, in2, in3);

			// Store results in place if we do not need to transpose
			// Otherwise, we have to store in sequence in the buffer, and copy back later
			if (sCount == NsFFT4 - 1) {
				dataIn[Idout0] = in0;
				dataIn[Idout1] = in1;
				dataIn[Idout2] = in2;
				dataIn[Idout3] = in3;
			} else {
				__local float2 block_cache[FFT4_LOCAL_SZ*4];
				block_cache[4*lid] = in0;
				block_cache[4*lid+1] = in1;
				block_cache[4*lid+2] = in2;
				block_cache[4*lid+3] = in3;
				barrier(CLK_LOCAL_MEM_FENCE);
				in0 = block_cache[lid];
				in1 = block_cache[lid + FFT4_LOCAL_SZ];
				in2 = block_cache[lid + 2*FFT4_LOCAL_SZ];
				in3 = block_cache[lid + 3*FFT4_LOCAL_SZ];
				buffer[Idout0] = in0;
				buffer[Idout0 + FFT4_LOCAL_SZ] = in1;
				buffer[Idout0 + 2*FFT4_LOCAL_SZ] = in2;
				buffer[Idout0 + 3*FFT4_LOCAL_SZ] = in3;
				barrier(CLK_GLOBAL_MEM_FENCE);
			}

			// Increment idx to go to next logical work item
			idx += FFT4_GLOBAL_SZ;
			angle_ *= 4.0f;
		}
		// If we needed to do tranpose, the FFT result is in buffer and not in dataIn
		// because we do not want to overwrite the input signals. So we need to copy
		// from buffer back to dataIn before the next FFT
		if (sCount < NsFFT4 - 1) {
			idx = itr;
			while (idx < 4*NFFT4) {
				float2 buf_val = buffer[idx];
				dataIn[idx] = buffer[idx];
				idx += FFT4_GLOBAL_SZ;
			}
		}
		sCount++;
	}
}
