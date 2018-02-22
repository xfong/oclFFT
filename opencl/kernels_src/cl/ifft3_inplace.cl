// C2C FFT kernel loop
// EnqueueNDRangeKernel should be given 1-D global (FFT3_GLOBAL_SZ) and local (FFT3_LOCAL_SZ) sizes
// Need to define the following during compile time:
// FFT3_LOCAL_SZ
// FFT3_GLOBAL_SZ
// NFFT3: total number of FFT3 to perform
// NsFFT3: total number of stages of FFT3 to perform
// K_W: actual length of FFT that needs to be performed
// LFFT3: length of FFT kernel performs (= 3 for FFT3)
// TPS: number of work-items working on a single batch
__kernel void
ifft3_c2c_long_interleaved_inp(__global float2* dataIn, __global float2* buffer) {
	int lid = get_local_id(0);
	int itr = get_group_id(0)*FFT3_LOCAL_SZ + lid;
	int sCount = 0;
	float angle_ = 2.0f * M_PI_F / (float)(K_W);

	while (sCount < NsFFFT3) {
		int idx = itr;
		int lti = idx % TPS;
		int iofs = idx - lti;
		iofs /= TPS;
		iofs *= K_W;
		while (idx < NFFT3) {
			// Load input into registers
			float2 in0, in1, in2;
			int Idout0 = idx + iofs;
			int Idout1 = Idout0+ K_W / LFFT3;
			int Idout2 = Idout0+2*K_W / LFFT3;

			in0 = dataIn[Idout0];
			in1 = dataIn[Idout1];
			in2 = dataIn[Idout2];

			// Scale values
			in0 *= 0.333333333333f;
			in1 *= 0.333333333333f;
			in2 *= 0.333333333333f;

			// Perform length-3 inverse FFT calculations
			if (sCount > 0) {
				float angle = angle_ * (float)(lti);
				twiddle_factor(1, angle, in1);
				twiddle_factor(2, angle, in2);
			}

			IFFT3(in0, in1, in2);

			// Store results in place if we do not need to transpose
			// Otherwise, we have to store in sequence in the buffer, and copy back later
			if (sCount == NsFFT3 - 1) {
				dataIn[Idout0] = in0;
				dataIn[Idout1] = in1;
				dataIn[Idout2] = in2;
			} else {
				__local float2 block_cache[FFT3_LOCAL_SZ*3];
				block_cache[3*lid] = in0;
				block_cache[3*lid+1] = in1;
				block_cache[3*lid+2] = in2;
				barrier(CLK_LOCAL_MEM_FENCE);
				in0 = block_cache[lid];
				in1 = block_cache[lid + FFT3_LOCAL_SZ];
				in2 = block_cache[lid + 2*FFT3_LOCAL_SZ];
				buffer[Idout0] = in0;
				buffer[Idout0 + FFT3_LOCAL_SZ] = in1;
				buffer[Idout0 + 2*FFT3_LOCAL_SZ] = in2;
				barrier(CLK_GLOBAL_MEM_FENCE);
			}

			// Increment idx to go to next logical work item
			idx += FFT3_GLOBAL_SZ;
			angle_ *= 3.0f;
		}
		// If we needed to do tranpose, the FFT result is in buffer and not in dataIn
		// because we do not want to overwrite the input signals. So we need to copy
		// from buffer back to dataIn before the next FFT
		if (sCount < NsFFT3 - 1) {
			idx = itr;
			while (idx < 3*NFFT3) {
				float2 buf_val = buffer[idx];
				dataIn[idx] = buffer[idx];
				idx += FFT3_GLOBAL_SZ;
			}
		}
		sCount++;
	}
}
