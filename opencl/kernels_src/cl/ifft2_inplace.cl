// C2C FFT kernel loop
// EnqueueNDRangeKernel should be given 1-D global (FFT2_GLOBAL_SZ) and local (FFT2_LOCAL_SZ) sizes
// Need to define the following during compile time:
// FFT2_LOCAL_SZ
// FFT2_GLOBAL_SZ
// NFFT2: total number of FFT2 to perform
// NsFFT2: total number of stages of FFT2 to perform
// K_W: actual length of FFT that needs to be performed
// LFFT2: length of FFT kernel performs (= 2 for FFT2)
// TPS: number of work-items working on a single batch
__kernel void
ifft2_c2c_long_interleaved_inp(__global float2* dataIn, __global float2* buffer) {
	int lid = get_local_id(0);
	int itr = get_group_id(0)*FFT2_LOCAL_SZ + lid;
	int sCount = 0;
	float angle_ = 2.0f * M_PI_F / (float)(K_W);

	while (sCount < Ns) {
		int idx = itr;
		int lti = idx % TPS;
		int iofs = idx - lti;
		iofs /= TPS;
		iofs *= K_W;
		while (idx < NFFT2) {
			// Load input into registers
			float2 in0, in1, V;
			int Idout0 = idx + iofs;
			int Idout1 = idx + iofs + K_W / LFFT2;

			in0 = dataIn[Idout0];
			in1 = dataIn[Idout1];

			// Scale values
			in0 *= 0.5000000000f;
			in1 *= 0.5000000000f;

			// Perform length-2 inverse FFT calculations
			if (sCount > 0) {
				float angle = angle_ * (float)(lti);
				twiddle_factor(1, angle, in1);
			}

			FFT2(in0, in1);

			// Store results in place if we do not need to transpose
			// Otherwise, we have to store in sequence in the buffer, and copy back later
			if (sCount == NsFFFT2 - 1) {
				dataIn[Idout0] = in0;
				dataIn[Idout1] = in1;
			} else {
				__local float2 block_cache[FFT2_LOCAL_SZ*2];
				block_cache[2*lid] = in0;
				block_cache[2*lid+1] = in1;
				barrier(CLK_LOCAL_MEM_FENCE);
				in0 = block_cache[lid];
				in1 = block_cache[lid + FFT2_LOCAL_SZ];
				buffer[Idout0] = in0;
				buffer[Idout0 + FFT2_LOCAL_SZ] = in1;
				barrier(CLK_GLOBAL_MEM_FENCE);
			}

			// Increment idx to go to next logical work item
			idx += FFT2_GLOBAL_SZ;
			angle_ *= 2.0f;
		}
		// If we needed to do tranpose, the FFT result is in buffer and not in dataIn
		// because we do not want to overwrite the input signals. So we need to copy
		// from buffer back to dataIn before the next FFT
		if (sCount < NsFFT2 - 1) {
			idx = itr;
			while (idx < 2*NFFT2) {
				float2 buf_val = buffer[idx];
				dataIn[idx] = buffer[idx];
				idx += FFT2_GLOBAL_SZ;
			}
		}
		sCount++;
	}
}
