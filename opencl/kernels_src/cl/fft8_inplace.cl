// C2C FFT kernel loop
// EnqueueNDRangeKernel should be given 1-D global (FFT8_GLOBAL_SZ) and local (FFT8_LOCAL_SZ) sizes
// Need to define the following during compile time:
// FFT8_LOCAL_SZ
// FFT8_GLOBAL_SZ
// NFFT8: total number of FFT8 to perform
// NsFFT8: total number of stages of FFT8 to perform
// K_W: actual length of FFT that needs to be performed
// LFFT8: length of FFT kernel performs (= 8 for FFT8)
// TPS: number of work-items working on a single batch (K_W / LFFT8)
__kernel void
fft8_c2c_long_interleaved_oop(__global float2* dataIn, __global float2* buffer) {
	int lid = get_local_id(0);
	int itr = get_group_id(0)*FFT8_LOCAL_SZ + lid;
	int sCount = 0;
	float angle_ = -2.0f * M_PI_F / (float)(K_W);

	while (sCount < Ns) {
		int idx = itr;
		int lti = idx % TPS;
		int iofs = idx - lti;
		iofs /= TPS;
		iofs *= K_W;
		while (idx < N) {
			// Load input into registers
			float2 in0, in1, in2, in3, in4, in5, in6, in7;
			int Idout0 = idx + iofs;
			int Idout1 = Idout0+K_W / LFFT8;
			int Idout2 = Idout0+2*K_W / LFFT8;
			int Idout3 = Idout0+3*K_W / LFFT8;
			int Idout4 = Idout0+4*K_W / LFFT8;
			int Idout5 = Idout0+5*K_W / LFFT8;
			int Idout6 = Idout0+6*K_W / LFFT8;
			int Idout7 = Idout0+7*K_W / LFFT8;

			in0 = dataIn[Idout0];
			in1 = dataIn[Idout1];
			in2 = dataIn[Idout2];
			in3 = dataIn[Idout3];
			in4 = dataIn[Idout4];
			in5 = dataIn[Idout5];
			in6 = dataIn[Idout6];
			in7 = dataIn[Idout7];

			// Perform length-8 forward FFT calculations
			if (sCount > 0) {
				float angle = angle_ * (float)(lti);
				twiddle_factor(1, angle, in1);
				twiddle_factor(2, angle, in2);
				twiddle_factor(3, angle, in3);
				twiddle_factor(4, angle, in4);
				twiddle_factor(5, angle, in5);
				twiddle_factor(6, angle, in6);
				twiddle_factor(7, angle, in7);
			}

			IFFT8(in0, in1, in2, in3, in4, in5, in6, in7);

			// Store results in place if we do not need to transpose
			// Otherwise, we have to store in sequence in the buffer, and copy back later
			if (sCount == NsFFT8 - 1) {
				buffer[Idout0] = in0;
				buffer[Idout1] = in1;
				buffer[Idout2] = in2;
				buffer[Idout3] = in3;
				buffer[Idout4] = in4;
				buffer[Idout5] = in5;
				buffer[Idout6] = in6;
				buffer[Idout7] = in7;
			} else {
				__local float2 block_cache[FFT8_LOCAL_SZ*8];
				block_cache[8*lid] = in0;
				block_cache[8*lid+1] = in1;
				block_cache[8*lid+2] = in2;
				block_cache[8*lid+3] = in3;
				block_cache[8*lid+4] = in4;
				block_cache[8*lid+5] = in5;
				block_cache[8*lid+6] = in6;
				block_cache[8*lid+7] = in7;
				barrier(CLK_LOCAL_MEM_FENCE);
				in0 = block_cache[lid];
				in1 = block_cache[lid + FFT8_LOCAL_SZ];
				in2 = block_cache[lid + 2*FFT8_LOCAL_SZ];
				in3 = block_cache[lid + 3*FFT8_LOCAL_SZ];
				in4 = block_cache[lid + 4*FFT8_LOCAL_SZ];
				in5 = block_cache[lid + 5*FFT8_LOCAL_SZ];
				in6 = block_cache[lid + 6*FFT8_LOCAL_SZ];
				in7 = block_cache[lid + 7*FFT8_LOCAL_SZ];
				buffer[Idout0] = in0;
				buffer[Idout0 + FFT8_LOCAL_SZ] = in1;
				buffer[Idout0 + 2*FFT8_LOCAL_SZ] = in2;
				buffer[Idout0 + 3*FFT8_LOCAL_SZ] = in3;
				buffer[Idout0 + 4*FFT8_LOCAL_SZ] = in4;
				buffer[Idout0 + 5*FFT8_LOCAL_SZ] = in5;
				buffer[Idout0 + 6*FFT8_LOCAL_SZ] = in6;
				buffer[Idout0 + 7*FFT8_LOCAL_SZ] = in7;
				barrier(CLK_GLOBAL_MEM_FENCE);
			}

			// Increment idx to go to next logical work-item
			idx += FFT8_GLOBAL_SZ;
			angle_ *= 8.0f;
		}
		// If we needed to do tranpose, the FFT result is in buffer and not in dataIn
		// because we do not want to overwrite the input signals. So we need to copy
		// from buffer back to dataIn before the next FFT
		if (sCount < NsFFT8 - 1) {
			idx = itr;
			while (idx < 8*NFFT8) {
				float2 buf_val = buffer[idx];
				dataIn[idx] = buffer[idx];
				idx += FFT8_GLOBAL_SZ;
			}
		}
		sCount++;
	}
}
