package oclFFT

import (
	"strings"

	"github.com/source/oclFFT/opencl/cl"
)

type oclFFTPlanBase struct {
	buffers					[]cl.MemObject
	kernels					[]cl.Kernel
	kernel_codes			[]string
	programs				[]cl.Program
	Transform_Length		int	=	1
	EntryKW					[...]int = {1, 1, 1, 1, 1, 1, 1}
	StageCnts				[...]int = {-1, -1, -1, -1, -1, -1}
}

func (s *oclFFTPlanBase) SetTransformSize(x int) {
	
}