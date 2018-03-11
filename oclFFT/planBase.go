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
	Transform_Length		= (int)(1)
	EntryKW					= [...]int { 1,  1,  1,  1,  1,  1, 1}
	StageCnts				= [...]int {-1, -1, -1, -1, -1, -1}
	Blocked					= false
}

func (s *oclFFTPlanBase) SetTransformSize(x int) error {
	if (Blocked == false) {
		s.Transform_Length = x
		s.StageCnts, s.EntryKW = countStages(x)
		return nil
	} else {
		return ErrUnknown
	}
}

func (s *oclFFTPlanBase) GetTransformSize() int {
	return s.Transform_Length
}

func (s *oclFFTPlanBase) GetKernelCode() []string {
	return s.kernel_codes
}

func (s *oclFFTPlanBase) GetEntryKW() [7]int {
	return s.EntryKW
}

func (s *oclFFTPlanBase) GetStageCnts() [6]int {
	return s.StageCnts
}
