package oclFFT

import (
	"github.com/source/oclFFT/opencl/cl"
)

var	FFTKernels	= [...]int {8, 7, 5, 4, 3, 2}

type oclFFTPlanBase struct {
	buffers					[]cl.MemObject
	kernels					[]cl.Kernel
	kernel_codes			[]string
	programs				[]cl.Program
	Transform_Length		int
	EntryKW					[7]int
	StageCounts				[6]int
	Blocked					bool
}

func NewOclFFTPlanBase() oclFFTPlanBase {
	tmpStruct := oclFFTPlanBase{EntryKW: [7]int{-1, -1, -1, -1, -1, -1, -1}, StageCounts: [6]int{0, 0, 0, 0, 0, 0}, Blocked: false}
	return tmpStruct
}

func (s *oclFFTPlanBase) SetTransformSize(x int) error {
	if (s.Blocked == false) {
		s.Transform_Length = x
		s.StageCounts, s.EntryKW = countStages(x)
		return nil
	} else {
		return cl.ErrUnknown
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

func (s *oclFFTPlanBase) GetStageCounts() [6]int {
	return s.StageCounts
}

func (s *oclFFTPlanBase) GetNumberOfBuffers() int {
	return len(s.buffers)
}

func (s *oclFFTPlanBase) GetNumberOfKernels() int {
	return len(s.kernels)
}

func (s *oclFFTPlanBase) GetNumberOfPrograms() int {
	return len(s.programs)
}

func (s *oclFFTPlanBase) CheckBlocked() bool {
	return s.Blocked
}

func (s *oclFFTPlanBase) destroy() {
	for _, vv := range s.buffers {
		vv.Release()
	}
	for _, vv := range s.kernels {
		vv.Release()
	}
	for _, vv := range s.programs {
		vv.Release()
	}
}
