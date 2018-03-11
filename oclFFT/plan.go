package oclFFT

import (
	"github.com/source/oclFFT/opencl/cl"
)

type OclFFTPlanType int

const (
	OCLFFTPLAN_1D_TRANSFORM	OclFFTPlanType		= 1
	OCLFFTPLAN_2D_TRANSFORM	OclFFTPlanType		= 2
	OCLFFTPLAN_3D_TRANSFORM	OclFFTPlanType		= 3
)

type OclFFTPlanPrecision int

const (
	OCLFFTPLAN_PRECISION_HALF	OclFFTPlanPrecision		= 1
	OCLFFTPLAN_PRECISION_FLOAT	OclFFTPlanPrecision		= 2
	OCLFFTPLAN_PRECISION_DOUBLE	OclFFTPlanPrecision		= 3
)

type OclFFTPlanDirection int

const (
	OCLFFTPLAN_FORWARD_TRANSFORM	OclFFTPlanDirection		=  1
	OCLFFTPLAN_BACKWARD_TRANSFORM	OclFFTPlanDirection		= -1
)

type plan struct {
	plType			OclFFTPlanType
	plPrecision		OclFFTPlanPrecision
	plDirection		OclFFTPlanDirection
	plans			[]oclFFTPlanBase
	num_plans		int
	dims			[3]int
	baked			bool
	queue			cl.CommandQueue
	platform		cl.Platform
	context			cl.Context
	device			cl.Device
}

func NewPlan() plan {
	s := plan{plType: OCLFFTPLAN_1D_TRANSFORM, plPrecision: OCLFFTPLAN_PRECISION_FLOAT, plDirection: OCLFFTPLAN_FORWARD_TRANSFORM, num_plans: 0, dims: [3]int{1, 1, 1}, baked: false}
	return s
}

func (s *plan) GetType() OclFFTPlanType {
	return s.plType
}

func (s *plan) GetPrecision() OclFFTPlanPrecision {
	return s.plPrecision
}

func (s *plan) GetDirection() OclFFTPlanDirection {
	return s.plDirection
}

func (s *plan) GetDimension() [3]int {
	return s.dims
}

func (s *plan) CheckExecution() bool {
	for _, vv := range s.plans {
		if vv.CheckBlocked() {
		} else {
			return true
		}
	}
	return false
}

func (s *plan) setType(inType OclFFTPlanType) {
	if (s.plType == inType) {
	} else {
		s.plType = inType
		s.baked = false
	}
}

func (s *plan) SetPrecision(inPrecision OclFFTPlanPrecision) {
	if (s.plPrecision == inPrecision) {
	} else {
		s.plPrecision = inPrecision
		s.baked = false
	}
}

func (s *plan) SetDirection(inDirection OclFFTPlanDirection) {
	if (s.plDirection == inDirection) {
	} else {
		s.plDirection = inDirection
		s.baked = false
	}
}

func (s *plan) SetDimension(inDimension [3]int) {
	testCond := true
	for idx, vv := range s.dims {
		if vv != inDimension[idx] {
			testCond = false
		}
	}

	if (testCond) {
	} else {
		s.dims = inDimension
		s.baked = false
		counts := 0
		for _, vv := range s.dims {
			if vv > 1 {
				counts++
			}
		}
		if s.num_plans > 0 {
			for _, vv := range s.plans {
				vv.destroy()
			}
			s.num_plans = 0
		}
		s.plans = make([]oclFFTPlanBase, counts)
		for ii := 0; ii < counts; ii++ {
			pl := NewOclFFTPlanBase()
			s.plans[ii] = pl
			pl.SetTransformSize(inDimension[ii])
		}
		switch counts {
		case 3:
			// 3D
			s.setType(OCLFFTPLAN_3D_TRANSFORM)
		case 2:
			// 2D
			s.setType(OCLFFTPLAN_2D_TRANSFORM)
		default:
			// 1D
			s.setType(OCLFFTPLAN_1D_TRANSFORM)
		}
	}
}

func (s *plan) destroy() {
	for _, vv := range s.plans {
		vv.destroy()
	}
}
