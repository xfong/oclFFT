package oclFFT

import (
	"github.com/source/oclFFT/opencl/cl"
)

type planType		uint

const (
	OCLFFT_1D_TRANSFORM		planType = 0
	OCLFFT_2D_TRANSFORM		planType = 1
	OCLFFT_3D_TRANSFORM		planType = 2
)

type planPrecision	uint

const (
	OCLFFT_PRECISION_HALF		planPrecision = 0
	OCLFFT_PRECISION_FLOAT		planPrecision = 1
	OCLFFT_PRECISION_DOUBLE		planPrecision = 2
)

type planDirection	int

const (
	OCLFFT_FORWARD_PLAN		planDirection =		 1
	OCLFFT_BACKWARD_PLAN		planDirection =		-1
)

type planAbs struct {
	type		planType
	precision	planPrecision
	direction	planDirection
	plans		[]*oclFFTPlanBase
	num_plans	int	=	-1
	queue		cl.CommandQueue
}

func (s *planAbs) SetTransformPrecision(inPrec planPrecision) {
	s.precision = inPrec
}

func (s *planAbs) SetTransformDirection(inDirection planDirection) {
	s.direction = inDirection
}

func (s *planAbs) SetTransformType(inType planType) {
	if s.num_plans > 1 {
		for (; s.num_plans > 0; s.num_plans--) {
			currIdx := s.num_plans - 1
			tmpPlanPtr := *s.plans[currIdx]
			ret := toError(tmpPlanPtr.Destroy())
			if ret == nil {
				s.plans[currIdx] = nil
			}
		}
	}
	
	s.type = inType
	num_plans = ((uint)(inType)) + 1
	s.plans = make([]*oclFFTPlanBase, num_plans)
}

