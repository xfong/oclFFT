package main
//
//import (
//	"math"
//	"strings"
//)

var	FFTKernels	= [...]int {8, 7, 5, 4, 3, 2}

func countStages(x int) ([6]int, [7]int) {
	tmpCnt := [6]int {-1, -1, -1, -1, -1, -1}
	initKW := [7]int {1, 1, 1, 1, 1, 1, 1}

	for idx := range tmpCnt {
		tmpCnt[idx] = -1
		initKW[idx] = 1
	}

	if (x > 1) {
		fft_sz := x
		initKW[0] = fft_sz
		for idx, v := range FFTKernels {
			tmpCnt[idx], fft_sz = countKernStage(fft_sz, v)
			initKW[idx+1] = fft_sz
		}
	}

	return tmpCnt, initKW
}

func countKernStage(x, y int) (int, int) {
	// y gives kernel length
	// x gives the full FFT length
	if (y == 1) {
		return x, x
	} else {
		if (y < 1) {
			return -1, x
		}
	}
	StageCnt	:=	(int)(0)
	OutLength	:=	x
	flag		:=	true
	for flag {
		tmpLength := OutLength / y
		if (tmpLength * y == OutLength) {
			StageCnt++
			OutLength = tmpLength
		} else {
			return StageCnt, OutLength
		}
	}
	return -1, x
}
