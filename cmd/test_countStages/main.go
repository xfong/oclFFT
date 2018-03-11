package main

import (
	"flag"
	"fmt"
)

var (
	Flag_length = flag.Int("count", -1, "fft length")
	initKW = [...]int {1, 1, 1, 1, 1, 1, 1}
)

func main() {
	flag.Parse()
	var tmpCnt [6]int
	fft_size := *Flag_length
	if (fft_size < 2) {
		fmt.Println("Cannot use count less than 2!")
		
	} else {
		fmt.Println("Initial fft length:\t", fft_size)
		tmpCnt, initKW = countStages(fft_size)
	}
	fmt.Println("Entry fft lengths:\t", initKW)
	fmt.Println("FFT Kernel lengths:\t", FFTKernels)
	fmt.Println("Stage counts:\t\t", tmpCnt)
}
