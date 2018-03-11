package main

import (
	"flag"
	"fmt"
	"github.com/source/oclFFT/oclFFT"
)

var (
	Flag_length = flag.Int("count", -1, "fft length")
)

func main() {
	flag.Parse()
	testVar := oclFFT.NewOclFFTPlanBase()
	testVar.SetTransformSize(*Flag_length)
	fmt.Println("Requested transform length:\t", testVar.GetTransformSize())
	fmt.Println("Kernel lengths:\t\t", oclFFT.FFTKernels)
	fmt.Println("Number of stages:\t", testVar.GetStageCounts())
	fmt.Println("KW entering stage:\t", testVar.GetEntryKW())
}
