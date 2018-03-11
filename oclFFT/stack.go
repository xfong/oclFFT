package oclFFT

var oclFFTStack []plan

func DestroyPlan(idx int) {
	targPlan := oclFFTStack[idx]
	targPlan.destroy()
	oclFFTStack = append(oclFFTStack[:idx], oclFFTStack[idx+1:]...)
}

func TearDownOclFFT() {
	for idx := range oclFFTStack {
		DestroyPlan(idx)
	}
	tmp := new([]plan)
	oclFFTStack = *tmp
}
