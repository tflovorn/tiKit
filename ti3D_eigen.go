package main

import (
	"fmt"
	"flag"
	"os"
	"errors"
)

// temp - 1d matrix
type Matrix float64

// temp - 1d vector
type Vector float64

// Function that takes a 3D k-vector and returns a matrix result
type KToMatrixFn func(k Vector) Matrix

// Parse command-line arguments calcType, kpointsFileName, outFileName
func parseArgs() (string, string, string, error) {
	flag.Parse()
	args := flag.Args()
	if len(args) != 3 {
		err := errors.New("Usage: ti3D_eigen [8band|4band|mnk12] kpointsFileName outFileName")
		return "", "", "", err
	}
	return args[0], args[1], args[2], nil
}

// Parse kpoints file
func getKpoints(kpointsFile string) ([]Vector, error) {
	//TODO open file or error
	kpoints := []Vector{}
	//TODO read file (append to kpoints) or error
	return kpoints, nil
}

// Return a Hamiltonian function (k -> H_k) of the type specified
func HamiltonianFn(calcType string) (KToMatrixFn, error) {
	// valid calcTypes: "8band", "4band", "mnk12"
	// TODO implement (8band first)
	fn := func(k Vector) Matrix {
		return 0.0
	}
	return fn, nil
}

// Return the eigenvalues and corresponding eigenvectors of M
func EigenDecompose(M Matrix) ([]float64, []Vector) {
	//TODO
	return []float64{}, []Vector{}
}

// Write the output for one kpoint
func writeOutput(k Vector, eigenvals []float64, eigenkets []Vector, outFile *os.File) error {
	return nil
}

func main() {
	// command line arguments
	calcType, kpointsFileName, outFileName, err := parseArgs()
	if err != nil {
		fmt.Println(err)
		os.Exit(2)
	}
	// read input file
	kpoints, err := getKpoints(kpointsFileName)
	if err != nil {
		fmt.Println(err)
		os.Exit(2)
	}
	// get appropriate Hamiltonian
	H, err := HamiltonianFn(calcType)
	if err != nil {
		fmt.Println(err)
		os.Exit(2)
	}
	// open output file
	outFile, err := os.Create(outFileName)
	if err != nil {
		fmt.Println(err)
		os.Exit(2)
	}

	// iterate over kpoints (parallelize later)
	for _, k := range kpoints {
		Hk := H(k)
		eigenvals, eigenkets := EigenDecompose(Hk)
		err := writeOutput(k, eigenvals, eigenkets, outFile)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
	}

	// clean up
	err = outFile.Close()
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
