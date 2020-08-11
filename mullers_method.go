

package main


import (
	
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"time"

)

type Matrix struct {

	Rows [][]complex128

}

type Equation struct {

    Items []complex128 `json:"eqtn"`

}



func main() {


	//TODO.... YOUR EQUATION GOES HERE
	//Represents: 3x^7 - 2x^6 + 20x^5 -4x^2 + x - 12
	baseEquation := [][]complex128{[]complex128{3, 7}, []complex128{-2, 6}, []complex128{20, 5}, []complex128{-4, 2}, []complex128{1, 1}, []complex128{-12, 0}}


	//TODO YOUR HIGHEST POWER GOES HERE
	highestPowerOfEquation := 7

	//technically do not need 2 variables but it is good to point out the highest power is total root count
	totalRoots := highestPowerOfEquation

	//get random xn xn - 1 and xn -2 values 

	//(the minuses are in the variable name nm1 = n minus 1.. nm2 = n minus 2... etc)

	//choose random seeds for the x values, between 0 and 100 

	var xn complex128

	xn = 1

	var xnm1 complex128

	xnm1 = 1

	var xnm2 complex128

	xnm2 = 1

	rootsFound := []complex128{}


	for len(rootsFound) < totalRoots-1 {

		xn = 0
	
		xnm1 = 0

		xnm2 = 0
		

		//x values must be unique
		for xn == xnm1 || xnm1 == xnm2 || xn == xnm2  {

			//pseudo randomally also need to check negative x inputs

			addANegativeTerm0 := false

			addANegativeTerm1 := false

			addANegativeTerm2 := false

			rand.Seed(time.Now().UnixNano())
			negative0 := rand.Int()%100
			rand.Seed(time.Now().UnixNano())
			negative1 := rand.Int()%100
			rand.Seed(time.Now().UnixNano())
			negative2 := rand.Int()%100

			if(negative0%2 == 0){
				addANegativeTerm0 = true
			}
			if(negative1%2 == 0){
				addANegativeTerm1 = true
			}
			if(negative2%2 == 0){
				addANegativeTerm2 = true
			}


			rand.Seed(time.Now().UnixNano())
			rxn := complex( float64(rand.Int()%100), float64(0))

			xn = rxn

			if(addANegativeTerm0){
				xn = xn * - 1
			}

			time.Sleep(time.Duration(2)*time.Millisecond)

			rand.Seed(time.Now().UnixNano())
			rxnm1 := complex(float64(rand.Int()%100), float64(0))

			xnm1 = rxnm1

			if(addANegativeTerm1){
				xnm1 = xnm1 * -1
			}

			time.Sleep(time.Duration(2)*time.Millisecond)

			rand.Seed(time.Now().UnixNano())
			rxnm2 := complex(float64(rand.Int()%100), float64(0))


			xnm2 = rxnm2 


			if(addANegativeTerm2){
				xnm2 = xnm2 * -1
			}		


		}

		// fmt.Println("Initial X vals", xn, xnm1, xnm2)

		yn := EvaluateEquationAtX(baseEquation, xn)
		ynm1 := EvaluateEquationAtX(baseEquation, xnm1)
		ynm2 := EvaluateEquationAtX(baseEquation, xnm2)


		topRow := []complex128{cmplx.Pow(xn, 2), xn, 1, yn}
		midRow := []complex128{cmplx.Pow(xnm1, 2), xnm1, 1, ynm1}
		botRow := []complex128{cmplx.Pow(xnm2, 2), xnm2, 1, ynm2}

		matrixOfSystem := Matrix{[][]complex128{topRow, midRow, botRow}}



		// fmt.Println("Initial Y Vals", yn, ynm1, ynm2)


		// a, b, c := GetEquationForParabola(xn, xnm1, xnm2, yn, ynm1, ynm2)


		a, b, c := CalculateABCForParabolaEquation(matrixOfSystem)



		// fmt.Println("a = ", a, "b = ", b, "c = ", c)



		//for a root not yet found
		for  !aboutEquals(EvaluateEquationAtX(baseEquation, xn), complex(0, 0)) {

			q := EvaluateQ(xn, xnm1, xnm2)

			BigA := EvaluateBigA(q, a, b, c, xn, xnm1, xnm2)

			BigB := EvaluateBigB(q, a, b, c, xn, xnm1, xnm2)

			BigC := EvaluateBigC(q, a, b, c, xn)

			nextTerm := EvaluateNextTerm(BigA, BigB, BigC, xn, xnm1)

			xnm2 = xnm1
			xnm1 = xn
			xn = nextTerm

			yn = EvaluateEquationAtX(baseEquation, xn)
			ynm1 = EvaluateEquationAtX(baseEquation, xnm1)
			ynm2 = EvaluateEquationAtX(baseEquation, xnm2)

			topRow := []complex128{cmplx.Pow(xn, 2), xn, 1, yn}
			midRow := []complex128{cmplx.Pow(xnm1, 2), xnm1, 1, ynm1}
			botRow := []complex128{cmplx.Pow(xnm2, 2), xnm2, 1, ynm2}

			matrixOfSystem := Matrix{[][]complex128{topRow, midRow, botRow}}



			// fmt.Println("Initial Y Vals", yn, ynm1, ynm2)



			a, b, c = CalculateABCForParabolaEquation(matrixOfSystem)


			// a, b, c = GetEquationForParabola(xn, xnm1, xnm2, yn, ynm1, ynm2)

			// fmt.Println("New X vals", xn, xnm1, xnm2)

			// fmt.Println("New Y Vals", yn, ynm1, ynm2)


		}

		// fmt.Println(xn, " is a root")

		// fmt.Println("f(", xn, ") = ", EvaluateEquationAtX(baseEquation, xn))

		rootsFound = CleanRootAddIfNotDuplicate(rootsFound, xn)

	// 	for i := 0; i < len(rootsFound); i++ {
	// 	fmt.Println(rootsFound[i])
	// }

	}


	for i := 0; i < len(rootsFound); i++ {
		fmt.Println(rootsFound[i])
	}



}


func EvaluateNextTerm(A complex128, B complex128, C complex128, xn complex128, xnm1 complex128) complex128 {

	numerator := 2 * C

	var finalDenominator complex128

	BSquared := cmplx.Pow(B, complex(2, 0))

	ACTimesFour := 4*A*C

	denominatorOption1 := B + cmplx.Pow((BSquared - ACTimesFour), complex(0.5, 0))
	denominatorOption2 := B - cmplx.Pow((BSquared - ACTimesFour), complex(0.5, 0))

	if(cmplx.Abs(denominatorOption1) > cmplx.Abs(denominatorOption2)){
		finalDenominator = denominatorOption1
	}else{
		finalDenominator = denominatorOption2
	}

	fractionEvaluated := numerator/finalDenominator

	multiplicationTermEvaluated := (xn - xnm1)*(fractionEvaluated)

	newTerm := xn - multiplicationTermEvaluated

	return newTerm


}

func EvaluateBigA(q complex128, a complex128, b complex128, c complex128, xn complex128, xnm1 complex128, xnm2 complex128) complex128 {

	pxn := EvaluateParabolaAtX(a, b, c, xn)
	pxnm1 := EvaluateParabolaAtX(a, b, c, xnm1)
	pxnm2 := EvaluateParabolaAtX(a, b, c, xnm2)

	firstTerm := (q) * (pxn) 
	secondTerm := (-1)*(q)*(1 + q)*(pxnm1)
	thirdTerm := (cmplx.Pow(q, 2))*(pxnm2)

	return firstTerm + secondTerm + thirdTerm

}

func EvaluateBigB(q complex128, a complex128, b complex128, c complex128, xn complex128, xnm1 complex128, xnm2 complex128) complex128 {

	pxn := EvaluateParabolaAtX(a, b, c, xn)
	pxnm1 := EvaluateParabolaAtX(a, b, c, xnm1)
	pxnm2 := EvaluateParabolaAtX(a, b, c, xnm2)

	firstTerm := ((2*q) + 1) * (pxn) 
	secondTerm := (-1)*(cmplx.Pow((1 + q), 2))*(pxnm1)
	thirdTerm := (cmplx.Pow(q, 2))* (pxnm2)

	return firstTerm + secondTerm + thirdTerm
	
}

func EvaluateBigC(q complex128, a complex128, b complex128, c complex128, xn complex128) complex128 {
	
	pxn := EvaluateParabolaAtX(a, b, c, xn)
	

	return (1 + q) * (pxn)

}

func EvaluateQ(xn complex128, xnm1 complex128, xnm2 complex128) complex128 {

	return (xn - xnm1)/(xnm1 - xnm2)

}



func EvaluateParabolaAtX(a complex128, b complex128, c complex128, x complex128) complex128 {

	aTerm := a * cmplx.Pow(x, 2)
	bTerm := b * x 
	cTerm := c

	return aTerm + bTerm + cTerm



} 

func EvaluateEquationAtX(equationSlice [][]complex128, x complex128) complex128 {

	summation := complex128(0)

	for i := 0; i < len(equationSlice); i++ {

		multiplier := equationSlice[i][0]
		exponent := equationSlice[i][1]

		summation = summation + (cmplx.Pow(x, exponent) * multiplier)

	}

	return summation



}


func CalculateABCForParabolaEquation(matrixInput Matrix) (complex128, complex128, complex128) {

	if(len(matrixInput.Rows[0]) != 4){
		panic("error, equation could not define a parabola length is not 4 for matrix")
	}

	// fmt.Println("Input matrix")
	// PrintMatrix(matrixInput)

	solutionColumn := []complex128{matrixInput.Rows[0][3], matrixInput.Rows[1][3], matrixInput.Rows[2][3]}

	// fmt.Println("solution column")
	// fmt.Println(solutionColumn)

	aDet := CalculateDeterminantOfMatrix(SwapSolutionColumnToXColumn(matrixInput, solutionColumn))
	

	bDet := CalculateDeterminantOfMatrix(SwapSolutionColumnToYColumn(matrixInput, solutionColumn))
	cDet := CalculateDeterminantOfMatrix(SwapSolutionColumnToZColumn(matrixInput, solutionColumn))

	mainDet := CalculateDeterminantOfMatrix(matrixInput)

//	panic("err")

	return (aDet/mainDet), (bDet/mainDet), (cDet/mainDet)


}


func SwapSolutionColumnToZColumn(matrixInput Matrix, solutionColumn []complex128) Matrix {

	if(len(solutionColumn) != 3 || len(matrixInput.Rows) != 3){
		panic("err solution column or matrix input != 3 SwapSolutionColumnToZColumn()")
	}  

	cleanCopyToReturn := CleanCopyMatrix(matrixInput)

	cleanCopyToReturn.Rows[0][2] = solutionColumn[0]
	cleanCopyToReturn.Rows[1][2] = solutionColumn[1]
	cleanCopyToReturn.Rows[2][2] = solutionColumn[2]

	return cleanCopyToReturn

}

func SwapSolutionColumnToYColumn(matrixInput Matrix, solutionColumn []complex128) Matrix {

	if(len(solutionColumn) != 3 || len(matrixInput.Rows) != 3){
		panic("err solution column or matrix input != 3 SwapSolutionColumnToYColumn()")
	}  

	cleanCopyToReturn := CleanCopyMatrix(matrixInput)

	cleanCopyToReturn.Rows[0][1] = solutionColumn[0]
	cleanCopyToReturn.Rows[1][1] = solutionColumn[1]
	cleanCopyToReturn.Rows[2][1] = solutionColumn[2]

	return cleanCopyToReturn

}


func SwapSolutionColumnToXColumn(matrixInput Matrix, solutionColumn []complex128) Matrix {

	if(len(solutionColumn) != 3 || len(matrixInput.Rows) != 3){
		panic("err solution column or matrix input != 3 SwapSolutionColumnToXColumn()")
	}  

	cleanCopyToReturn := CleanCopyMatrix(matrixInput)

	cleanCopyToReturn.Rows[0][0] = solutionColumn[0]
	cleanCopyToReturn.Rows[1][0] = solutionColumn[1]
	cleanCopyToReturn.Rows[2][0] = solutionColumn[2]

	return cleanCopyToReturn

}

//this is the scale factor of the 3D transform
func CalculateDeterminantOfMatrix(matrixInput Matrix) complex128{


	if(len(matrixInput.Rows) != 3){
		panic("not a 3D transform CalculateDeterminantOfMatrix()")
	}

	rows := matrixInput.Rows

	topRow := rows[0]
	middleRow := rows[1]
	bottomRow := rows[2]

	summation := complex(0, 0)

	for i := 0; i < 3; i++ {

		multiplier := topRow[i]



		if(i == 0){

			positiveTerm := middleRow[1]*bottomRow[2]
			negativeTerm := (middleRow[2]*bottomRow[1])*-1

			summation = summation+(multiplier*(positiveTerm + negativeTerm))



		}else if(i == 1){

			multiplier := multiplier * -1

			positiveTerm := middleRow[0]*bottomRow[2]
			negativeTerm := (middleRow[2]*bottomRow[0])*-1

			summation = summation+(multiplier*(positiveTerm + negativeTerm))


		}else if(i == 2){
			positiveTerm := middleRow[0]*bottomRow[1]
			negativeTerm := (middleRow[1]*bottomRow[0])*-1

			summation = summation+(multiplier*(positiveTerm + negativeTerm))
		}

	}

	return summation


}




func PrintMatrix(matrixInput Matrix) {

	fmt.Println("PRINTING MATRIX")

	for i := 0; i < len(matrixInput.Rows); i++ {
		fmt.Println(matrixInput.Rows[i])
	}
}







func aboutEquals(checkVal complex128, result complex128) bool {
	
	//fmt.Println(checkVal, result)

	differenceReal  := math.Abs(real(checkVal) - real(result))

	//fmt.Println(differenceReal)

	differenceImag := math.Abs(imag(checkVal) - imag(result))

	//fmt.Println(differenceImag)

	if(differenceReal  <  float64(0.001) && differenceImag  <  float64(0.001)) {
		return true
	}else{
		return false
	}
}





func CleanCopyMatrix(matrixInput Matrix) Matrix {

	rowsToCopy := matrixInput.Rows

	copiedReturnRows := [][]complex128{}

	for i := 0; i < len(rowsToCopy); i++ {

		newRow := make([]complex128, len(rowsToCopy[i]))

		itemsCopied := copy(newRow, rowsToCopy[i])

		if(itemsCopied != len(rowsToCopy[i])){
			panic("error copying CleanCopyMatrix()")
		}

		copiedReturnRows = append(copiedReturnRows, newRow)

	}

	if(len(copiedReturnRows) != len(rowsToCopy)){
		panic("not all rows copied CleanCopyMatrix()")
	}

	return Matrix{copiedReturnRows}


}



func CleanRootAddIfNotDuplicate(solutions []complex128, newSolution complex128) []complex128 {


	newSolutionCleaned := newSolution

	//check if the real term is 0
	if(aboutEquals(complex(real(newSolution), 0), complex(0, 0))) {
		newSolutionCleaned = complex(0, imag(newSolutionCleaned))
	}

	//check if the imaginary term is 0
	if(aboutEquals(complex(0, imag(newSolution)), complex(0, 0))) {
		newSolutionCleaned = complex(real(newSolutionCleaned), 0)
	}


	solutionIsDuplicate := false

	for i := 0; i < len(solutions); i++ {
		
		if(aboutEquals(newSolutionCleaned, solutions[i]) || aboutEquals(solutions[i], complex(real(newSolutionCleaned), (imag(newSolutionCleaned) * -1)))) {

			
			solutionIsDuplicate = true
			break
		}
	}

	if(!solutionIsDuplicate){

		solutions = append(solutions, newSolutionCleaned)

		//add both the solution and it's conjugate if the imaginary part is not null
		if(imag(newSolutionCleaned) != 0){
			solutions = append(solutions, complex(real(newSolutionCleaned), (imag(newSolutionCleaned) * -1)))
		}

	}


	return solutions

}

































