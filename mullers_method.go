

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



	//Represents: 3x^7 - 2x^6 + 20x^5 -4x^2 + x - 12
	baseEquation := [][]complex128{[]complex128{3, 7}, []complex128{-2, 6}, []complex128{20, 5}, []complex128{-4, 2}, []complex128{1, 1}, []complex128{-12, 0}}


	//get random xn xn - 1 and xn -2 values 

	//(the minuses are in the variable name nm1 = n minus 1.. nm2 = n minus 2... etc)

	//choose random seeds for the x values, between 0 and 100 

	var xn complex128

	xn = 1

	var xnm1 complex128

	xnm1 = 1

	var xnm2 complex128

	xnm2 = 1

	//x values must be unique
	for xn == xnm1 || xnm1 == xnm2 || xn == xnm2  {




		rand.Seed(time.Now().UnixNano())
		rxn := complex( float64(rand.Int()%7), float64(0))

		xn = rxn

		time.Sleep(time.Duration(2)*time.Millisecond)

		rand.Seed(time.Now().UnixNano())
		rxnm1 := complex(float64(rand.Int()%10), float64(0))

		xnm1 = rxnm1

		time.Sleep(time.Duration(2)*time.Millisecond)

		rand.Seed(time.Now().UnixNano())
		rxnm2 := complex(float64(rand.Int()%20), float64(0))

		xnm2 = rxnm2 
	


	}

	fmt.Println("Initial X vals", xn, xnm1, xnm2)

	yn := EvaluateEquationAtX(baseEquation, xn)
	ynm1 := EvaluateEquationAtX(baseEquation, xnm1)
	ynm2 := EvaluateEquationAtX(baseEquation, xnm2)


	topRow := []complex128{cmplx.Pow(xn, 2), xn, 1, yn}
	midRow := []complex128{cmplx.Pow(xnm1, 2), xnm1, 1, ynm1}
	botRow := []complex128{cmplx.Pow(xnm2, 2), xnm2, 1, ynm2}

	matrixOfSystem := Matrix{[][]complex128{topRow, midRow, botRow}}



	fmt.Println("Initial Y Vals", yn, ynm1, ynm2)


	// a, b, c := GetEquationForParabola(xn, xnm1, xnm2, yn, ynm1, ynm2)


	a, b, c := CalculateABCForParabolaEquation(matrixOfSystem)



	fmt.Println("a = ", a, "b = ", b, "c = ", c)

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



		fmt.Println("Initial Y Vals", yn, ynm1, ynm2)



		a, b, c = CalculateABCForParabolaEquation(matrixOfSystem)


		// a, b, c = GetEquationForParabola(xn, xnm1, xnm2, yn, ynm1, ynm2)

		fmt.Println("New X vals", xn, xnm1, xnm2)

		fmt.Println("New Y Vals", yn, ynm1, ynm2)


	}

	fmt.Println(xn, " is a root")

	fmt.Println("f(", xn, ") = ", EvaluateEquationAtX(baseEquation, xn))




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


func GetEquationForParabola(xn complex128, xnm1 complex128, xnm2 complex128, yn complex128, ynm1 complex128, ynm2 complex128) (complex128, complex128, complex128){

	IsSolution := false


	equation1 := []complex128{-1*cmplx.Pow(xn, 2),   -1*xn,   -1,  yn}
	equation2 := []complex128{-1*cmplx.Pow(xnm1, 2),   -1*xnm1,   -1,  ynm1}
	equation3 := []complex128{-1*cmplx.Pow(xnm2, 2),   -1*xnm2,   -1,  ynm2}

	systemMatrix := Matrix{[][]complex128{equation1, equation2, equation3}}


	var originalEquationsToSatisfy []Equation

	for i := 0; i < len(systemMatrix.Rows); i++ {
		originalEquationsToSatisfy = append(originalEquationsToSatisfy, Equation{systemMatrix.Rows[i]})
	}


	cleanCopyOriginal := CleanCopyMatrix(systemMatrix)

	newMatrix := systemMatrix

	repeatedFailureCount := 0

	for !IsSolution {

		// PrintMatrix(newMatrix)

		for i := 0; i < 30; i++ {

			rand.Seed(time.Now().UnixNano())
			//random number between rows 0 and 5
			randomNumber := rand.Int()%(len(newMatrix.Rows))

			newMatrix = ZeroOutAllButSpecificIndex(randomNumber, newMatrix)
			newMatrix = ZeroOutAllButSpecificIndex(randomNumber, newMatrix)
			newMatrix = ZeroOutAllButSpecificIndex(randomNumber, newMatrix)


			for j := 0; j < len(newMatrix.Rows); j++ {
				newMatrix = ZeroOutAllButSpecificIndex(j, newMatrix)
				newMatrix = ScaleEachRowInMatrixByLowestNumberInRow(newMatrix)
				
			}

		}




		for i := 0; i < 3; i++ {

			rand.Seed(time.Now().UnixNano())
			//random number between rows 0 and 5
			randomNumber := rand.Int()%(len(newMatrix.Rows))

			newMatrix = ZeroOutAllButSpecificIndex(randomNumber, newMatrix)

			newMatrix = GiveARandomRowSubstitutionRegardlessOfOutcome(newMatrix)

			for j := 0; j < len(newMatrix.Rows); j++ {
				newMatrix = ZeroOutAllButSpecificIndex(j, newMatrix)
				newMatrix = ScaleEachRowInMatrixByLowestNumberInRow(newMatrix)
			}

		}


		// PrintMatrix(newMatrix)

		


		// fmt.Println("looping ok")
		allRows1 := true

		for i := 0; i < len(newMatrix.Rows); i++ {
			if(NonZeroVariableCount(newMatrix.Rows[i]) != 1){
				allRows1 = false
			}
		}

		// fmt.Println("all rows are 1", allRows1)

		if(AllRowsAre1Or0(newMatrix) && !allRows1){
			
			repeatedFailureCount++

				if(repeatedFailureCount == 60){
					return complex(0, 0), complex(0, 0), complex(0, 0)
				}

			//if it was all 1's and 0's in the matrix 
			//(meaning not every variable got a solution, and the matrix cannot be further simplified)
			//reset to original matrix and try solving again
			// panic("resets at least oce")
			newMatrix = CleanCopyMatrix(cleanCopyOriginal)
		}else if(allRows1){
	
			//newMatrix = NullOutImaginaryParts(newMatrix)

			// fmt.Println("matrix with nulled")
			// PrintMatrix(newMatrix)

		//	fmt.Println("POTENTIAL SOLUTION FOUND")
		//	PrintMatrix(newMatrix)

			solutionGuess := make([]complex128, len(newMatrix.Rows))

			

			for i := 0; i < len(newMatrix.Rows); i++ {

				currentRow := newMatrix.Rows[i]

				
				for j := 0; j < (len(currentRow)-1); j++ {

					if(currentRow[j] != complex(0, 0)){
						solutionGuess[j] = currentRow[len(currentRow)-1]
					}

				}


			}

			fmt.Println("Guess", solutionGuess)

			areSame := TwoMatricesEvaluateTheSame(cleanCopyOriginal, solutionGuess)

			


				

			IsSolution = CheckOriginalEquationsAreSatisfied(solutionGuess, originalEquationsToSatisfy)


			if(IsSolution && areSame){

				fmt.Println("GUESS IS A SOLUTION!", IsSolution)

				for i := 0; i < len(solutionGuess); i++ {
					fmt.Println(VariableNameAlphabetIndex(i), " = ", solutionGuess[i]*-1)
				}
				

				return solutionGuess[0], solutionGuess[1], solutionGuess[2]

			}else{

				repeatedFailureCount++

				if(repeatedFailureCount == 60){
					return complex(0, 0), complex(0, 0), complex(0, 0)
				}
				//if it was all 1's in the matrix but not a solution,
				//reset to original matrix and try solving again
				newMatrix = CleanCopyMatrix(cleanCopyOriginal)
			}
		}



	}

	return complex128(-1), complex128(-1), complex128(-1)

	

}

func NullOutImaginaryParts(matrixInput Matrix) Matrix {

	rows := matrixInput.Rows 

	for i := 0; i < len(rows); i++ {

		currentRow := rows[i]

		complexZero := complex(0, 0)

		imaginaryPartToNull := imag(complexZero)

		indexOfNull := -1

		for j := 0; j < (len(currentRow)-1); j++ {
			if(currentRow[j] != complex(0, 0)){
				imaginaryPartToNull = imag(currentRow[j])
				indexOfNull = j
			}
		}

		if(indexOfNull == -1){
			panic("could not get non zero index for row NullOutImaginaryParts")
		}

		matrixInput.Rows[i][indexOfNull] = matrixInput.Rows[i][indexOfNull] - complex(0, imaginaryPartToNull)
		matrixInput.Rows[i][len(rows[i])-1] = matrixInput.Rows[i][len(rows[i])-1] - complex(0, imaginaryPartToNull) 

	}

	return matrixInput


}

func AllPossibleOutcomesTwoRows(row1 []complex128, row2 []complex128) ([][]complex128, int) {

	var worker []complex128

	var workedOn []complex128

	row1NonZeroCount := NonZeroVariableCount(row1)
	row2NonZeroCount := NonZeroVariableCount(row2)

	returnIntCodeToFlagWhichIsWorker := -1

	if(row1NonZeroCount > row2NonZeroCount){
		worker = row2
		workedOn = row1
		returnIntCodeToFlagWhichIsWorker = 1
	}else{
		worker = row1
		workedOn = row2
		returnIntCodeToFlagWhichIsWorker = 0
	}


	if(returnIntCodeToFlagWhichIsWorker == -1){
		panic("invalid int code AllPossibleOutcomesTwoRows()")
	}


	resultsToReturn := [][]complex128{}

	

	for i := 0; i < (len(row1)-1); i++ {

		if(worker[i] == complex(0, 0) || workedOn[i] == complex(0, 0)){
			continue
		}else{

			negativeWorkedOn := workedOn[i]*-1
			eliminationScalar := negativeWorkedOn/worker[i] 


	

			newResult := []complex128{}

			for j := 0; j < len(row1); j++ {
				workerValScaled := worker[j]*eliminationScalar
				workedOnVal := workedOn[j]

				

				newVal := workerValScaled + workedOnVal

				

				newResult = append(newResult, newVal)

			}

			resultsToReturn = append(resultsToReturn, newResult)

		}

	}

	return resultsToReturn, returnIntCodeToFlagWhichIsWorker

}

func NonZeroVariableCount(row []complex128) int {

	nonZeroCount := 0

	for i := 0; i < (len(row) - 1); i++ {

		if(real(row[i]) != 0){
			nonZeroCount++
		}

	}

	return nonZeroCount



}

func ZeroOutAllButSpecificIndex(indexNotToZero int, matrixInput Matrix) Matrix{

	c1 := 0

	c2 := 1 

	foundSolution := false

	rows := matrixInput.Rows

	currentRow := rows[c1]

	beginningToEndNotEvenOneGoodChange := true

	for !foundSolution {

		//fmt.Println("good change", beginningToEndNotEvenOneGoodChange)


		validMovingCursorFound := false

		for !validMovingCursorFound {

//		fmt.Println("stuck in loop")	
//			fmt.Println("loop c1", c1, "c2", c2)

			if(c2 != len(rows) && c2 != c1){
				validMovingCursorFound = true
				break
			}else{
				c2++ 
			}

			if(c2 >= len(rows)){
				c2 = 0
				c1++
				if(c1 == len(rows)){
					if(!beginningToEndNotEvenOneGoodChange){
						c1 = 0 
						currentRow = rows[c1]
						c2 = 1
						beginningToEndNotEvenOneGoodChange = true
		//				PrintMatrix(matrixInput)
					}else{
						//PrintMatrix(matrixInput)



						//panic("out of bounds3")
						return matrixInput
					}
				}else{
					currentRow = rows[c1]
				}
			}

		}

		nonZeroCountCurrentRow := NonZeroVariableCount(currentRow)

		currentMoving := rows[c2]

		outcomes, intCode := AllPossibleOutcomesTwoRows(currentRow, currentMoving)


		if(intCode != 1 || len(outcomes) == 0){
			c2++
			continue

		}

		bestOutcome := []complex128{}
		bestNonZeroCount := nonZeroCountCurrentRow

		
		for i := 0; i < len(outcomes); i++ {
			if(!(RowOfInterestIsZero(indexNotToZero, outcomes[i])) ){
				nonZeroCountCurrentOutcome := NonZeroVariableCount(outcomes[i])
				
				if(nonZeroCountCurrentOutcome < bestNonZeroCount){
					bestNonZeroCount = nonZeroCountCurrentOutcome
					bestOutcome = outcomes[i]
					beginningToEndNotEvenOneGoodChange = false
				}
			}
		}


		if(len(bestOutcome) == 0){
			c2++
			continue
		}else{
			currentRow = bestOutcome
			matrixInput.Rows[c1] = bestOutcome
		}

		

	}


	return matrixInput

}



func TwoMatricesEvaluateTheSame(matrix1 Matrix, testVals []complex128) bool {

	if(len(testVals) != (len(matrix1.Rows[0]) - 1)){
		panic("incorrent number of test vals TwoMatricesEvaluateTheSame()")
	}

	// fmt.Println("matrix to evaluate")
	// PrintMatrix(matrix1)

	fmt.Println("test vals", testVals)

	for i := 0; i < len(matrix1.Rows); i++ {

		summationRow := complex(0, 0)

		currentRow := matrix1.Rows[i]

		for j := 0; j < len(currentRow)-1; j++ {

			summationRow = summationRow + currentRow[j] * testVals[j]

		}

		if(!aboutEquals(summationRow, currentRow[len(currentRow) - 1])) {
			fmt.Println(summationRow, currentRow[len(currentRow) - 1])
			panic("matrix does not evaluate the same TwoMatricesEvaluateTheSame()")
		}


	}

	return true





}

func CalculateABCForParabolaEquation(matrixInput Matrix) (complex128, complex128, complex128) {

	if(len(matrixInput.Rows[0]) != 4){
		panic("error, equation could not define a parabola length is not 4 for matrix")
	}

	fmt.Println("Input matrix")
	PrintMatrix(matrixInput)

	solutionColumn := []complex128{matrixInput.Rows[0][3], matrixInput.Rows[1][3], matrixInput.Rows[2][3]}

	fmt.Println("solution column")
	fmt.Println(solutionColumn)

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


func RowOfInterestIsZero(indexToZero int, row []complex128) bool {

	if(row[indexToZero] == complex(0, 0)){
		return true
	}else{
		return false
	}


}




func PrintMatrix(matrixInput Matrix) {

	fmt.Println("PRINTING MATRIX")

	for i := 0; i < len(matrixInput.Rows); i++ {
		fmt.Println(matrixInput.Rows[i])
	}
}

func ScaleEachRowInMatrixByLowestNumberInRow(matrixInput Matrix) Matrix {


	fmt.Println("pre scaled")
	PrintMatrix(matrixInput)

	rows := matrixInput.Rows

	
	for i := 0; i < len(rows); i++ {

		currentRow := rows[i]

		lowestMultiplier := float64(0)

		for j := 0; j < (len(currentRow)-1); j++ {

			if(real(currentRow[j]) != 0 && lowestMultiplier == 0){
				lowestMultiplier = real(currentRow[j])
			}else if(real(currentRow[j]) != 0 && math.Abs(real(currentRow[j])) < math.Abs(lowestMultiplier)) {
				lowestMultiplier = real(currentRow[j])
			}

		}

		if(lowestMultiplier != 0){
			for k := 0; k < len(matrixInput.Rows[i]); k++ {
				matrixInput.Rows[i][k] = matrixInput.Rows[i][k]/(complex(lowestMultiplier, 0)) 
			}
		}		

	}

	fmt.Println("post scaled")
	PrintMatrix(matrixInput)

	return matrixInput

}



// func ScaleEachRowInMatrixByLowestNumberInRow(matrixInput Matrix) Matrix {

// 	rows := matrixInput.Rows

	
// 	for i := 0; i < len(rows); i++ {

// 		currentRow := rows[i]

// 		newScaled := []complex128{}

// 		lowestNumberForRowAbs := float64(0)

// 		lowestNumberForRowActual := float64(0)

// 		fmt.Println("pre scale")
// 		PrintMatrix(matrixInput)


// 		for j := 0; j < (len(currentRow)-1); j++ {
				
// 			number := real(currentRow[j])

// 			if(lowestNumberForRowAbs == 0 && number != 0){
// 				lowestNumberForRowAbs = number
// 				lowestNumberForRowActual = real(currentRow[j])
// 			}else if(lowestNumberForRowAbs != 0 && number != 0){
// 				if(number < lowestNumberForRowAbs){
// 					lowestNumberForRowAbs = number
// 					lowestNumberForRowActual = real(currentRow[j])
// 				}
// 			}
// 		}


// 		// for j := 0; j < (len(currentRow)-1); j++ {
				
// 		// 	number := cmplx.Abs(currentRow[j])

// 		// 	if(lowestNumberForRowAbs == complex(1, 0) && number != complex(1, 0)){
// 		// 		lowestNumberForRowAbs = complex(number, 0)
// 		// 		lowestNumberForRowActual = currentRow[j]
// 		// 	}else if(lowestNumberForRowAbs != 0 && number != 0){
// 		// 		if(cmplx.Abs(complex(number, 0)) < cmplx.Abs(lowestNumberForRowAbs)){
// 		// 			lowestNumberForRowAbs = complex(number, 0)
// 		// 			lowestNumberForRowActual = currentRow[j]
// 		// 		}
// 		// 	}
// 		// }

// 		fmt.Println("post scale")
// 		PrintMatrix(matrixInput)

// 		if(lowestNumberForRowActual != 0){

// 			for k := 0; k < len(currentRow); k++ {

// 				// if(real(lowestNumberForRowActual) == -1) {
// 				// 	newScaled = append(newScaled, currentRow[k]/ complex(real(-1*lowestNumberForRowActual), 0))
// 				// }else{
// 				// 	newScaled = append(newScaled, currentRow[k]/ complex(real(lowestNumberForRowActual), 0))
// 				// }
				
// 				// newScaled = append(newScaled, currentRow[k]/ complex(lowestNumberForRowActual, 0) )

// 				newScaled = append(newScaled,  ((currentRow[k]/complex(lowestNumberForRowActual, 0))/-1))

// 			}

// 			matrixInput.Rows[i] = newScaled
// 		}


// 		fmt.Println("post scale 2")
// 		PrintMatrix(matrixInput)

		
// 		realMultiplier := float64(0)

// 		for l := 0; l < (len(matrixInput.Rows[i]) - 1); l++ {
			
// 			if(matrixInput.Rows[i][l] != complex(0, 0) && realMultiplier == 0){
// 				realMultiplier = real(matrixInput.Rows[i][l])
// 			}else if(matrixInput.Rows[i][l] != complex(0, 0) && math.Abs(real(matrixInput.Rows[i][l])) < math.Abs(realMultiplier)){
// 				realMultiplier = real(matrixInput.Rows[i][l])
// 			}
// 		}

// 		if(realMultiplier != 0){
			
// 			for m := 0; m < (len(matrixInput.Rows[i])); m++ {
// 				matrixInput.Rows[i][m] = matrixInput.Rows[i][m]/(complex(realMultiplier, 0))
// 			}

// 		}
		

// 		fmt.Println("post scale 3")
// 		PrintMatrix(matrixInput)

// 		panic("exit")

// 	}

// 	// fmt.Println("post scale")
// 	// PrintMatrix(matrixInput)

// 	return matrixInput

// }



func GiveARandomRowSubstitutionRegardlessOfOutcome(matrixInput Matrix) Matrix {

	if(AllRowsAre1Or0(matrixInput)){
		return matrixInput
	}

	substitutionMade := false

	counter := 0

	for !substitutionMade {

		

		rand.Seed(time.Now().UnixNano())

		//random number between rows 0 and 5
		randomNumber1 := rand.Int()%len(matrixInput.Rows)

		rand.Seed(time.Now().UnixNano())
		randomNumber2 := rand.Int()%len(matrixInput.Rows)

		for randomNumber1 == randomNumber2{
			rand.Seed(time.Now().UnixNano())
			randomNumber2 = rand.Int()%len(matrixInput.Rows)
		}

		rows := matrixInput.Rows

		outcomes, worker := AllPossibleOutcomesTwoRows(rows[randomNumber1], rows[randomNumber2])

		if(len(outcomes) != 0){
			if(worker == 0){
				matrixInput.Rows[randomNumber1] = outcomes[0] 
			}else{
				matrixInput.Rows[randomNumber2] = outcomes[0] 
			}

			substitutionMade = true		
		}


		counter++

		if(counter == 1000){
			return matrixInput
		}

	}


	return matrixInput


}




func CheckOriginalEquationsAreSatisfied(guessToSolution []complex128, originalEquations []Equation) bool {


	for i := 0; i < len(originalEquations); i++ {
		
		equation := originalEquations[i]

		equationResult := SubstituteValuesIntoEquationFinalGuessReturnSummation(equation, guessToSolution)

		if(len(guessToSolution) != (len(originalEquations[i].Items)-1) ){
			panic("length error")
		}

		if(!aboutEquals(equationResult, equation.Items[len(equation.Items)-1])){
			fmt.Println(equation)
			fmt.Println("didnt satisfy equation", (i+1), " ", equationResult, " ", equation.Items[len(equation.Items)-1])
			return false
		}		
	}

	return true

}






func SubstituteValuesIntoEquationFinalGuessReturnSummation(equationReceiver Equation, itemsToSubstitute []complex128) complex128 {

	substitutedFloatSlice := []complex128{}

	fmt.Println(itemsToSubstitute)

	//ok to index both at i since length has been checked
	for i := 0; i < len(equationReceiver.Items) - 1; i++ {
		substitutedFloatSlice = append(substitutedFloatSlice, itemsToSubstitute[i] * equationReceiver.Items[i])
	}

	summationOfSubstitutedFloatSlice := GetSummationcomplex128Slice(substitutedFloatSlice)

	return summationOfSubstitutedFloatSlice

}



func GetSummationcomplex128Slice(items []complex128) complex128 {

	summation := complex128(0)

	for i := 0; i < len(items); i++ {

		summation = summation + items[i]

	}

	return summation




}


func aboutEquals(checkVal complex128, result complex128) bool {
	
fmt.Println(checkVal, result)

	differenceReal  := math.Abs(real(checkVal) - real(result))

	fmt.Println(differenceReal)

	differenceImag := math.Abs(imag(checkVal) - imag(result))

	fmt.Println(differenceImag)

	if(differenceReal  <  float64(0.001) && differenceImag  <  float64(0.001)) {
		return true
	}else{
		return false
	}
}





func VariableNameAlphabetIndex(index int) string {

	if(index < 0){
		panic("impossible index for alphabet")
	}


	switch index {

		case 0:
			return "A"
		case 1:
			return "B"
		case 2:
			return "C"
		case 3:
			return "D"
		case 4:
			return "E"
		case 5:
			return "F"
		case 6:
			return "G"
		case 7:
			return "H"
		case 8:
			return "I"
		case 9:
			return "J"
		case 10:
			return "K"
		case 11:
			return "L"
		case 12:
			return "M"
		case 13:
			return "N"
		case 14:
			return "O"
		case 15:
			return "P"
		case 16:
			return "Q"
		case 17:
			return "R"
		case 18:
			return "S"
		case 19:
			return "T"
		case 20:
			return "U"
		case 21:
			return "V"
		case 22:
			return "W"
		case 23:
			return "X"
		case 24:
			return "Y"
		case 25:
			return "Z"
		default:


			newNumber := index - 25

			return "A"+ VariableNameAlphabetIndex(newNumber)
			
			
	}

	return "-1*z*error"



}


func AllRowsAre1Or0(matrixInput Matrix) bool {

	

	for i := 0; i < len(matrixInput.Rows); i++ {

		currentRow := matrixInput.Rows[i]

		for j := 0; j < len(currentRow)-1; j++ {
			if(currentRow[j] != complex(0, 0) && currentRow[j] != complex(1, 0)){
				
				return false
			}
		}

	}

	

	return true
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




































