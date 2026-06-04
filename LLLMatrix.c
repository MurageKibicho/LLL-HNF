#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
//clear && gcc LLLMatrix.c -lm -o m.o && ./m.o
//Double vector that permits negative indices
double *NegativeIndexArray(long nl, long nh)
{
	int vectorEnd = 1;
	double *v = calloc((size_t) ((nh - nl + 1 + vectorEnd)), sizeof(double));
	assert(v != NULL);
	return v - nl + vectorEnd;
}

void FreeNegativeIndexArray(long nl, long nh,double *v)
{
	int vectorEnd = 1;
	free((char*)(v + nl - vectorEnd));
}

double **Matrix2D(int rows, int cols)
{
	double **matrix2D = malloc(rows * sizeof(double*));
	for(int i = 0; i < rows; i++)
	{
		matrix2D[i] = calloc(cols, sizeof(double));
	}
	return matrix2D;
}

void FreeMatrix2D(int rows, double **matrix2D)
{
	for(int i = 0; i < rows; i++)
	{
		free(matrix2D[i]);
	}
	free(matrix2D);
}

void PrintAugmentedMatrix(double **augmentedMatrix, int rows, int cols, int inputCols)
{
    // Calculate column widths for pretty printing
    int col_widths[cols];
    for(int j = 0; j < cols; j++) {
        col_widths[j] = 1;
        for(int i = 0; i < rows; i++) {
            char buf[32];
            snprintf(buf, sizeof(buf), "%.0f", fabs(augmentedMatrix[i][j]));
            int len = strlen(buf);
            if(len > col_widths[j]) col_widths[j] = len;
        }
        // Add space for sign
        col_widths[j] += 1;
    }
    
    // Print column headers
    printf("Row\\Col ");
    for(int j = 0; j < cols; j++) {
        if(j == inputCols + 1) printf(" | ");
        else if(j == cols-1) printf(" |e ");
        else printf(" ");
        
        printf("%*d", col_widths[j], j);
    }
    printf("\n");
    
    // Print separator line
    printf("--------");
    for(int j = 0; j < cols; j++) {
        printf("-");
        for(int w = 0; w < col_widths[j]; w++) printf("-");
        if(j == inputCols) printf("--|--");
        else if(j == cols-1) printf("--|--");
    }
    printf("\n");
    
    // Print matrix rows
    for(int i = 0; i < rows; i++) {
        printf("Row %2d: ", i);
        
        for(int j = 0; j < cols; j++) {
            // Add separators
            if(j == inputCols + 1) printf(" | ");
            else if(j == cols-1) printf(" | ");
            else printf(" ");
            
            // Print value
            double val = augmentedMatrix[i][j];
            if(fabs(val) < 1e-10) val = 0.0; // Treat near-zero as 0
            
            printf("%*.0f", col_widths[j], val);
        }
        
        // Add row description
        printf("  ");
        if(i < inputCols) {
            printf("(I row %d)", i);
        } else if(i == inputCols) {
            printf("(b row)");
        } else if(i == rows-1) {
            printf("(e_m row)");
        }
        
        printf("\n");
    }
}

void LLLMinus(int previousRow, int inputRows, int inputCols, double **augmentedMatrix, double **gramLambdaSquareMatrix)
{
	for(int col = 0; col <= inputRows + inputCols + 1; col++)
	{
		augmentedMatrix[previousRow][col] = -augmentedMatrix[previousRow][col];
	}
	for(int row = 1; row <= inputCols; row++)
	{
		for(int col = 0; col <= row - 1; col++)
		{
			if(row == previousRow || col == previousRow)
			{
				gramLambdaSquareMatrix[row][col] = -gramLambdaSquareMatrix[row][col];
			}
		}
	}
}

void LLLReduce(int *pivot0, int *pivot1,int previousRow, int currentRow, int augmentedRows, int augmentedCols, int inputRows, int inputCols, double **augmentedMatrix, double **gramLambdaSquareMatrix, double *gramDeterminants)
{
	int noPivotFound = inputRows + inputCols + 2;int pivotColumn0 = noPivotFound;int pivotColumn1 = noPivotFound;
	//Find Pivots
	for(int col = inputCols + 1; col <= inputRows+inputCols+1; col++)
	{
		if(augmentedMatrix[previousRow][col] != 0)
		{
			pivotColumn0 = col;
			break;	
		}
	}
	for(int col = inputCols + 1; col <= inputRows+inputCols+1; col++)
	{
		if(augmentedMatrix[currentRow][col] != 0)
		{
			pivotColumn1 = col;
			break;	
		}
	}
	
	double quotient = 0;
	if(pivotColumn0 < noPivotFound)
	{
		//Handle negatives
		if(augmentedMatrix[previousRow][pivotColumn0] < 0)
		{
			//Minus
			LLLMinus(previousRow, inputRows, inputCols, augmentedMatrix, gramLambdaSquareMatrix);
		}
		quotient = floor(augmentedMatrix[currentRow][pivotColumn0] / augmentedMatrix[previousRow][pivotColumn0]); 
	}
	else
	{
		double currentLambda = gramLambdaSquareMatrix[currentRow][previousRow];
		if(2 * fabs(currentLambda) > gramDeterminants[previousRow])
		{
			quotient = round(currentLambda / gramDeterminants[previousRow]);
		}
	}
	
	if(quotient != 0)
	{
		int noOfColumnsToReduce = pivotColumn0 == (noPivotFound) ? inputCols : inputRows + inputCols + 1;
		//Reduce columns
		for(int col = 0; col <= noOfColumnsToReduce; col++)
		{
			augmentedMatrix[currentRow][col] -= quotient * augmentedMatrix[previousRow][col];
		}
		//Update lambda
		gramLambdaSquareMatrix[currentRow][previousRow] -= quotient * gramDeterminants[previousRow];
		for(int col = 0; col <= previousRow -1; col++)
		{
			gramLambdaSquareMatrix[currentRow][col] -= quotient * gramLambdaSquareMatrix[previousRow][col];
		}
	}
	//Return pivot data
	*pivot0 = pivotColumn0;*pivot1 = pivotColumn1;
}

void swap_array_pointers(double *a, double *b)
{
	double tmpd = *a;
	*a = *b;
	*b = tmpd;
}

void SwapRow(int currentRow, int inputRows, int inputCols, double **augmentedMatrix, double **gramLambdaSquareMatrix, double *gramDeterminants)
{
	int previousRow = currentRow - 1;
	for(int col = 0; col <= inputRows + inputCols + 1; col++)
	{
		swap_array_pointers(&augmentedMatrix[currentRow][col], &augmentedMatrix[previousRow][col]);		
	}
	for(int col = 0; col <= previousRow - 1; col++)
	{
		swap_array_pointers(&gramLambdaSquareMatrix[currentRow][col], &gramLambdaSquareMatrix[previousRow][col]);
	}
	
	//Update GramSchmidt coefficients
	double currentLambda = gramLambdaSquareMatrix[currentRow][previousRow];
    	double gramTemp = (gramDeterminants[previousRow - 1] * gramDeterminants[currentRow] + currentLambda * currentLambda) / gramDeterminants[previousRow] ;			
	for(int row = currentRow + 1; row <= inputCols; row++)
	{
		double rowLambda = gramLambdaSquareMatrix[row][currentRow];
		gramLambdaSquareMatrix[row][currentRow] = (gramDeterminants[currentRow] * gramLambdaSquareMatrix[row][previousRow] - currentLambda * rowLambda) / gramDeterminants[previousRow];		
		gramLambdaSquareMatrix[row][previousRow] = (gramTemp * rowLambda + currentLambda * gramLambdaSquareMatrix[row][currentRow]) / gramDeterminants[currentRow];
	}
	gramDeterminants[previousRow] = gramTemp;
}

#define aln   80
#define ald   81
void SolveAugmentedMatrix(int augmentedRows, int augmentedCols, int inputRows, int inputCols, double **augmentedMatrix, double **gramLambdaSquareMatrix, double *gramDeterminants)
{
	int currentRow = 1;int previousRow = 0;int noPivotFound = inputRows + inputCols + 2;
	int pivotColumn0 = 0;int pivotColumn1 = 0;int loopCounter = 0;
	while(currentRow <= inputCols)
	{
		previousRow = currentRow - 1;
		LLLReduce(&pivotColumn0, &pivotColumn1, previousRow, currentRow, augmentedRows, augmentedCols, inputRows, inputCols,augmentedMatrix, gramLambdaSquareMatrix, gramDeterminants);
		int bothRowsHaveNoPivot = (pivotColumn0 == noPivotFound && pivotColumn1 == noPivotFound);
		int lovaszConditionFails = 0;int shouldSwapRows = 0;
		if(bothRowsHaveNoPivot)
		{
			double currentLambda = gramLambdaSquareMatrix[currentRow][previousRow];
    			double gramTemp = gramDeterminants[previousRow - 1] * gramDeterminants[currentRow] + currentLambda * currentLambda;
			bothRowsHaveNoPivot = (gramTemp * ald) < (gramDeterminants[previousRow] * gramDeterminants[previousRow] * aln);
		}
		if(bothRowsHaveNoPivot || (pivotColumn0 <= pivotColumn1 && pivotColumn0 < noPivotFound))
		{
			//Swap and decrease current row
			SwapRow(currentRow, inputRows, inputCols, augmentedMatrix, gramLambdaSquareMatrix, gramDeterminants);
			if(currentRow > 1){currentRow -= 1;}
		}
		else
		{
			//Complete reduction
			for(int i = previousRow - 1; i >= 0; i--)
			{
				LLLReduce(&pivotColumn0, &pivotColumn1, i, currentRow, augmentedRows, augmentedCols, inputRows, inputCols,augmentedMatrix, gramLambdaSquareMatrix, gramDeterminants);
			}currentRow += 1;
		}loopCounter += 1;
	}	
}

void ProcessMatrix(int inputRows, int inputCols, double inputMatrix[inputRows][inputCols], double *target)
{
	int gramDeterminantsLength = inputCols + 3;
	int gramLambdaDimension    = inputCols + 1;
	int augmentedRows = inputCols + 1;
	int augmentedCols = inputRows + inputCols + 2;
	double **augmentedMatrix = Matrix2D(augmentedRows, augmentedCols);
	double **gramLambdaSquareMatrix = Matrix2D(gramLambdaDimension, gramLambdaDimension);
	double *gramDeterminants = NegativeIndexArray(-1, inputCols + 1);
	//Store input matrix in transposed form
	for(int row = 0; row < inputRows; row++)
	{
	for(int col = 0; col < inputCols; col++){augmentedMatrix[col][row + inputCols + 1] = inputMatrix[row][col];}
	}
	//Store desired target
	for(int row = 0; row < inputRows; row++)
	{
		augmentedMatrix[inputCols][inputCols + 1 + row] = target[row];
	}
	//Set identity matrix
	for(int col = 0; col < inputCols+1; col++)
	{
		augmentedMatrix[col][col] = 1.0;
	}
	//Set bottom right entry to 1
	augmentedMatrix[inputCols][inputRows+inputCols+1] = 1.0;
    	//Initialize Gram determinants
	for(int col = -1; col < inputCols+1; col++){gramDeterminants[col] = 1.0;}  	
	PrintAugmentedMatrix(augmentedMatrix, augmentedRows, augmentedCols, inputCols);
	SolveAugmentedMatrix(augmentedRows, augmentedCols, inputRows, inputCols, augmentedMatrix, gramLambdaSquareMatrix, gramDeterminants);
	printf("\nSolved\n");
	PrintAugmentedMatrix(augmentedMatrix, augmentedRows, augmentedCols, inputCols);
	FreeNegativeIndexArray(-1, inputCols + 1, gramDeterminants);
	FreeMatrix2D(augmentedRows, augmentedMatrix);
	FreeMatrix2D(gramLambdaDimension, gramLambdaSquareMatrix);
}

void TestMatrix0()
{
	printf("\n=== Test 1: No integral solution ===\n");
	int inputRows = 2, inputCols = 2;
	double testMatrix[2][2] =
	{
		{2, 0},
		{2, 1}
	};
	double target[2] = {1, 2};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix1()
{
	printf("\n=== Test 2: Indeterminate ===\n");
	int inputRows = 2, inputCols = 3;
	double testMatrix[2][3] =
	{
		{1, 3, 5},
		{4, 6, 8}
	};
	double target[2] = {0, 0}; 
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix2()
{
	printf("\n=== Test 3: Singular square ===\n");
	int inputRows = 3, inputCols = 3;
	double testMatrix[3][3] = 
	{
		{1, 7, 4},
		{2, 8, 5},
		{3, 9, 6}
	};
	double target[3] = {0, 0, 0};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix3()
{
	printf("\n=== Test 4: Overdetermined ===\n");
	int inputRows = 3, inputCols = 2;
	double testMatrix[3][2] = 
	{
		{2, 1},
		{6, 5},
		{7, 6}
	};
	double target[3] = {2, 2, 2};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix4() 
{
	printf("\n=== Test 5: Square ===\n");
	int inputRows = 3, inputCols = 3;
	double testMatrix[3][3] =
	{
		{2, -3, 4},
		{5, 6, 7},
		{8, 9, 10}
	};
	double target[3] = {9, 3, 3};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

// Test Case 6: Modular inverse
void TestMatrix5()
{
	printf("\n=== Test 6: Modular inverse ===\n");
	int inputRows = 1, inputCols = 2;
	double testMatrix[1][2] = 
	{
		{42, -2017}
	};
	double target[1] = {1};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix6()
{
	printf("\n=== Test 7: Indian kuttaka problem ===\n");
	int inputRows = 1, inputCols = 2;
	double testMatrix[1][2] =
	{
		{195, -221}
	};
	double target[1] = {65};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix7()
{
	printf("\n=== Test 8: Bachet banquet problem ===\n");
	int inputRows = 2, inputCols = 3;
	double testMatrix[2][3] =
	{
		{1, 1, 1},
		{12, 9, 1}
	};
	double target[2] = {41, 120};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix8()
{
	printf("\n=== Test 9: Malm ===\n");
	int inputRows = 2, inputCols = 4;
	double testMatrix[2][4] =
	{
		{1, 1, 1, 1},
		{165, 235, 85, 389}
	};
	double target[2] = {80, 16324};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix9() 
{
	printf("\n=== Test 10: Sunzi Suanjing ===\n");
	int inputRows = 3, inputCols = 4;
	double testMatrix[3][4] =
	{
		{1, 3, 0, 0},
		{1, 0, 5, 0},
		{1, 0, 0, 7}
	};
	double target[3] = {2, 3, 2};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix10()
{
	printf("\n=== Test 11: Sunzi-type 2 ===\n");
	int inputRows = 2, inputCols = 3;
	double testMatrix[2][3] = 
	{
		{17, 7, 0},
		{11, 0, 15}
	};
	double target[2] = {-1, -2};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix11() 
{
	printf("\n=== Test 12: Shushu jiuzhang ===\n");
	int inputRows = 8, inputCols = 9;
	double testMatrix[8][9] = 
	{
		{1, 130, 0, 0, 0, 0, 0, 0, 0},
		{1, 0, 110, 0, 0, 0, 0, 0, 0},
		{1, 0, 0, 120, 0, 0, 0, 0, 0},
		{1, 0, 0, 0, 60, 0, 0, 0, 0},
		{1, 0, 0, 0, 0, 25, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 100, 0, 0},
		{1, 0, 0, 0, 0, 0, 0, 50, 0},
		{1, 0, 0, 0, 0, 0, 0, 0, 20}
	};
	double target[8] = {-60, -30, -10, -10, 10, 10, 10, 10};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix12() 
{
	printf("\n=== Test 13: 5 sailor coconut puzzle ===\n");
	int inputRows = 6, inputCols = 7;
	double testMatrix[6][7] =
	{
		{1, -5, 0, 0, 0, 0, 0},
		{0, 4, -5, 0, 0, 0, 0},
		{0, 0, 4, -5, 0, 0, 0},
		{0, 0, 0, 4, -5, 0, 0},
		{0, 0, 0, 0, 4, -5, 0},
		{0, 0, 0, 0, 0, 4, -5}
	};
	double target[6] = {1, 1, 1, 1, 1, 0};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix13() 
{
	printf("\n=== Test 14: Unbounded knapsack ===\n");
	int inputRows = 3, inputCols = 6;
	double testMatrix[3][6] = 
	{
		{3000, 1800, 2500, 1, 0, 0},
		{3, 2, 20, 0, 1, 0},
		{25, 15, 2, 0, 0, 1}
	};
	double target[3] = {54500, 250, 250};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix14() 
{
	printf("\n=== Test 15: Subset sum ===\n");
	int inputRows = 1, inputCols = 9;
	double testMatrix[1][9] = 
	{
		{575, 436, 1586, 1030, 1921, 569, 721, 1183, 1570}
	};
	double target[1] = {6665};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix15()
{
	printf("\n=== Test 16: HMM extended gcd ===\n");
	int inputRows = 1, inputCols = 10;
	double testMatrix[1][10] = 
	{
		{763836, 1066557, 113192, 1785102, 1470060, 3077752, 114793, 3126753, 1997137, 2603018}
	};
	double target[1] = {1};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix16() 
{
	printf("\n=== Test 17: Fibonacci segment ===\n");
	int inputRows = 1, inputCols = 8;
	double testMatrix[1][8] = 
	{
		{13, 21, 34, 55, 89, 144, 233, 377}
	};
	double target[1] = {1};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


void TestMatrix17()
{
	printf("\n=== Test 18: Fibonacci to Lucas ===\n");
	int inputRows = 8, inputCols = 8;
	double testMatrix[8][8] = 
	{
		{1, 0, 0, 0, 0, 0, 18, -7},
		{1, 1, 0, 0, 0, 0, -11, 4},
		{-1, 1, 1, 0, 0, 0, 7, -3},
		{0, -1, 1, 1, 0, 0, -4, 1},
		{0, 0, -1, 1, 1, 0, 3, -1},
		{0, 0, 0, -1, 1, 1, -1, 1},
		{0, 0, 0, 0, -1, 1, 1, 0},
		{0, 0, 0, 0, 0, -1, -1, 0}
	};
	double target[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix18()
{
	printf("\n=== Test 19: Hnf example ===\n");
	int inputRows = 10, inputCols = 10;
	double testMatrix[10][10] = 
	{
		{3, 11, 31, 69, 131, 223, 351, 521, 739, 1011},
		{7, 36, 113, 262, 507, 872, 1381, 2058, 2927, 4012},
		{13, 77, 249, 583, 1133, 1953, 3097, 4619, 6573, 9013},
		{21, 134, 439, 1032, 2009, 3466, 5499, 8204, 11677, 16014},
		{31, 207, 683, 1609, 3135, 5411, 8587, 12813, 18239, 25015},
		{43, 296, 981, 2314, 4511, 7788, 12361, 18446, 26259, 36016},
		{57, 401, 1333, 3147, 6137, 10597, 16821, 25103, 35737, 49017},
		{73, 522, 1739, 4108, 8013, 13838, 21967, 32784, 46673, 64018},
		{91, 659, 2199, 5197, 10139, 17511, 27799, 41489, 59067, 81019},
		{111, 812, 2713, 6414, 12515, 21616, 34317, 51218, 72919, 100020}
	};
	double target[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix19()
{
	printf("\n=== Test 20: Gauss arctan (simplified) ===\n");
	int inputRows = 1, inputCols = 4;
	double testMatrix[1][4] = 
	{
		{-1, 41841, 175421, 554985} 
	};
	double target[1] = {7853982}; 
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}

void TestMatrix20()
{
	printf("\n=== Example 7.1: Extended GCD ===\n");
	int inputRows = 1, inputCols = 4;
	double testMatrix[1][4] = 
	{
		{116085838, 181081878, 314252913, 10346840} 
	};
	double target[1] = {0}; 
	ProcessMatrix(inputRows, inputCols, testMatrix, target);
}


int main()
{
	/*TestMatrix0();
	TestMatrix1();
	TestMatrix2();
	TestMatrix3();
	TestMatrix4();
	TestMatrix5();
	TestMatrix6();
	TestMatrix7();
	TestMatrix8();
	TestMatrix9();
	TestMatrix10();
	TestMatrix11();
	TestMatrix12();
	TestMatrix13();
	TestMatrix14();
	TestMatrix15();
	TestMatrix16();
	TestMatrix17();
	TestMatrix18();
	TestMatrix19();*/
	TestMatrix20();
	return 0;
}
