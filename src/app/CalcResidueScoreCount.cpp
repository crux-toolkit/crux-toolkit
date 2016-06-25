/* calcResidueScoreCount.cpp
 *
 * Calculates counts of peptides with various residue evidence scores, given a preprocessed
 * residue evidence matrix, using dynamic programming.
 *
 * This version:
 *  - calculates most of dimension and indexing variables required for dynamic
 *      programming inside function, instead of externally in MATLAB
 * - uses uniform amino acid probabilities for all positions in peptide
 *  
 * This is a MEX-file to be called from MATLAB.
 *
 * Written by Jeff Howbert, August, 2015.
 */

/* Added by Andy Lin March 2016
 * Edited to work within Crux code instead of with orginal Matlab code
 */
         
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "CalcResidueScoreCount.h"

void calcResidueScoreCount(
  int nAa,
  int pepMassInt,
  vector<vector<double> >& residueEvidenceMatrix,
  vector<int>& aaMass, 
  double* probN, //not being used at present
  double* probI, //not being used at present
  double* probC, //not being used at present
  int NTermMass, //this is NTermMassBin
  int CTermMass, //this is CTermMassBin
  int minAaMass,
  int maxAaMass,
  int maxEvidence,
  int maxScore
)
//void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    int minEvidence     = 0;
    int minScore        = 0;
    
    // copy residue evidence matrix to C++ memory
    int** residEvid = new int* [ nAa ];
    for ( int row = 0; row < nAa; row++ ) {
        residEvid[ row ] = new int[ pepMassInt ];
        for ( int col = 0; col < pepMassInt; col++ ) {
            residEvid[row][col] = residueEvidenceMatrix[row][col];
        }
    }
    
    // internal variables
    int row;
    int col;
    int ma;
    int evid;
    int de;
    int evidRow;
    double sumScore;

    int bottomRowBuffer = maxEvidence;
    int topRowBuffer = -minEvidence;
    int colBuffer = maxAaMass;
    int colStart = NTermMass;
    int nRow = bottomRowBuffer - minScore + 1 + maxScore + topRowBuffer;
    int nCol = colBuffer + pepMassInt;
    int rowFirst = bottomRowBuffer + 1;
    int rowLast = rowFirst - minScore + maxScore;
    int colFirst = colStart + 1;
    int colLast = pepMassInt - CTermMass;
    int initCountRow = bottomRowBuffer - minScore + 1;
    int initCountCol = maxAaMass + colStart;
    // convert to zero-based indexing
    rowFirst = rowFirst - 1;
    rowLast = rowLast - 1;
    colFirst = colFirst - 1;
    colLast = colLast - 1;
    initCountRow = initCountRow - 1;
    initCountCol = initCountCol - 1;
    
    double** dynProgArray = 0;
    dynProgArray = new double* [ nRow ];
    for ( row = 0; row < nRow; row++ ) {
        dynProgArray[ row ] = new double[ nCol ];
        for ( col = 0; col < nCol; col++ ) {
            dynProgArray[ row ][ col ] = 0.0;
        }
    }
    dynProgArray[ initCountRow ][ initCountCol ] = 1.0;    // initial count of peptides with mass = NTermMass
    int* aaMassCol = new int[ nAa ];
    // populate matrix with scores for first (i.e. N-terminal) amino acid in sequence
    for ( de = 0; de < nAa; de++ ) {
        ma = aaMass[ de ];
//         row = initCountRow + evidence[ ma + NTermMass - 1 ];    //&& -1 is to account for zero-based indexing in evidence vector
        row = initCountRow + residEvid[ de ][ ma + NTermMass - 1 ];    //&& -1 is to account for zero-based indexing in evidence vector
        col = initCountCol + ma;
        if ( col <= maxAaMass + colLast ) {
//             dynProgArray[ row ][ col ] += dynProgArray[ initCountRow ][ initCountCol ] * probN[ de ];
            dynProgArray[ row ][ col ] += dynProgArray[ initCountRow ][ initCountCol ]; //&& ignore aa probs for now
        }
    }
    dynProgArray[ initCountRow ][ initCountCol ] = 0.0;    // set to zero now that score counts for first amino acid are in matrix
    // populate matrix with score counts for non-terminal amino acids in sequence 
    for ( ma = colFirst; ma < colLast; ma++ ) {
        col = maxAaMass + ma;
//         evid = evidence[ ma ];
        for ( de = 0; de < nAa; de++ ) {
            aaMassCol[ de ] = col - aaMass[ de ];
        }
        for ( row = rowFirst; row <= rowLast; row++ ) {
//             evidRow = row - evid;
            sumScore = dynProgArray[ row ][ col ];
            for ( de = 0; de < nAa; de++ ) {
                evidRow = row - residEvid[ de ][ ma ];
//                 sumScore += dynProgArray[ evidRow ][ aaMassCol[ de ] ] * probI[ de ];
                sumScore += dynProgArray[ evidRow ][ aaMassCol[ de ] ]; //&& ignore aa probs for now
            }
            dynProgArray[ row ][ col ] = sumScore;
        }
    }
    // populate matrix with score counts for last (i.e. C-terminal) amino acid in sequence
    ma = colLast;
    col = maxAaMass + ma;
//     evid = evidence[ ma ];    //&& wrong
    evid = 0;                       // no evidence should be added for last amino acid in sequence
    for ( de = 0; de < nAa; de++ ) {
        aaMassCol[ de ] = col - aaMass[ de ];
    }
    for ( row = rowFirst; row <= rowLast; row++ ) {
        evidRow = row - evid;
        sumScore = 0.0;
        for ( de = 0; de < nAa; de++ ) {
//             sumScore += dynProgArray[ evidRow ][ aaMassCol[ de ] ] * probC[ de ];  // C-terminal residue
            sumScore += dynProgArray[ evidRow ][ aaMassCol[ de ] ];  //&& ignore aa probs for now
        }
        dynProgArray[ row ][ col ] = sumScore;
    }

  	plhs[ 0 ] = mxCreateNumericMatrix( nRow, 1, mxDOUBLE_CLASS, mxREAL );
   	plhs[ 1 ] = mxCreateNumericMatrix( 1, 1, mxINT32_CLASS, mxREAL );
    double* scoreCount = mxGetPr( plhs[ 0 ] );
    int* scoreOffset = ( int* )mxGetPr( plhs[ 1 ] );
    
    int colScoreCount = maxAaMass + colLast;
    for ( int row = 0; row < nRow; row++ ) {
        scoreCount[ row ] = dynProgArray[ row ][ colScoreCount ];
    }
    *scoreOffset = initCountRow + 1;     // back to one-based indexing
    
    // clean up
    for( int row = 0; row < nRow; row++ ) {
        delete [] dynProgArray[ row ];
    }
    delete [] dynProgArray;
    delete [] aaMassCol;
    for ( int row = 0; row < nAa; row++ ) {
        delete [] residEvid[ row ];
    }
    delete [] residEvid;
}
