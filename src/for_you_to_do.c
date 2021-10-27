#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, j, k;
    double *temprow = (double*) malloc(sizeof(double) * n);
    
    for (i = 0; i < (n - 1); i ++)
    {
        int jp = i;
        pivot = fabs(A[i * n + i]);
        
        // find pivoting
        for (j = i + 1; j < n; j ++)
        {
            if (fabs(A[j * n + i]) > pivot){
                pivot = fabs(A[j * n + i]);
                jp = j;
            }
        }
        
        // swap
        if (pivot == 0)
        {
            printf("Warning: Matrix A is singular.")
            return -1;
        }
        
        else
        {
            if (jp != i)
            {
                int temp = ipiv[i];
                ipiv[i] = ipiv[jp];
                ipiv[jp] = temp;
                
                memcpy(temprow, A + i * n, n * sizeof(double));
                memcpy(A + i * n, A + jp * n, n * sizeof(double));
                memcpy(A + jp * n, temprow, n * sizeof(double));
            }    
        }
        
        // factorization
        for (j = i + 1; j < n; j ++)
            A[j * n + i] = A[j * n + i] / A[i * n + i];
        for (j = i + 1; j < n; j ++)
        {
            for (k = i + 1; k < n; k ++)
            {
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
            }
        }
    }
    free(temprow);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    double *y = (double*) malloc(n * sizeof(double));
    int i, j;
    double sum;
    
    switch(UPLO)
    {
        case 'L':
            y[0] = B[ipiv[0]];
            for (i = 1; i < n; i ++)
            {
                sum = 0.0;
                for (j = 0; j < i; j ++)
                {
                    sum += y[j] * A[i * n + j];
                }
                y[i] = B[ipiv[i]] - sum;
            }
            break;
        case 'U':
            y[n - 1] = B[n - 1] / A[(n - 1) * n + n - 1];
            for (i = n - 2; i >= 0; i --)
            {
                for (j = n - 1; j > i; j --)
                {
                    sum += y[j] * A[i * n + j];
                }
                y[i] = (B[i] - sum) / A[i * n + i];
            }
            break;
        default:
            printf("\nCheck your UPLO input, it must be L or U.\n");
            break;
    }
    
    memcpy(B, y, sizeof(double) * n);
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib;
    
    for (ib = 0; ib < n; ib += b)
    {
        end = ib + b - 1;
        
        int i, j, k;
        double *temprow = (double*) malloc(sizeof(double) * b);

        for (i = ib; i < (n - 1); i ++)
        {
            int jp = i;
            pivot = fabs(A[i * n + i]);

            // find pivoting
            for (j = i + 1; j < n; j ++)
            {
                if (fabs(A[j * n + i]) > pivot){
                    pivot = fabs(A[j * n + i]);
                    jp = j;
                }
            }

            // swap
            if (pivot == 0)
            {
                printf("Warning: Matrix A is singular.")
                return -1;
            }

            else
            {
                if (jp != i)
                {
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[jp];
                    ipiv[jp] = temp;
                    
                    memcpy(temprow, A + i * n, n * sizeof(double));
                    memcpy(A + i * n, A + jp * n, n * sizeof(double));
                    memcpy(A + jp * n, temprow, n * sizeof(double));
                }    
            }

            // factorization
            for (j = i + 1; j < n; j ++)
                A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (j = i + 1; j < n; j ++)
            {
                for (k = i + 1; k < ib + b; k ++)
                {
                    A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
                }
            }
        }
        free(temprow);
        
        // built LL ^ (-1)
        double *ll = (double*) calloc(sizeof(double) (b * b));
        
        for (int i = 0; i < b; i ++)
        {
            ll[i * b + i] = 0;
            for (int j = i + 1; j < b; j ++)
            {
                for (k = i; k < j; k ++)
                {
                    ll[j * b + i] = ll[j * b + i] - A[ib * n + ib + j * n + k] * ll[k * b + i];
                }
            }
        }
        
        // update A(ib:end, end+1:n)
        for (i = ib; i <= end; i ++)
        {
            for (j = end + 1; i < n
        }
        
    }
    
    
    return 0;
}

