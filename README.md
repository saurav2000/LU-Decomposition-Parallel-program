# LU-Decomposition-Parallel-program
## About LU Decomposition
**Lower–Upper** (**LU**) **decomposition** factors a matrix as the product of a lower triangular matrix and an upper triangular matrix. The product sometimes includes a permutation matrix as well. LU decomposition can be viewed as the matrix form of [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination "Gaussian elimination"). Source: [Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)

## Running the code
The code can be compiled using makefile. It has 3 modes
- Sequential
- Pthreads
- OpenMP

The code is compiled by
> For sequential code:
> `make seq`
> 
> For code parallelised by Pthreads:
> `make pthread`
> 
> For code parallelised by OpenMP:
> `make omp`

The executable generated is a.out and takes 2 parameters for Pthreads and OpenMP and only 1 parameter for the sequential code.
> For Pthread and OpenMP:
> `./a.out [SIZE_OF_MATRIX] [NO_OF_THREADS]`
> 
> For sequential:
> `./a.out [SIZE_OF_MATRIX]`
