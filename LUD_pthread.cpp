#include <iostream>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <pthread.h>

using namespace std;

//Maximum number of threads
#define T 32

// struct for the thread argument
struct arg_t
{
	int low, high, id;
};

// k is a global variable for all threads
int n, t, k;

//Pointers for all the data required
int *pi;
double **a, **orig, **u, **l, **res;

//Thread array for pthreads
pthread_t threads[T];
double max_arr[T];
int k_arr[T];

/**
* This function initialises all the matrices
* and pi appropriately using drand48
*/
void initialise()
{
	// Memory Allocation of all pointers
	pi = (int*)malloc(sizeof(int) * (n+1));
	a = (double**)malloc(sizeof(double*) * (n+1));
	orig = (double**)malloc(sizeof(double*) * (n+1));
	l = (double**)malloc(sizeof(double*) * (n+1));
	u = (double**)malloc(sizeof(double*) * (n+1));
	res = (double**)malloc(sizeof(double*) * (n+1));
	for(int i=0;i<n+1;++i)
	{
		a[i] = (double*)(malloc(sizeof(double) * (n+1)));
		orig[i] = (double*)(malloc(sizeof(double) * (n+1)));
		l[i] = (double*)(malloc(sizeof(double) * (n+1)));
		u[i] = (double*)(malloc(sizeof(double) * (n+1)));
		res[i] = (double*)(malloc(sizeof(double) * (n+1)));
	}

	// Initialising a, l and u
	srand48(time(NULL));
	for(int i=1;i<=n;++i)
	{
		for(int j=1;j<=n;++j)
		{
			double x = drand48();
			orig[i][j] = x;
			a[i][j] = x;
			l[i][j] = 0;
			u[i][j] = 0;
		}

		pi[i] = i;
		l[i][i] = 1;
	}
}

/**
* This loop is the thread routine which runs the 
* step where the matrix a is decomposed using l and u
*/
void* n2_loop(void *argv)
{
	struct arg_t *arg = ((struct arg_t*)argv);
	int low = arg->low, high=arg->high;

	for(int i=low;i<=high;++i)
	{
		for(int j=k+1;j<=n;++j)
			a[i][j] -= (l[i][k] * u[k][j]);
	}
	free(arg);
}

/**
* This is the main function where the entire 
* steps of LU decomposition happens
*/
int LUD()
{
	for(k=1; k<=n; ++k)
	{
		// Computing maximum
		double max = 0;
		int k_ = 0;
		for(int i=k; i<=n; ++i)
		{
			if(max<fabs(a[i][k]))
			{
				max = fabs(a[i][k]);
				k_ = i;
			}
		}

		if(!k_)
			return 1;

		// Swapping the appropriate rows of pi, a and l
		swap(pi[k], pi[k_]);
		swap(a[k], a[k_]);

		for(int i=1;i<=k-1;++i)
			swap(l[k][i], l[k_][i]);

		u[k][k] = a[k][k];

		for(int i=k+1; i<=n; ++i)
		{
			l[i][k] = a[i][k] / u[k][k];
			u[k][i] = a[k][i];
		}
		
		// x is the size of each thread
		// n^2 loop
		int x = (n-k)/t;
		for(int i=0;i<t;++i)
		{
			// Creating the struct with lower and higher 
			// bound of the thread
			struct arg_t *arg = (struct arg_t*)malloc(sizeof(struct arg_t));
			arg->low = (i*x) + (k+1);
			if(i!=t-1)
				arg->high = (i+1)*x +k;
			else
				arg->high = n;
			arg->id = i;
			// Creating the thread
			pthread_create(&threads[i], NULL, n2_loop, (void*)arg);
			
		}
		// Joining all threads
		for(int i=0;i<t;++i)
			pthread_join(threads[i], NULL);
	}

	return 0;
}

/**
* This function calculates the L21 norm
* of the matrix and returns it
*/
double verify()
{
	// Res = PA
	for(int i=1;i<=n;++i)
		for(int j=1;j<=n;++j)
				res[i][j] = orig[pi[i]][j];

	// Res-=LU
	for(int i=1; i<=n; ++i)
		for(int j=1; j<=n; ++j)
			for(int k=1; k<=n; ++k)
				res[i][j]-= l[i][k] * u[k][j];

	// Calculating Norm
	double l21 = 0.0;
	for(int i=1; i<=n; ++i)
	{
		double term = 0;
		for(int j=1; j<=n; ++j)
			term+= (res[i][j] * res[i][j]);
		l21+= sqrt(term);
	}

	return l21;
}

int main(int argc, char const *argv[])
{
	n = atoi(argv[1]);
	t = atoi(argv[2]);
	initialise();

	auto t1 = chrono::high_resolution_clock::now();
	int err = LUD();
	if(err)
	{
		cout<<" Singular Matrix\n";
		return 0;
	}
	//Printing time
	auto t2 = chrono::high_resolution_clock::now();
	auto count = std::chrono::duration_cast<std::chrono::duration<double> >(t2-t1).count();
	cout<<count<<"\n";
	// cout<<"L21 Norm: "<<verify()<<"\n";
	return 0;
}