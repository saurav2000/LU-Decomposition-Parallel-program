#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <pthread.h>

using namespace std;

# define T 32

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

void initialise()
{
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

void* max_loop(void *argv)
{
	struct arg_t *arg = ((struct arg_t*)argv);
	int l = arg->low, h=arg->high, id = arg->id;
	double max = 0;
	int k_ = 0;
	for(int i=l;i<=h;++i)
	{
		if(max<fabs(a[i][k]))
		{
			max = fabs(a[i][k]);
			k_ = i;
		}
	}
	max_arr[id] = max;
	k_arr[id] = k_;

	free(arg);
}

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

void* swap_loop(void *argv)
{
	struct arg_t *arg = ((struct arg_t*)argv);
	int low = arg->low, high=arg->high, k_=arg->id;

	for(int i=low;i<=high;++i)
		swap(l[k][i], l[k_][i]);
	free(arg);	
}

void* div_loop(void *argv)
{
	struct arg_t *arg = ((struct arg_t*)argv);
	int low = arg->low, high=arg->high;

	for(int i=low;i<=high;++i)
	{
		l[i][k] = a[i][k] / u[k][k];
		u[k][i] = a[k][i];
	}
	free(arg);
}

int LUD()
{
	for(k=1; k<=n; ++k)
	{
		// Loop to find out maximum
		int x = (n-k+1)/t;
		for(int i=0;i<t;++i)
		{
			struct arg_t *arg = (struct arg_t*)malloc(sizeof(struct arg_t));
			arg->low = i*x + k;
			if(i!=t-1)
				arg->high = (i+1)*x + k - 1;	
			else
				arg->high = n;
			arg->id = i;
			pthread_create(&threads[i], NULL, max_loop, (void*)arg);
		}
		for(int i=0;i<t;++i)
			pthread_join(threads[i], NULL);
		
		double max = 0;
		int k_ = 0;
		for(int i=0;i<t;++i)
		{
			if(max<max_arr[i])
			{
				max = max_arr[i];
				k_ = k_arr[i];
			}
		}

		if(!k_)
			return 1;
		//Swapping pi and a:

		// double max = 0;
		// int k_ = 0;
		// for(int i=k; i<=n; ++i)
		// {
		// 	if(max<fabs(a[i][k]))
		// 	{
		// 		max = fabs(a[i][k]);
		// 		k_ = i;
		// 	}
		// }

		// if(!k_)
		// 	return 1;

		swap(pi[k], pi[k_]);
		swap(a[k], a[k_]);

		x = (k-1)/t;
		for(int i=0;i<t;++i)
		{
			struct arg_t *arg = (struct arg_t*)malloc(sizeof(struct arg_t));
			arg->low = i*x + 1;
			if(i!=t-1)
				arg->high = (i+1)*x;	
			else
				arg->high = k-1;
			arg->id = k_;
			pthread_create(&threads[i], NULL, swap_loop, (void*)arg);
		}

		for(int i=0;i<t;++i)
			pthread_join(threads[i], NULL);

		u[k][k] = a[k][k];

		x = (n-k)/t;
		for (int i=0;i<t;++i)
		{
			struct arg_t *arg = (struct arg_t*)malloc(sizeof(struct arg_t));
			arg->low = (i*x) + (k+1);
			if(i!=t-1)
				arg->high = (i+1)*x +k;
			else
				arg->high = n;
			pthread_create(&threads[i], NULL, div_loop, (void*)arg);
		}
		
		for(int i=0;i<t;++i)
			pthread_join(threads[i], NULL);
		
		// n^2 loop
		
		for(int i=0;i<t;++i)
		{
			struct arg_t *arg = (struct arg_t*)malloc(sizeof(struct arg_t));
			arg->low = (i*x) + (k+1);
			if(i!=t-1)
				arg->high = (i+1)*x +k;
			else
				arg->high = n;
			arg->id = i;
			pthread_create(&threads[i], NULL, n2_loop, (void*)arg);
			
		}
		for(int i=0;i<t;++i)
			pthread_join(threads[i], NULL);
	}

	return 0;
}

double verify()
{
	for(int i=1;i<=n;++i)
		for(int j=1;j<=n;++j)
				res[i][j] = orig[pi[i]][j];

	for(int i=1; i<=n; ++i)
		for(int j=1; j<=n; ++j)
			for(int k=1; k<=n; ++k)
				res[i][j]-= l[i][k] * u[k][j];

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
	auto t2 = chrono::high_resolution_clock::now();
	auto count = std::chrono::duration_cast<std::chrono::duration<double> >(t2-t1).count();
	cout<<count<<"\n";
	// cout<<"L21 Norm: "<<verify()<<"\n";
	return 0;
}