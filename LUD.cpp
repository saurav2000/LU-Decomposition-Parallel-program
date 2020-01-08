#include <iostream>
#include <vector>
#include <chrono>
using namespace std;

#ifndef vd
#define vd vector<double>
#endif
#ifndef vi
#define vi vector<int>
#endif

#define N 8001

int n, t;

vi pi(N);
vector<vd> a(N, vd(N));
vector<vd> u(N, vd(N));
vector<vd> l(N, vd(N));

int LUD()
{
	++n;
	
	--n;

	for(int i=1;i<=n;++i)
	{
		pi[i] = i;
		l[i][i] = 1;
	}


	for(int k=1;k<=n;++k)
	{
		double max = 0;
		int k_ = 0;
		for (int i=k; i<=n; ++i)
		{
			if(max<fabs(a[i][k]))
			{
				max = fabs(a[i][k]);
				k_ = i;
			}
		}

		if(!max)
			return 1;

		swap(pi[k], pi[k_]);
		a[k].swap(a[k_]);

		for(int i=1;i<=k-1;++i)
			swap(l[k][i], l[k_][i]);

		u[k][k] = a[k][k];

		for(int i=k+1; i<=n; ++i)
		{
			l[i][k] = a[i][k] / u[k][k];
			u[k][i] = a[k][i];
		}

		for(int i=k+1; i<=n; ++i)
		{
			for(int j=k+1; j<=n; ++j)
				a[i][j] = a[i][j] - (l[i][k] * u[k][j]);
		}
	}


	return 0;
}

double verify()
{
	vector<vd> res(1, vd(n+1));

	for(int i=1;i<=n;++i)
	{
		for(int j=1;j<=n;++j)
			res.push_back(a[pi[i]]);
	}

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
	cin >> n >> t;
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	LUD();
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
	
	cout<<"L21 Norm: "<<verify()<<"\n";
	return 0;
}