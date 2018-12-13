#include "defs.h"

using namespace std;

numvector<double, dimPh> inverseRotate(const numvector<double, dimPh>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] - n.y() * sol[2],  n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}
numvector<double, dimPh> rotate(const numvector<double, dimPh>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] + n.y() * sol[2], - n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}

///
/// Overloading of vector<numvector<double,dimPhys * dimShapes>>.
///

vector<numvector<double, dimS>> operator * (const vector<numvector<double, dimS>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = dimS;

	vector<numvector<double, dimS>> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			m[cell][val] *= b;

	return m;
};
vector<numvector<double, dimS>> operator + (const vector<numvector<double, dimS>>& b, const vector<numvector<double, dimS>>& a)
{
	size_t dimx = a.size();
	size_t dimy = dimS;

	vector<numvector<double, dimS>> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			m[cell][val] = a[cell][val] + b[cell][val];

	return m;
};
void sum(const vector<numvector<double, dimS>>& a, const vector<numvector<double, dimS>>& b, vector<numvector<double, dimS>>& res)
{
	size_t dimx = a.size();
	size_t dimy = dimS;

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			res[cell][val] = a[cell][val] + b[cell][val];
};


///
/// Overloading for vector<vector<double>>.
///

vector<vector<double>>& operator += (vector<vector<double>>& a, const vector<vector<double>>& b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			a[cell][val] += b[cell][val];
	return a;
};
vector<vector<double>>& operator -= (vector<vector<double>>& a, const vector<vector<double>>& b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			a[cell][val] -= b[cell][val];
	return a;
};
vector<vector<double>>& operator *= (vector<vector<double>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			a[cell][val] *= b;
	return a;
};
vector<vector<double>> operator * (const vector<vector<double>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

	vector<vector<double>> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			m[cell][val] *= b;
	return m;
};
vector<vector<double>> operator * (const double b, const vector<vector<double>>& a)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

	vector<vector<double>> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
#pragma omp simd
		for (size_t val = 0; val < dimy; ++val)
			m[cell][val] *= b;
	return m;
};

///
/// Overloading for vector<double>.
///

vector<double>& operator *= (vector<double>& a, const double b)
{
	size_t dimx = a.size();
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		a[cell] *= b;
	return a;
};
vector<double>& operator += (vector<double>& a, const vector<double>& b)
{
	size_t dimx = a.size();
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		a[cell] += b[cell];
	return a;
};
vector<double>& operator -= (vector<double>& a, const vector<double>& b)
{
	size_t dimx = a.size();
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		a[cell] -= b[cell];
	return a;
};
vector<double> operator + (const vector<double>& a, const vector<double>& b)
{
	size_t dimx = a.size();
	vector<double> res; res.resize(dimx);
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		res[cell] = a[cell] + b[cell];
	return res;
};
vector<double> operator - (const vector<double>& a, const vector<double>& b)
{
	size_t dimx = a.size();
	vector<double> res; res.resize(dimx);
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		res[cell] = a[cell] - b[cell];
	return res;
};
vector<double> operator * (const vector<double>& a, const double b)
{
	size_t dimx = a.size();

	vector<double> m(a);
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		m[cell] *= b;
	return m;
};
vector<double> operator / (const vector<double>& a, const double b)
{
	size_t dimx = a.size();

	vector<double> m(a);
#pragma omp simd
	for (size_t cell = 0; cell < dimx; ++cell)
		m[cell] /= b;
	return m;
};
