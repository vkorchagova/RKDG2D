#include "defs.h"

using namespace std;

numvector<double, 5> inverseRotate(const numvector<double, 5>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] - n.y() * sol[2],  n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}

numvector<double, 5> rotate(const numvector<double, 5>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] + n.y() * sol[2], - n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}


vector<numvector<double, dim>> operator * (const vector<numvector<double, dim>>& a, const double b)
{
    size_t dimx = a.size();
    size_t dimy = dim;

    vector<numvector<double, dim>> m(a);

    for (size_t cell = 0; cell < dimx; ++cell)
    for (size_t val = 0; val < dimy; ++val)
        m[cell][val] *= b;

    return m;
};

vector<numvector<double, dim>> operator + (const vector<numvector<double, dim>>& b, const vector<numvector<double, dim>>& a)
{
    size_t dimx = a.size();
    size_t dimy = dim;

    vector<numvector<double, dim>> m(a);

    for (size_t cell = 0; cell < dimx; ++cell)
    for (size_t val = 0; val < dimy; ++val)
        m[cell][val] = a[cell][val] + b[cell][val];

    return m;
};


void sum (const vector<numvector<double, dim>>& a, const vector<numvector<double, dim>>& b, vector<numvector<double, dim>>& res)
{
    size_t dimx = a.size();
    size_t dimy = dim;

    for (size_t cell = 0; cell < dimx; ++cell)
    for (size_t val = 0; val < dimy; ++val)
        res[cell][val] = a[cell][val] + b[cell][val];
};

// return LU-factorization of matrix (+ transform of rhs)

vector<vector<double> > forwardGauss(const vector<vector<double> >& data, const bool partChoice)
{
    int n = data.size();

    vector<vector<double> > LU = data;

    for (int i = 0; i < n; ++i) //for all rows
    {
        if (partChoice)
        {
            int maxI = maxAbsPosition(LU,i);
            ChangeRows(LU,i,maxI);
        }

        if (fabs(LU[i][i]) < 1e-12)
        {
            cout << "Bad matrix: element " << i << " in diag = " << LU[i][i] << endl;
            exit(0);
        }

        for (int k = i+1; k < n; ++k) //just see elements lower than matrix diag
        {
            double c = LU[k][i]/LU[i][i];

            for (int j = i; j < n; ++j)
            {
                LU[k][j] = LU[k][j] - c*LU[i][j];
            }

            //for right part
            LU[k][n] = LU[k][n] - c*LU[i][n];
        }
    }

    return LU;
}

vector<double > reverseGauss(const vector<vector<double> >& data)
{
    int n = data.size();

    vector<double> solution(n,0);

    for (int i = n-1; i >=0 ; --i)
    {
        solution[i] = data[i][n] / data[i][i];

        for (int j = n-1; j > i; --j)
        {
            solution[i] -= solution[j]*data[i][j]/data[i][i];
        }
    }

    return solution;
}

int maxAbsPosition(const vector<vector<double > >& data, const int i)
{
    int n = data.size();

    int maxI = i;
    double maxAbs = 0;

    for (int k = i+1; k < n; ++k)
    {
        if (maxAbs < fabs(data[k][i]))
        {
            maxAbs = fabs(data[k][i]);
            maxI = k;
        }
    }

    return maxI;
}

void ChangeRows(vector<vector<double > >& data, const int p, const int q)
{
    if (p == q)
    {
        return;
    }

    int n = data[p].size();

    for (int j = 0; j < n; ++j)
    {
        data[p][j] -= data[q][j];
        data[q][j] += data[p][j];
        data[p][j] = data[q][j] - data[p][j];
    }
}


//-----------------------------------------------------------------------------------------

//Прибавление к одному тройному массиву другого тройного массива
vector<vector<vector<double>>>& operator += (vector<vector<vector<double>>>& a, const vector<vector<vector<double>>>& b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();
	size_t dimz = a[0][0].size();


	for (size_t cell = 0; cell < dimx; ++cell)
	for (size_t shape = 0; shape < dimy; ++shape)
	for (size_t val = 0; val < dimz; ++val)
	
		a[cell][shape][val] += b[cell][shape][val];
	return a;
};

//Прибавление к одному двумерному массиву другого двумерного массива
vector<vector<double>>& operator += (vector<vector<double>>& a, const vector<vector<double>>& b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();


    for (size_t cell = 0; cell < dimx; ++cell)
    for (size_t val = 0; val < dimy; ++val)
        a[cell][val] += b[cell][val];
    return a;
};

//Вычитание из одного двумерного массива другого двумерного массива
vector<vector<double>>& operator -= (vector<vector<double>>& a, const vector<vector<double>>& b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();


	for (size_t cell = 0; cell < dimx; ++cell)
	for (size_t val = 0; val < dimy; ++val)
		a[cell][val] -= b[cell][val];
	return a;
};

//Домножение двумерного массива на число
vector<vector<double>>& operator *= (vector<vector<double>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

    for (size_t cell = 0; cell < dimx; ++cell)
    for (size_t val = 0; val < dimy; ++val)
        a[cell][val] *= b;
    return a;
};

//Домножение трехмерного массива на число
vector<vector<vector<double>>>& operator *= (vector<vector<vector<double>>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();
	size_t dimz = a[0][0].size();


	for (size_t cell = 0; cell < dimx; ++cell)
	for (size_t shape = 0; shape < dimy; ++shape)
	for (size_t val = 0; val < dimz; ++val)
		a[cell][shape][val] *= b;
	return a;
}

vector<vector<double>> operator * (const vector<vector<double>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();

	vector<vector<double>> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
	for (size_t val = 0; val < dimy; ++val)
		m[cell][val] *= b;
	return m;
};


vector<vector<vector<double>>> operator * (const vector<vector<vector<double>>>& a, const double b)
{
	size_t dimx = a.size();
	size_t dimy = a[0].size();
	size_t dimz = a[0][0].size();

	vector<vector<vector<double>>> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
	for (size_t shape = 0; shape < dimy; ++shape)
	for (size_t val = 0; val < dimz; ++val)
		m[cell][shape][val] *= b;
	return m;

}

//Домножение вектора на число
vector<double>& operator *= (vector<double>& a, const double b)
{
	size_t dimx = a.size();

	for (size_t cell = 0; cell < dimx; ++cell)
		a[cell] *= b;
	return a;
};

vector<double>& operator += (vector<double>& a, const vector<double>& b)
{
	size_t dimx = a.size();

	for (size_t cell = 0; cell < dimx; ++cell)
		a[cell] += b[cell];
	return a;
};

vector<double>& operator -= (vector<double>& a, const vector<double>& b)
{
	size_t dimx = a.size();

	for (size_t cell = 0; cell < dimx; ++cell)
		a[cell] -= b[cell];
	return a;
};

vector<double> operator * (const vector<double>& a, const double b)
{
	size_t dimx = a.size();

	vector<double> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
		m[cell] *= b;
	return m;
};

vector<double> operator / (const vector<double>& a, const double b)
{
	size_t dimx = a.size();

	vector<double> m(a);

	for (size_t cell = 0; cell < dimx; ++cell)
		m[cell] /= b;
	return m;
};


//Умножение матрицы на вектор
void prodMatrVec(const vector<vector<double>>& A, \
    const vector<double>& b, \
    vector<double>& c)
{
	size_t dim = A.size();
	for (size_t row = 0; row < dim; ++row)
    {
        c[row] = 0.0;
        for (size_t col = 0; col < dim; ++col)
            c[row] += A[row][col] * b[col];
    }
}

//Умножение матриц из собственных векторов и собств. чисел (для КИР)
void prodWrAbsLWl(const vector<vector<double>>& Wr, \
    const vector<vector<double>>& Wl, \
    const vector<double>& L, \
    vector<vector<double>>& Prod)
{
	size_t dim = Wr.size();

	for (size_t row = 0; row < dim; ++row)
	for (size_t col = 0; col < dim; ++col)
	{
		Prod[row][col] = 0.0;
	}
	
	for (size_t k = 0; k < dim; ++k)
	for (size_t row = 0; row < dim; ++row)
    for (size_t col = 0; col < dim; ++col) 
    {         
            Prod[row][col] += Wr[row][k] * fabs(L[k]) * Wl[k][col];
    }
}

