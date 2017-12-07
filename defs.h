#ifndef DEFS_H_
#define DEFS_H_

#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <functional>

#include "numvector.h"
#include "Problem.h"

using namespace std;

const int dim = 5 * nShapes;

vector<numvector<double, dim>> operator * (const vector<numvector<double, dim>>& a, const double b);


vector<numvector<double, dim>> operator + (const vector<numvector<double, dim>>& b, const vector<numvector<double, dim>>& a);


void sum (const vector<numvector<double, dim>>& a, const vector<numvector<double, dim>>& b, vector<numvector<double, dim>>& res);


// ---------------------------------------------

//Прибавление к одному тройному массиву другого тройного массива
vector<vector<vector<double>>>& operator += (vector<vector<vector<double>>>& a, const vector<vector<vector<double>>>& b);

//Прибавление к одному двумерному массиву другого двумерного массива
vector<vector<double>>& operator += (vector<vector<double>>& a, const vector<vector<double>>& b);

//Вычитание из одного двумерного массива другого двумерного массива
vector<vector<double>>& operator -= (vector<vector<double>>& a, const vector<vector<double>>& b);

//Домножение трехмерного массива на число
vector<vector<vector<double>>>& operator *= (vector<vector<vector<double>>>& a, const double b);
vector<vector<vector<double>>> operator * (const vector<vector<vector<double>>>& a, const double b);

//Домножение двумерного массива на число
vector<vector<double>>& operator *= (vector<vector<double>>& a, const double b);
vector<vector<double>> operator * (const vector<vector<double>>& a, const double b);

//Домножение вектора на число
vector<double>& operator *= (vector<double>& a, const double b);
vector<double> operator * (const vector<double>& a, const double b);
vector<double> operator / (const vector<double>& a, const double b);

//Прибавление к одному вектору второго
vector<double>& operator += (vector<double>& a, const vector<double>& b);

//Вычитание из одного вектора второго
vector<double>& operator -= (vector<double>& a, const vector<double>& b);


//Умножение матрицы на вектор
void prodMatrVec(const vector<vector<double>>& A, \
    const vector<double>& b, \
    vector<double>& c);

//Умножение матриц из собственных векторов и собств. чисел (для КИР)
void prodWrAbsLWl(const vector<vector<double>>& Wr, \
    const vector<vector<double>>& Wl, \
    const vector<double>& L, \
    vector<vector<double>>& Prod);

#endif
