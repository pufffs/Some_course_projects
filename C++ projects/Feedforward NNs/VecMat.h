#pragma once
#include <vector>
#include<iostream>
#include<ostream>
using namespace std;
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double& operator[](int index)
	{
		return v[index];
	}

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		return v[index];
	}

	int size() const { return v.size(); } // number of elements

private:
	std::vector<double> v;
};

class MMatrix
{
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {} //default constructor.
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n* m, x) {}

	// set all matrix entries equal to a double
	MMatrix& operator=(double x)
	{
		for (unsigned i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[ j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double& operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }
	int size() const { return A.size(); }
private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};
