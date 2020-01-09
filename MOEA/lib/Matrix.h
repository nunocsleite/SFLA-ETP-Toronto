#ifndef MATRIX
#define MATRIX

#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;


class Matrix {

private:
	vector<int> vec;
	int nlines, ncols;

public:
	Matrix();
	Matrix(int nlines, int ncols);
	vector<int>& getVec();
	int getVal(int i, int j) const;
	int nLines() const;
	int nCols() const;
	void setDimensions(int nlines, int ncols);
    void reserve(int nlines, int ncols);
	friend ostream& operator<<(ostream& os, const Matrix& Matrix);

};

#endif
