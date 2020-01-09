

#include "Matrix.h"


Matrix::Matrix() : nlines(0), ncols(0), vec(0) { }

//Matrix::Matrix(int nlines, int ncols) : nlines(nlines), ncols(ncols), vec(nlines*ncols) {
Matrix::Matrix(int nlines, int ncols) : nlines(nlines), ncols(ncols) {
//    cout << *this << endl;
}

vector<int>& Matrix::getVec() { return vec; }

// Exams's indexes start at 1
int Matrix::getVal(int i, int j) const { 
//    int i0 = i, j0 = j;
    int i0 = i-1, j0 = j-1;
//    cout << "[i, j] = [" << i0 << ", " << j0 << "]" << endl;
//    assert(i0 >= 0 && i0 < nlines && j0 >= 0 && j0 < ncols);
    return vec[i0*ncols+j0];
}

int Matrix::nLines() const { return nlines; }

int Matrix::nCols() const { return ncols; }

void Matrix::setDimensions(int nlines, int ncols) {
	this->nlines = nlines;
	this->ncols = ncols;
	//vec.resize(nlines*ncols);
}

void Matrix::reserve(int nlines, int ncols) {
    this->nlines = nlines;
    this->ncols = ncols;
    vec.resize(nlines*ncols);
}

ostream& operator<<(ostream& os, const Matrix& matrix) {
	os << endl << "Matrix" << endl;
	os << "nlines = " << matrix.nLines() << ", ncols = " << matrix.nCols() << endl;
	os << "size = " << matrix.nLines()*matrix.nCols() << endl;
	os << "vec.size = " << matrix.vec.size() << endl;

    for (int i = 0; i < matrix.nLines(); ++i) {
        for (int j = 0; j < matrix.nCols(); ++j) {
            os << matrix.getVal(i, j) << " ";
        }
        os << endl;
    }
	return os;
}

