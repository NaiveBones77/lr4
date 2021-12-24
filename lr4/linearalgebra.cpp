#include "linearalgebra.h"

vector::vector(){
    this->resize(0);
}

vector::vector(int size){
    this->resize(size);
}

vector::vector(const vector & rVector) {
    this->data.assign(rVector.data.begin(), rVector.data.end());
}

void vector::resize(int size) {
    assert(size >= 0);
    this->data.resize(size);
}

int vector::size() const {
    return this->data.size();
}

void vector::print() const {
    std::cout << "["; std::copy(this->data.begin(), this->data.end(), std::ostream_iterator<long double>(std::cout, " ")); std::cout << "]" << std::endl;
}

long double vector::length() const
{
    long double length{0.0L};
        for (int i = 0; i < this->size(); ++i)
            length += pow(this->data[i], 2.0L);
        return sqrt(length);
}

vector & vector::operator = (const vector & rVector)
{
    this->data = rVector.data;
    return *this;
}
vector vector::operator + (const vector & rVector) const
{
    vector resVector;
    std::transform(this->data.begin(), this->data.end(),
                   rVector.data.begin(), std::back_inserter(resVector.data),
                   std::plus<long double>());
    return resVector;
}
vector vector::operator - (const vector & rVector) const
{
    vector resVector;
    std::transform(this->data.begin(), this->data.end(),
                   rVector.data.begin(), std::back_inserter(resVector.data),
                   std::minus<long double>());
    return resVector;
}

vector vector::operator^(const vector & rVector) const {
    assert(rVector.size() == 3);
    vector resVector(3);
    resVector[0] = this->operator[](1) * rVector[2] - rVector[1] * this->operator[](2);
    resVector[1] = this->operator[](2) * rVector[0] - rVector[2] * this->operator[](0);
    resVector[2] = this->operator[](0) * rVector[1] - rVector[0] * this->operator[](1);
    return resVector;
}

long double vector::operator*(const vector & rVector) const {
    assert(rVector.size() == size());
    long double res = 0;
    for (int i = 0; i < rVector.size(); i++)
    {
        res += this->operator[](i) * rVector[i];
    }
    return res;
}

vector vector::operator * (long double rValue) const {
    vector resVector(this->size());
    std::vector<long double>::const_iterator it;
    it = this->data.begin();
    int i{0};
    while (it != this->data.end()) {
        resVector[i] = *it * rValue;
        ++i;
        ++it;
    }
    return resVector;
}


vector::~vector(){
}

matrix::matrix(){
    this->resize(0, 0);
}

matrix::matrix(int rows, int columns){
    this->resize(rows, columns);
}

matrix::matrix(const matrix & rMatrix){
    this->data.assign(rMatrix.data.begin(), rMatrix.data.end());
}

void matrix::resize(int rows, int columns){
    this->data.resize(rows);
    for (int i = 0; i < this->rowsCount(); ++i) {
        this->data.at(i).resize(columns);
    }
}

int matrix::rowsCount() const{
    return this->data.size();
}

int matrix::colsCount() const{
    return this->data.front().size();
}

long double matrix::det()
{
    const long double EPS {1e-9};
    assert(this->rowsCount() == this->colsCount());
    long double det{1.0};
    for (int i = 0; i < this->rowsCount(); ++i) {
        int k = i;
        for (int j = i + 1; j < this->rowsCount(); ++j)
            if (this->operator()(i,j) > this->operator()(k, i))
                k = j;
        if (std::abs(this->operator()(k, i)) < EPS) {
            det = 0;
            break;
        }
        this->swapRows(i, k);
        if (i != k)
            det = -det;
        det *= this->operator()(i, i);
        for (int j = i + 1; j < this->rowsCount(); ++j)
            this->operator()(i, j) /= this->operator()(i, i);
        for (int j = 0; j < this->rowsCount(); ++j)
            if (j != i && std::abs(this->operator()(j, i)) > EPS)
                for (int k = i + 1; k < this->rowsCount(); ++k)
                    this->operator()(j, k) -= this->operator()(i, k) * this->operator()(j, i);
    }

    return det;
}

matrix & matrix::operator = (const matrix & rMatrix)
{
    this->data.assign(rMatrix.data.begin(), rMatrix.data.end());
    return *this;
}

matrix matrix::operator + (const matrix & rMatrix) const
{
    assert( ( rowsCount() == rMatrix.rowsCount()) && (colsCount() == rMatrix.colsCount()));
    matrix resMatrix;
    for (int i = 0; i < this->rowsCount(); ++i) {
        for (int j = 0; j < this->colsCount(); ++j) {
            resMatrix(i, j) = this->operator() (i, j) + rMatrix(i, j);
        }
    }
    return resMatrix;
}

matrix matrix::operator-(const matrix & rMatrix) const
{
    assert( (rowsCount() == rMatrix.rowsCount()) && (colsCount() == rMatrix.colsCount()) );
    matrix resMatrix;
    for (int i = 0; i < this->rowsCount(); ++i) {
        for (int j = 0; j < this->colsCount(); ++j) {
            resMatrix(i, j) = this->operator() (i, j) + rMatrix(i, j);
        }
    }
    return resMatrix;
}

matrix matrix::operator*(const matrix & rMatrix) const
{
    assert(this->colsCount() == rMatrix.rowsCount());
    matrix resMatrix(this->rowsCount(), rMatrix.colsCount());
    for (int i = 0; i < resMatrix.rowsCount(); i++)
    {
        for (int j = 0; j < resMatrix.colsCount(); j++)
        {
            resMatrix(i, j) = 0;
            for (int k = 0; k < this->colsCount(); k++)
            {
                resMatrix(i, j) += this->operator()(i, k) * rMatrix(k, j);
            }
        }
    }
    return resMatrix;
}

matrix matrix::operator*(long double rValue) const
{

    matrix resMatrix(this->rowsCount(), this->colsCount());
    for (int i = 0; i < this->rowsCount(); ++i) {
        for (int j = 0; j < this->colsCount(); ++j) {
            resMatrix(i, j) = this->operator()(i, j) * rValue;
        }
    }
    return resMatrix;
}

matrix matrix::operator!()
{
    matrix R(this->rowsCount(), this->rowsCount());
    matrix M;
    M = *this;
    matrix I(this->rowsCount(), this->rowsCount());
    I = I.I();
    long double element = 0;
    for (int i = 0; i < this->rowsCount(); i++) {
        element = this->operator()(i, i);
        if (element == 0) {
            bool swapFlag = false;
            for (int s = i; s < this->rowsCount(); s++) {
                if (this->operator()(s, i) != 0) {
                    swapFlag = true;
                    M.swapRows(i, s);
                    I.swapRows(i, s);
                }
            }
            if (swapFlag == false)
                break;
        }
        element = this->operator()(i, i);
        for (int k = 0; k < this->rowsCount(); k++) {
            this->operator()(i, k) = this->operator()(i, k) / element;
            I(i, k) = I(i, k) / element;
        }
        long double c = 0;
        for (int k = 0; k < this->rowsCount(); k++) {
            if (k == i) continue;
            c = this->operator()(k, i);
            for (int p = 0; p < this->rowsCount(); p++) {
                this->operator()(k, p) = this->operator()(k, p) - c * this->operator()(i, p);
                I(k, p) = I(k, p) - c * I(i, p);
                R(k, p) = I(k, p);
            }
        }
    }
    for (int i = 0; i < this->rowsCount(); i++) {
        for (int j = 0; j < this->rowsCount(); j++) {
            R(i, j) = I(i, j);
        }
    }
    return R;
}

matrix  matrix::T()
{
    if (this->rowsCount() > this->colsCount()) {
        int temp = this->rowsCount() - this->colsCount();
        matrix resMatrix(*this);
        resMatrix.resize(rowsCount(), rowsCount());
        for (int i = 0; i < resMatrix.rowsCount(); i++) {
            for (int j = i; j < resMatrix.colsCount(); j++) {
                std::swap(resMatrix(i, j), resMatrix(j, i));
            }
        }

        resMatrix.resize(rowsCount() - temp, rowsCount());
        return resMatrix;
    }
    if (this->rowsCount() < this->colsCount()) {
        int temp = this->colsCount() - this->rowsCount();
        matrix resMatrix(*this);
        resMatrix.resize(colsCount(), colsCount());
        for (int i = 0; i < resMatrix.rowsCount(); i++) {
            for (int j = i; j < resMatrix.colsCount(); j++) {
                std::swap(resMatrix(i, j), resMatrix(j, i));
            }
        }
        resMatrix.resize(colsCount(), colsCount() - temp);
        return resMatrix;
    }
    else {
        matrix resMatrix(rowsCount(), rowsCount());
        for (int i = 0; i < resMatrix.rowsCount(); i++) {
            for (int j = 0; j < resMatrix.colsCount(); j++) {
                std::swap(resMatrix(i, j), resMatrix(j, i));
            }
        }
        return resMatrix;
    }
}

matrix & matrix::swapRows(int pos, int newPos)
{
    long double tmp_cell;
    for (int k = 0; k < this->rowsCount(); k++)
    {
        tmp_cell = this->operator()(pos, k);
        this->operator()(pos, k) = this->operator()(newPos, k);
        this->operator()(newPos, k) = tmp_cell;
    }
    return *this;
}

matrix matrix::I()
{
    matrix I(this->rowsCount(), this->rowsCount());
    for (int i = 0; i < this->rowsCount(); i++) {
        I(i, i) = 1;
        for (int j = i + 1; j < this->rowsCount(); j++)
        {
            I(i, j) = I(j, i) = 0;
        }
    }
    return I;
}

vector matrix::operator*(const vector & rVector) const
{
    assert(rVector.size() == this->colsCount());
    vector V(this->rowsCount());
    for (int i = 0; i < this->rowsCount(); i++)
    {
        V[i] = 0;
        for (int j = 0; j < this->colsCount(); j++)
            V[i] += this->operator()(i, j) * rVector[j];
    }
    return V;
}


matrix::~matrix()
{
}
