#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>

class vector
{
private:
    std::vector<long double> data;
    int elements;
public:
    vector();
    vector(int size);
    vector(const vector & rVector);
    void resize(int size);
    int size() const;
    void print() const;
    long double length() const;
    friend void print(const vector & prinvector)
    {
        std::cout << "["; std::copy(prinvector.data.begin(), prinvector.data.end(), std::ostream_iterator<long double>(std::cout, " ")); std::cout << "]"<< std::endl;
    }
    long double operator [] (int place) const { return this->data.at(place); }
    long double & operator [] (int place) { return this->data.at(place); }
    vector & operator = (const vector & rVector);
    vector operator + (const vector & rVector) const;
    vector operator - (const vector & rVector) const;
    vector operator ^ (const vector & rVector) const;
    vector operator * (long double rValue) const;
    long double operator * (const vector & rVector) const;
    ~vector();
};
class matrix
{
private:
    std::vector<std::vector<long double>> data;
    int rows, cols;
public:
    matrix();
    matrix(int rows, int columns);
    matrix(const matrix & rMatrix);
    void resize(int rows, int columns);
    int rowsCount() const;
    int colsCount() const;
    long double det();
    long double operator () (int row, int col) const { return this->data.at(row).at(col); }
    long double & operator () (int row, int col) { return this->data.at(row).at(col); }
    int getHighRow() { return this->rowsCount() - 1; }
    int getHighCol() { return this->colsCount() - 1; }
    matrix & operator = (const matrix & rMatrix);
    matrix operator + (const matrix & rMatrix) const;
    matrix operator - (const matrix & rMatrix) const;
    matrix operator * (const matrix & rMatrix) const;
    matrix operator * (long double rValue) const;
    matrix operator ! ();
    matrix  T();
    matrix& swapRows(int pos, int newPos);
    matrix I();
    vector operator * (const vector& rVector) const;
    friend void print(const matrix & printMatrix)
    {
        std::cout << "[ ";
        for (int i = 0; i < printMatrix.rowsCount(); ++i) {
            if (i != 0) std::cout << "  ";
            std::copy(printMatrix.data.at(i).begin(), printMatrix.data.at(i).end(), std::ostream_iterator<long double>(std::cout, " "));
            if(i != printMatrix.rowsCount() - 1) std::cout << std::endl;
        }
        std::cout << "]" << std::endl;
    }
    ~matrix();
};

#endif // LINEARALGEBRA_H
