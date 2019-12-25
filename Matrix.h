#pragma once

#include <iostream>

using std::cout;
using std::string;

class Matrix {
	int rows, cols;
	double** ptr;

public:

	friend void operator<< (std::ostream& out, const Matrix& x);
	friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
	friend Matrix operator+(const Matrix& lhs, Matrix&& rhs);
	friend Matrix operator+(Matrix&& lhs, const Matrix& rhs);
	friend Matrix operator-(const Matrix& lhs, const Matrix& rhs); // перегрузить на мув 
	friend Matrix operator+(const Matrix& lhs, const std::initializer_list < std::initializer_list <double> >& rhs);

	Matrix();

	void allocate_memory(int _rows, int _cols);
	void delete_memory();

	Matrix(int _rows, int _cols);

	void read_from_file(const string& file_name);

	Matrix(const string& file_name);

	Matrix(const std::initializer_list < std::initializer_list <double> >& lst);

	Matrix(const Matrix& rhs);

	Matrix(Matrix&& rhs);

	void init();

	void operator=(const Matrix& rhs);

	void operator=(Matrix&& rhs);

	double operator() (int i, int j) const;

	double& operator() (int i, int j);

	double& operator() (int i);

	double operator() (int i) const;

	Matrix operator-() const;

	void operator-= (const Matrix& rhs);
	void operator+= (const Matrix& rhs);
	Matrix operator*(const double a) const;

	Matrix operator*(const Matrix& rhs) const;

	void get_minor(int _row, int _col);

	/* геттеры */
	int get_rows() const;
	int get_cols() const;

	/* преобразования строк*/
	void sub(int i, int j);
	void scale(int i, int _factor);
	void replace(int i, int j);

	void transpose();
	double norma_l2() const;

	double norma_l1();
	double norma_inf();

	void invert2x2();
	//Matrix inversed();
	//Matrix inversed_via_QR();

	void write_to_file(std::ostream& out);

	~Matrix();
};

Matrix eye(int size);
Matrix zeros(int _rows, int _cols);
Matrix transpose_v(const Matrix& V);
double dot(const Matrix& V1, const Matrix& V2);

//Matrix Gauss(Matrix& const A, Matrix& const b);
//void Hessenberg(Matrix& A);
//Matrix QR(Matrix& A, bool Hess, int& rotations);

