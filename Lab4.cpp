#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <ctime>

using std::cout;
using std::cin;
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

	Matrix() : rows(0), cols(0), ptr(nullptr) {
		cout << "default ctor running for " << this << '\n';
	}

	void allocate_memory(int _rows, int _cols) {
		ptr = new double* [_rows];
		// cout << "Memory at " << ptr << " allocated\n";
		for (int i = 0; i < _rows; ++i) {
			ptr[i] = new double[_cols];
			// cout << "Memory at " << ptr[i] << " allocated\n";
			for (int j = 0; j < _cols; ++j)
				ptr[i][j] = 0;
		}
		rows = _rows; cols = _cols;
	}
	void delete_memory() {
		// cout << "Destr running for " << this << '\n';
		if (this->ptr != nullptr) {
			for (int i = 0; i < rows; ++i) {
				// cout << "Memory at " << ptr[i] << " deleted\n";
				delete[] ptr[i];
			}
			delete[] ptr;
			// cout << "Memory at " << ptr << " deleted\n";
		}
	}

	Matrix(int _rows, int _cols) {
		// cout << "Construtor (_rows, _cols) running for " << this << '\n';
		rows = _rows;
		cols = _cols;
		ptr = new double* [rows];
		// cout << "Memory at " << ptr << " allocated\n";

		for (int i = 0; i < rows; ++i) {
			ptr[i] = new double[cols];
			// cout << "Memory at " << ptr[i] << " allocated\n";
			for (int j = 0; j < cols; ++j)
				ptr[i][j] = 0;
		}
	}

	void read_from_file(const string& file_name) {

		std::ifstream fin(file_name);
		if (fin.is_open()) {
			string line, word;

			getline(fin, line, ' ');
			int _rows = std::stod(line);
			getline(fin, line, '\n');
			int _cols = std::stod(line);
			if (!((rows == _rows) && (cols == _cols)))  // если с новыми размерами, то
			{
				this->delete_memory();
				this->allocate_memory(_rows, _cols);
			}

			int i = 0;
			int j = 0;
			while (getline(fin, line, '\n')) {
				std::istringstream stream(line);
				while (getline(stream, word, ' ')) {
					ptr[i][j] = std::stod(word);
					j++;
				}
				i++;
				j = 0;
			}
		}
		else
			cout << "Error: file wasn't opened\n";
	}

	Matrix(const string& file_name) { // иницализировать матрицу файлом
		cout << "Constructor(file) running for " << this << '\n';
		std::ifstream fin(file_name);
		if (fin.is_open()) {
			string line, word;

			/* размеры матрицы */
			getline(fin, line, ' ');
			rows = std::stod(line);
			getline(fin, line, '\n');
			cols = std::stod(line);
			//cout << "Rows = " << rows << ", Cols = " << cols << '\n';

			ptr = new double* [rows];
			//cout << "Memory at " << ptr << " allocated\n";
			for (int i = 0; i < rows; ++i) {
				ptr[i] = new double[cols];
				//cout << "Memory at " << ptr[i] << " allocated\n";
				for (int j = 0; j < cols; ++j)
					ptr[i][j] = 0;
			}

			int i = 0;
			int j = 0;
			while (getline(fin, line, '\n')) {
				std::istringstream stream(line);
				//cout << "Reading line: " << line << '\n';  // разделитель входит или нет?
				while (getline(stream, word, ' ')) {
					//cout << "Reading word: " << word << '\n';
					ptr[i][j] = std::stod(word);
					j++;
				}
				i++;
				j = 0;
			}
		}
		else
			cout << "Error: file wasn't opened\n";
	}

	Matrix(const std::initializer_list < std::initializer_list <double> >& lst) {
		// cout << "Construtor with init list running for " << this << '\n';
		rows = lst.size();
		ptr = new double* [rows];
		// cout << "Memory at " << ptr << " allocated\n";
		int i = 0;
		for (auto str : lst) {
			ptr[i] = new double[str.size()];
			// cout << "Memory at " << ptr[i] << " allocated\n";
			int j = 0;
			for (auto el : str) {
				cols = str.size();
				ptr[i][j] = el;
				j++;
			}
			i++;
		}
	}


	Matrix(const Matrix& rhs) {
		//cout << "Copy constructor running for " << this << ", rhs = " << &rhs << '\n';
		rows = rhs.rows;
		cols = rhs.cols;
		ptr = new double* [rows];
		//cout << "Memory at " << ptr << " allocated\n";

		for (int i = 0; i < rows; ++i) {
			ptr[i] = new double[cols];
			//cout << "Memory at " << ptr[i] << " allocated\n";
			for (int j = 0; j < cols; j++)
				ptr[i][j] = rhs.ptr[i][j];
		}
	}

	Matrix(Matrix&& rhs) {
		//cout << "Copy move constructor running for " << this << ", rhs = " << &rhs << '\n';
		ptr = rhs.ptr;
		rows = rhs.rows;
		cols = rhs.cols;
		rhs.ptr = nullptr;
	}

	void init() {
		for (int i = 0; i < rows; ++i) {
			cout << "Enter " << i << " string:\n";
			for (int j = 0; j < cols; ++j)
				cin >> ptr[i][j];
		}
	}

	void operator=(const Matrix& rhs) {
		// cout << "Im in op= for " << this << ", rhs = " << &rhs << '\n';
		/* !удаление старой памяти*/
		for (int i = 0; i < rows; ++i) {
			//cout << "Memory at " << ptr[i] << " deleted\n";
			delete[] ptr[i];
		}
		//cout << "Memory at " << ptr << " deleted\n";
		delete[] ptr;

		rows = rhs.rows;
		cols = rhs.cols;
		ptr = new double* [rows];
		//cout << "Memory at " << ptr << " allocated\n";
		for (int i = 0; i < rows; ++i) {
			ptr[i] = new double[cols];
			//cout << "Memory at " << ptr[i] << " allocated\n";
			for (int j = 0; j < cols; j++)
				ptr[i][j] = rhs.ptr[i][j];
		}
	}

	void operator=(Matrix&& rhs) {
		// cout << "Im in move for " << this << ", rhs = " << &rhs << '\n';
		/* удалить старую память*/
		for (int i = 0; i < rows; ++i) {
			// cout << "Memory at " << ptr[i] << " deleted\n";
			delete[] ptr[i];
		}
		// cout << "Memory at " << ptr << " deleted\n";
		delete[] ptr;

		ptr = rhs.ptr;
		rows = rhs.rows; cols = rhs.cols;
		rhs.ptr = nullptr;
	}

	double operator() (int i, int j) const {
		return ptr[i][j];
	}
	double& operator() (int i, int j) {
		return ptr[i][j];
	}

	double& operator() (int i) {
		if (rows == 1)
			return ptr[0][i];
		if (cols == 1)
			return ptr[i][0];
		else
		{
			throw std::invalid_argument("Think what ya doing when you do it!");
		}
	}
	double operator() (int i) const {
		if (rows == 1)
			return ptr[0][i];
		if (cols == 1)
			return ptr[i][0];
		else
		{
			throw std::invalid_argument("Think what ya doing when you do it!");
		}
	}

	Matrix operator-() const {
		Matrix res = *this;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				res.ptr[i][j] *= -1;
		return res;
	}

	void operator-= (const Matrix& rhs) {
		if ((rows != rhs.rows) or (cols != rhs.cols))
			throw std::invalid_argument("Dimensions error!");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				ptr[i][j] -= rhs.ptr[i][j];
		}
	}
	void operator+= (const Matrix& rhs) {
		if ((rows != rhs.rows) or (cols != rhs.cols))
			throw std::invalid_argument("Dimensions error!");
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				ptr[i][j] += rhs.ptr[i][j];
		}
	}

	Matrix operator*(const double a) const {
		Matrix res(rows, cols);
		for (int i = 0; i < res.rows; i++) {
			for (int j = 0; j < res.cols; j++) {
				res.ptr[i][j] = a * ptr[i][j];
			}
		}
		return res;
	}

	Matrix operator*(const Matrix& rhs) const {
		if (cols != rhs.rows)
			throw std::invalid_argument("Dimensions should correlate!");
		else {
			Matrix res(rows, rhs.cols);
			for (int i = 0; i < res.rows; i++) {
				for (int j = 0; j < res.cols; j++) {
					for (int t = 0; t < cols; t++)
						res.ptr[i][j] += this->ptr[i][t] * rhs.ptr[t][j];
				}
			}
			return res;
		}
	}

	void get_minor(int _row, int _col) {
		delete [] ptr[_row];
		for (int i = _row; i < rows - 1; i++) ptr[i] = ptr[i + 1];
		rows--;
		for (int i = 0; i < rows; i++) {
			for (int j = _col; j < cols - 1; j++) ptr[i][j] = ptr[i][j + 1];
		}
		cols--;
	}

	/* геттеры */
	int get_rows() const {
		return rows;
	}
	int get_cols() const {
		return cols;
	}

	/* преобразования строк*/
	void sub(int i, int j) {   // не забудь правильно проиндексировать!
		for (int t = 0; t < cols; ++t) {
			ptr[j][t] = ptr[j][t] - ptr[i][t];
		}
	}
	void scale(int i, int _factor) {
		for (int j = 0; j < cols; ++j)
			ptr[i][j] *= _factor;
	}
	void replace(int i, int j) {
		double* tmp;
		tmp = ptr[i];
		ptr[i] = ptr[j];
		ptr[j] = tmp;
	}

	void transpose() {
		if ((rows != cols) or (cols == 0) or (rows == 0))
			throw std::invalid_argument("Transpose works for square matrices only");
		for (int i = 0; i < rows; i++)
			for (int j = i + 1; j < cols; j++)
				std::swap(ptr[i][j], ptr[j][i]);
	}

	double norma_l2() const {
		double sum(0);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				sum += ptr[i][j] * ptr[i][j];
		}
		return sqrt(sum);
	}

	double norma_l1() {
		double max = 0;
		double sum = 0;
		for (int i = 0; i < rows; i++) { // цикл по строчкам
			for (int j = 0; j < cols; j++) {
				sum += abs(this->ptr[i][j]);
			}
			if (sum > max)
				max = sum;
			sum = 0;
		}
		return max;
	}
	double norma_inf() {
		double max = 0;
		double sum = 0;
		for (int j = 0; j < cols; j++) {
			for (int i = 0; i < rows; i++)
				sum += abs(this->ptr[i][j]);
			if (sum > max)
				max = sum;
			sum = 0;
		}
		return max;
	}



	//Matrix inversed() {
	//	if ((rows != cols) or (cols == 0) or (rows == 0))
	//		throw std::invalid_argument("Inverse works for square matrices only");
	//	int size = rows;
	//	Matrix res(size, size);
	//	Matrix rhs(1, size); rhs(0) = 1;
	//	//Matrix _this = *this; 
	//   // Matrix slv(1, size);

	//	for (int i = 0; i < size; i++) {
	//		for (int j = 0; j < size; j++)
	//			rhs(j) = 0;
	//		rhs(i) = 1;
	//		//_this = *this; 
	//		rhs = Gauss_by_value(*this, rhs, 1e-9);
	//		//cout << "Slv = " << rhs;
	//		for (int j = 0; j < size; j++)
	//			res(j, i) = rhs(j);
	//	}
	//	return res;
	//}
	// Matrix inversed_via_QR(){
	//     if ( (rows != cols) or (cols == 0) or (rows == 0))
	//         throw std::invalid_argument("Inverse works for square matrices only");
	//     int size = rows;
	//     Matrix res(size, size);
	//     Matrix rhs(1, size); rhs(0) = 1;

	//     Matrix Q(size, size);
	//     Matrix _this(*this);
	//     QR(_this, Q);
	// }

	void write_to_file(std::ostream& out) {
		if (cols == 1) {
			out << ptr[0][0];
			for (int k = 1; k < rows; k++)
				out << ";" << ptr[k][0];
			out << '\n';
			return;
		}
		for (int i = 0; i < rows; i++) {
			out << ptr[i][0];
			for (int j = 1; j < cols; j++)
				out << ';' << ptr[i][j];
			out << '\n';
		}
	}

	~Matrix() {
		// cout << "Destr running for " << this << '\n';
		if (this->ptr != nullptr) {
			for (int i = 0; i < rows; ++i) {
				// cout << "Memory at " << ptr[i] << " deleted\n";
				delete[] ptr[i];
			}
			delete[] ptr;
			// cout << "Memory at " << ptr << " deleted\n";
		}
	}

};

Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
	if ((lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
		std::invalid_argument("Dimensions error!");
	Matrix res = rhs;
	for (int i = 0; i < lhs.rows; i++) {
		for (int j = 0; j < lhs.cols; j++)
			res.ptr[i][j] += lhs.ptr[i][j];
	}
	return res;
}
Matrix operator+(const Matrix& lhs, Matrix&& rhs) {
	if ((lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
		throw std::invalid_argument("Dimensions error!");
	Matrix res = std::move(rhs);
	for (int i = 0; i < lhs.rows; i++) {
		for (int j = 0; j < lhs.cols; j++)
			res.ptr[i][j] += lhs.ptr[i][j];
	}
	return res;
}
Matrix operator+(Matrix&& lhs, const Matrix& rhs) {
	if ((lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
		throw std::invalid_argument("Dimensions error!");
	Matrix res = std::move(lhs);
	for (int i = 0; i < lhs.rows; i++) {
		for (int j = 0; j < lhs.cols; j++)
			res.ptr[i][j] += rhs.ptr[i][j];
	}
	return res;
}
Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
	return lhs + (-rhs);
}
Matrix operator+(const Matrix& lhs, const std::initializer_list < std::initializer_list <double> >& rhs) {
	Matrix res = lhs;
	int i(0), j(0);
	for (auto str : rhs) {
		for (auto elem_of_str : str) {
			res.ptr[i][j] += elem_of_str;
			j++;
		}
		i++;
	}
	return res;
}
void operator<<(std::ostream& out, const Matrix& x) {
	out << '\n';
	for (int i = 0; i < x.rows; ++i) {
		for (int j = 0; j < x.cols; ++j)
			out << x(i, j) << " ";
		out << '\n';
	}
	out << '\n';
}
Matrix eye(int size) {
	Matrix res(size, size);
	for (int i = 0; i < size; ++i)
		res(i, i) = 1;
	return res;
}
Matrix zeros(int _rows, int _cols) {
	Matrix res(_rows, _cols);
	return res;
}

Matrix transpose_v(const Matrix& V) { //proverka!!
	if (V.get_rows() == 1) {
		int size = V.get_cols();
		Matrix v(size, 1);
		for (int i = 0; i < size; i++) v(i) = V(i);
		return v;
	}
	if (V.get_cols() == 1) {
		int size = V.get_rows();
		Matrix v(1, size);
		for (int i = 0; i < size; i++) v(i) = V(i);
		return v;
	}
	throw std::invalid_argument("Think what ya doing when you do it!");
}

double dot(const Matrix& V1, const Matrix& V2) { //proverka!!
	int size;
	if (V1.get_cols() > V1.get_rows()) size = V1.get_cols();
	else size = V1.get_rows();
	double scal(0);
	for (int i = 0; i < size; i++) scal += V1(i) * V2(i);
	return scal;
}

const double eps = 1e-6;
const double eps_hess = eps;
const double eps_gauss = eps;
const double eps_QR = eps;
const double eps_inv_it = eps;

int find_leading_element(const Matrix& A, int k) {
	double max = abs(A(k, k));
	int row = k;
	for (int i = k + 1; i < A.get_rows(); i++) {
		if (abs(A(i, k)) > max) {
			max = abs(A(i, k));
			row = i;
		}
	}
	return row;
}

void gauss_back_run(const Matrix& A, const Matrix& b, Matrix& X) {
	double sum(0);
	int rows_num = A.get_rows();
	for (int j = rows_num - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
		sum = 0;
		for (int k = j + 1; k < rows_num; k++)
			sum += X(k) * A(j, k);
		X(j) = (b(j) - sum) / A(j, j);
	}
}

Matrix Gauss(Matrix& A, Matrix& b, const double eps) {
	// cout << "In Gauss\nArguments passed:\nA = " << A;
	// cout << "b = " << b;

	int rows_num = A.get_rows();  // а может всегда вызывать get_...()? сколько занимает вызов метода (функции?)
	int cols_num = A.get_cols();
	double tmp1(0), tmp2(0);  // инициализировать нулями?

	/*  Прямой ход (c частичным выбором) */
	for (int k = 0; k < (rows_num - 1); k++) {
		//tmp1 = A(k,k);
		int row_of_leading_element = find_leading_element(A, k);
		/* проверк на вырожденность с заданной точностью eps */
		if (abs(A(row_of_leading_element, k)) < eps)
			throw std::invalid_argument("Given Matrix is non-invertable");

		A.replace(k, row_of_leading_element);
		/* не забыть поменять строки вектора правой части */
		std::swap(b(k), b(row_of_leading_element));

		tmp1 = A(k, k);
		for (int i = k + 1; i < rows_num; i++) {
			tmp2 = A(i, k) / tmp1;

			for (int j = k; j < cols_num; j++) {
				A(i, j) = A(i, j) - A(k, j) * tmp2;
			}

			/* соответствующим образом меняем вектор правой части */
			b(i) = b(i) - b(k) * tmp2;
		}
	}
	// cout << "После прямого хода:\nA= " << A;
	// cout << "b = " << b;

	/* Обратный ход */
	Matrix X(rows_num, 1);
	gauss_back_run(A, b, X);
	return X;
}

void Hessenberg(Matrix& A) {
	int size = A.get_rows();
	for (int k = 1; k < size - 1; k++) {
		for (int l = k + 1; l < size; l++) {
			double root = sqrt(A(k, k - 1) * A(k, k - 1) + A(l, k - 1) * A(l, k - 1));
			double alpha = A(k, k - 1) / root;
			double beta = A(l, k - 1) / root;
			if (alpha * alpha + beta * beta - 1 > eps_hess) {
				cout << "alpha^2 + beta^2 = " << alpha * alpha + beta * beta << "\n";
				throw std::invalid_argument("Trigonometric identity doesn't work!");
			}
			for (int j = 0; j < size; j++) {
				double akj = A(k, j);
				A(k, j) = alpha * A(k, j) + beta * A(l, j);
				if (j != k - 1) {
					A(l, j) = alpha * A(l, j) - beta * akj;
				}
			}
			for (int i = 0; i < size; i++) {
				double aik = A(i, k);
				A(i, k) = alpha * A(i, k) + beta * A(i, l);
				A(i, l) = alpha * A(i, l) - beta * aik;
			}
			A(l, k - 1) = 0;
		}
	}
}

Matrix QR(Matrix& A, bool Hess, int& rotations) {
	int size = A.get_rows();
	Matrix str(1, size);
	Matrix Q = eye(size);
	int lim = size;
	rotations = 0;
	for (int k = 0; k < size - 1; k++) {
		if (Hess) lim = k + 2;
		for (int i = k + 1; i < lim; i++) {
			double norm = sqrt(A(k, k) * A(k, k) + A(i, k) * A(i, k));
			double c = A(k, k) / norm;
			double s = A(i, k) / norm;

			for (int t = 0; t < size; t++)
				str(t) = Q(k, t);
			for (int j = 0; j < size; j++)
				Q(k, j) = c * str(j) + s * Q(i, j);
			for (int j = 0; j < size; j++)
				Q(i, j) = -s * str(j) + c * Q(i, j);

			for (int t = k; t < size; t++)
				str(t) = A(k, t);
			for (int j = k; j < size; j++)
				A(k, j) = c * str(j) + s * A(i, j);
			for (int j = k; j < size; j++)
				A(i, j) = -s * str(j) + c * A(i, j);
			rotations += 1;
		}
	}
	Q.transpose();
	return Q;
} 

bool condition_for_reduction_shift(const Matrix& A) { //l2 modified for row k left from diag
	double sum(0);
	int row = A.get_rows() - 1;
	for (int i = 0; i < row; i++) {
		sum += A(row, i) * A(row, i);
	}
	if (sqrt(sum) >= eps_QR) return true;
	return false;
}

bool condition_for_reduction(const Matrix& A, const Matrix& A_prev, int iter) {
	if (iter == 0) return true;
	if ((abs(A(0,0)) - abs(A_prev(0,0))) >= eps_QR) return true;
	return false;
}

Matrix eigenvalues_QR_with_reduction(Matrix A, bool shift, bool Hess)
{
	int total_rotations(0);
	int total_iter(0);
	const int max_iter = 1e+3;
	int size = A.get_rows();
	Matrix eigenvalues(1, size);
	if (Hess) Hessenberg(A);
	Matrix A_prev = A;
	for (int k = size - 1; k > 0; --k) {
		Matrix E = eye(k + 1);
		int sum_rotations(0);
		int iter(0);
		if (shift) {
			while (condition_for_reduction_shift(A)) {
				double sigma = A(k, k);
				A = A - E * sigma;
				int rotations(0);
				Matrix Q = QR(A, Hess, rotations);
				sum_rotations += rotations;
				A = A * Q;
				A = A + E * sigma;
				if (iter == max_iter) throw std::invalid_argument("Something went wrong in QR_eigenvalues...");
				//cout << "Q: " << Q;
				//cout << "\n";
				iter++;
			}
		}
		else {
			while (condition_for_reduction(A, A_prev, iter)) {
				A_prev = A;
				int rotations(0);
				Matrix Q = QR(A, Hess, rotations);
				sum_rotations += rotations;
				A = A * Q;
				if (iter == max_iter) throw std::invalid_argument("Something went wrong in QR_eigenvalues...");
				//cout << "Q: " << Q;
				//cout << "\n";
				iter++;
			}
		}
		//cout << "A: " << A;
		total_rotations += sum_rotations;
		total_iter += iter;
		if (shift) eigenvalues(0, k) = A(k, k);
		else eigenvalues(0, k) = A(0, 0);
		if (shift) A.get_minor(k, k);
		else A.get_minor(0, 0);
	}
	cout << "The number of iterations: " << total_iter << "\n";
	cout << "The number of total rotations: " << total_rotations << "\n";
	eigenvalues(0, 0) = A(0, 0);
	return eigenvalues;
}

//double calc_delta(const Matrix& A) { //l2 modified for all rows left from diag
//	double sum(0);
//	for (int row = A.get_rows() - 1; row > 0; --row) {
//		for (int i = 0; i < row; i++) sum += A(row, i) * A(row, i);
//	}
//	return sqrt(sum);
//}

//double calc_sigma(const Matrix& A, const Matrix& A_prev, double sigma_pred) {
//	double max(0);
//	int row_index(0);
//	for (int i = A.get_rows() - 1; i > 0; --i) {
//		for (int j = 0; j < i; j++) {
//			double n = A(i, j) / abs(abs((A(i, j)) - abs(A_prev(i, j))));
//			if (n > max) {
//				max = n;
//				row_index = i;
//			}
//		}
//	}
//	if (max < 1) return sigma_pred;
//	return A(row_index, row_index);
//}

// DO NOT USE
//Matrix eigenvalues_QR(Matrix A, bool shift, bool Hess)
//{
//	int total_rotations(0);
//	const int max_iter = 1e+3;
//	int size = A.get_rows();
//	Matrix eigenvalues(1, size);
//	if (Hess) Hessenberg(A);
//	Matrix E = eye(size);
//	int iter = 0;
//	Matrix A_prev = A;
//	double sigma(0);
//	while (calc_delta(A) > eps_QR) {
//		if (shift) A_prev = A;
//		if (shift) A = A - E * sigma;
//		int rotations(0);
//		Matrix Q = QR(A, Hess, rotations);
//		total_rotations += rotations;
//		A = A * Q;
//		if (shift) A = A + E * sigma; 
//		if (shift) sigma = calc_sigma(A, A_prev, sigma);
//		if (iter == max_iter) throw std::invalid_argument("Something went wrong in QR_eigenvalues...");
//		iter++;
//		//cout << "Q: " << Q;
//		//cout << "\n";
//	}
//	cout << "The number of iterations: " << iter << "\n";
//	cout << "The number of total rotations: " << total_rotations << "\n";
//	cout << "A: " << A;
//	for (int j = 0; j < size; j++) eigenvalues(j) = A(j, j);
//	return eigenvalues;
//}

Matrix sort(Matrix& v) {
	int size = v.get_cols();
	for (int i = 0; i < size - 1; i++)
	for (int j = i + 1; j < size; j++) {
		if (v(0, i) > v(0, j)) {
			double t = v(0, i);
			v(0, i) = v(0, j);
			v(0, j) = t;
		}
	}
	return v;
}

double relay_func(const Matrix& A, const Matrix& x) {
	double lambda(0);
	lambda = dot(A * x, x);
	return lambda;
}

Matrix inverse_iterations(const Matrix& A, const double eigenvalue_QR, bool relay, Matrix& x) {
	int size = A.get_rows();
	Matrix E = eye(size);
	x = x * (1 / x.norma_l2());
	double eigenvalue;
	if (relay) eigenvalue = relay_func(A, x);
	else eigenvalue = eigenvalue_QR;
	Matrix lambda_E = E * eigenvalue;
	Matrix M = A - lambda_E;
	while ((M * x).norma_l2() > eps_inv_it) {
		x = Gauss(M, x, eps_gauss);
		x = x * (1 / x.norma_l2());
		if (relay) lambda_E = E * relay_func(A, x);
		M = A - lambda_E;
	}
	if (relay) cout << "Eigenvalue: " << relay_func(A, x);
	return x;
}


int main() {

	Matrix A { {1.5,0,-0.43,-0.75}, {0,3,0.87,-0.5}, {-0.43,0.87,2.9,-0.22}, {-0.75,-0.5,-0.22,2.6}};
	int size = A.get_rows();
	Matrix e {{ 0.997313, 2.00425, 2.98702, 4.01142 }}; //right answer: eigenvalues
	cout << "A: " << A;

	cout << "\n------------QR EIGENVALUES--------------\n\n";//
	cout << "\n----------without Hessenberg-----------\n";// 
	Matrix eigenvalues = eigenvalues_QR_with_reduction(A, true, false);
	cout << "\nEigenvalues whithout Hessenberg: " << sort(eigenvalues);
	cout << "Accuracy: " << (e - sort(eigenvalues)).norma_l2();
	cout << "\n";

	cout << "\n------------with Hessenberg------------\n";// 
	eigenvalues = eigenvalues_QR_with_reduction(A, true, true);
	cout << "\nEigenvalues with Hessenberg: " << sort(eigenvalues);
	cout << "Accuracy: " << (e - sort(eigenvalues)).norma_l2();
	cout << "\n\n\n";

	Matrix x(size,1); // first approximation

	cout << "\n-------------EIGENVECTORS-------------\n\n";
	cout << "\n----------of QR eigenvalues-----------\n";
	for (int i = 0; i < size; i++) x(i) = 1;
	for (int i = 0; i < A.get_rows(); i++) {
		Matrix eigenvector = inverse_iterations(A, eigenvalues(i), false, x);
		cout << "Eigenvalue: " << eigenvalues(i);
		cout << "\nEigenvector: " << transpose_v(eigenvector);
	}

	cout << "\n---------of Relay eigenvalues---------\n";
	x = { {-0.6}, {0.3}, {-0.1}, {-0.2} }; //1
	//x = { {0.1}, {0.6}, {-0.5}, {0.2} }; //2
	//x = { {0.4}, {0}, {-0.3}, {-0.6} }; //3
	//x = { {-0.1}, {0.5}, {0.3}, {-0.2} }; //4
	cout << "\nFor the first approximation: " << transpose_v(x);
	Matrix eigenvector = inverse_iterations(A, 0, true, x);
	cout << "\nEigenvector: " << transpose_v(eigenvector);

	return 0;
}