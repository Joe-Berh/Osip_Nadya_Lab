#include <vector>
#include <stdio.h>
#include <fstream>

#include "Matrix.h"

using std::cout;
using std::vector;

typedef double(*Function)(double x);
typedef double(*Function2)(double x1, double x2);

struct Segment {
	double a;
	double b;
};

vector<Segment> locate_roots(Function f, Segment initial_segment, double step) 
{
	vector<Segment> segments;
	double x = initial_segment.a;
	double y_prev = f(x);
	double x_prev = x;
	x += step;
	bool last_step = false;
	do {
		if (x >= initial_segment.b) {
			x = initial_segment.b;
			last_step = true;
		}
		double y = f(x);
		if (y * y_prev <= 0) segments.push_back({ x_prev, x });
		y_prev = y;
		x_prev = x;
		x += step;
	} while (!last_step);

	return segments;
}

double solve_bisection(Function f, Segment segment, double eps, int& iter)
{
	double a = segment.a;
	double b = segment.b;
	double fa = f(a);
	double fb = f(b);

	if (a > b) throw std::invalid_argument("Incorrect order of segment ends: a > b !");

	iter = 0;

	while ((b - a) / 2 > eps) {
		if (fa * fb > 0) {
			cout << "Required accuracy cannot be achieved - \nthe function on both ends of segment has the same sign!\n";
			break;
		}
		double c = (b + a) / 2;
		double fc = f(c);
		if (fa * fc <= 0) {
			b = c;
			fb = fc;
		} 
		else {
			a = c;
			fa = fc;
		}
		++iter;
	}
	
	cout << "iterations: " << iter << "\n";
	
	return (b + a) / 2;
}

void CheckAndFit(double& x, Segment segment)
{
	if (x < segment.a) x = segment.a;
	if (x > segment.b) x = segment.b;
}

void CheckAndFit(Matrix& x, double l1, double l2)
{
	if (x(0) < -l1) x(0) = -l1;
	if (x(0) > l1) x(0) = l1;
	if (x(1) < -l2) x(1) = -l2;
	if (x(1) > l2) x(1) = l2;
}

double horde(Function f, Segment segment)
{
	double fa = f(segment.a);
	double fb = f(segment.b);
	return (fa * segment.b - fb * segment.a) / (fa - fb);
}

double diff(Function f, double x)
{
	double const dx = 1e-04;
	return (f(x + dx) - f(x)) / dx;
}

double solve_newton(Function f, Function df, Segment segment, double x0, double eps, int& iter)
{
	const int max_iter = 30;
	iter = 0;
	double x;
	double x_prev = x0;
	do {
		double dfx_prev;
		if (df == 0) dfx_prev = diff(f, x_prev);
		else dfx_prev = df(x_prev);

		x = x_prev - f(x_prev) / dfx_prev;
		++iter;

		cout << iter << ": " << x << "\n";

		if (abs(x - x_prev) <= eps) break;

		// avoid running out of segment
		CheckAndFit(x, segment);

		x_prev = x;
	} while (iter < max_iter);

	cout << "iterations: " << iter << "\n";

	return x;
}

Matrix Func(Function2 f1, Function2 f2, const Matrix& x) //vector-function
{
	Matrix result(2,1);
	result(0) = f1(x(0), x(1));
	result(1) = f2(x(0), x(1));
	return result;
}

Matrix Jacobi(Function2* J, Function2 f1, Function2 f2, const Matrix& x)
{
	Matrix result(2, 2);
	if (J != 0)	{
		result(0, 0) = J[0](x(0), x(1));
		result(0, 1) = J[1](x(0), x(1));
		result(1, 0) = J[2](x(0), x(1));
		result(1, 1) = J[3](x(0), x(1));
	}
	else {
		double const dx = 1e-04;
		result(0, 0) = (f1(x(0) + dx, x(1)) - f1(x(0), x(1))) / dx;
		result(0, 1) = (f1(x(0), x(1) + dx) - f1(x(0), x(1))) / dx;
		result(1, 0) = (f2(x(0) + dx, x(1)) - f2(x(0), x(1))) / dx;
		result(1, 1) = (f2(x(0), x(1) + dx) - f2(x(0), x(1))) / dx;
	}
	return result;
}

Matrix solve_newton(Function2 f1, Function2 f2, Function2* J, double l1, double l2, const Matrix& x0, double eps, int& iter)
{
	const int max_iter = 30;
	iter = 0;
	Matrix x(2,1);
	Matrix x_prev = x0;
	do {
		Matrix jacobi = Jacobi(J, f1, f2, x_prev);
		jacobi.invert2x2();

		//cout << "J:\n" << jacobi;
		//cout << "\n";

		x = x_prev - jacobi * Func(f1, f2, x_prev);
		++iter;

		//cout << iter << ": " << "\n" << x;
		//cout << "\n";

		if ((x - x_prev).norma_l2() <= eps) break;

		// avoid running out of area
		CheckAndFit(x, l1, l2);

		x_prev = x;
	} while (iter < max_iter);
	return x;
}


void convergence_area(Function2 f1, Function2 f2, Function2* J, double l1, double l2, int N, double eps, string filename)
{
	double step1 = 2 * l1 / N;
	double step2 = 2 * l2 / N;
	Matrix cells(N + 1, N + 1);
	Matrix X0(2, 1);
	int iter;
	for (int i = 0; i <= N; ++i) {
		for (int j = 0; j <= N; ++j) {
			X0(0) = -l1 + i * step1;
			X0(1) = -l2 + j * step2;
			solve_newton(f1, f2, J, l1, l2, X0, eps, iter);
			cells(i, j) = iter;
		}
	}
	cout << "convergence map:";
	cout << cells;
	cout << "\n";
	std::ofstream file(filename);
	cells.write_to_file(file);
	file.close();
}


const double eps = 1e-6;

double f_mytest(double x) { return (x - 0.06) * (x - 0.06) - 1; }
double df_mytest(double x) { return 2 * x - 0.12; }
Segment s_mytest = { -2, 2 };
// roots:  -0.94 1.06

double f_mytest2(double x) { return x * x * x * x; }
double df_mytest2(double x) { return 4 * x * x * x; }
Segment s_mytest2 = { -1, 2 };
// roots:  -0.94 1.06

double f_test1(double x) { return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75); }
double df_test1(double x) { return (x - 0.75) * (x - 0.7) * (x - 0.55) * (x - 0.22) + (x - 0.75) * (x - 0.7) * (x - 0.55) * (x - 0.1) + (x - 0.75) * (x - 0.7) * (x - 0.22) * (x - 0.1) + (x - 0.75) * (x - 0.55) * (x - 0.22) * (x - 0.1) + (x - 0.7) * (x - 0.55) * (x - 0.22) * (x - 0.1); }
Segment s_test1 = { 0, 1 };
// roos:  0.1  0.22  0.55  0.7  0.75

double f_test2(double x) { return sqrt(x + 1) - 1; }
double df_test2(double x) { return 1 / (2 * sqrt(x + 1)); }
Segment s_test2 = { -1, 10 };
// roots:  0

double f_test3(double x) { return 35 * x * x * x - 67 * x * x - 3 * x + 3; }
double df_test3(double x) { return 35 * 3 * x * x - 67 * 2 * x - 3; }
Segment s_test3 = { 0, 1 };
// roots:  0.2

//SLAY - my test
double f1_mytest(double x1, double x2) { return x1 + x2 - 3; }
double f2_mytest(double x1, double x2) { return x1 - x2 + 2; }
double L1_mytest = 10;
double L2_mytest = 10;
// roots:  {-3,1}  {1,-3}  {1, 2}  {2, 1}

double f1_test4(double x1, double x2) { return x1 * x1 - x2 * x2 - 15; }
double f2_test4(double x1, double x2) { return x1 * x2 + 4; }
double df1dx1_test4(double x1, double x2) { return 2 * x1; }
double df1dx2_test4(double x1, double x2) { return -2 * x2; }
double df2dx1_test4(double x1, double x2) { return x2; }
double df2dx2_test4(double x1, double x2) { return x1; }
Function2 j_test4[] = { df1dx1_test4, df1dx2_test4, df2dx1_test4, df2dx2_test4 };
double L1_test4 = 10;
double L2_test4 = 10;
// roots:  {4,-1}  {-4,1}

double f1_test5(double x1, double x2) { return x1 * x1 + x2 * x2 + x1 + x2 - 8; }
double f2_test5(double x1, double x2) { return x1 * x1 + x2 * x2 + x1 * x2 - 7; }
double df1dx1_test5(double x1, double x2) { return 2 * x1 + 1; }
double df1dx2_test5(double x1, double x2) { return 2 * x2 + 1; }
double df2dx1_test5(double x1, double x2) { return 2 * x1 + x2; }
double df2dx2_test5(double x1, double x2) { return 2 * x2 + x1; }
Function2 j_test5[] = { df1dx1_test5, df1dx2_test5, df2dx1_test5, df2dx2_test5 };
double L1_test5 = 10;
double L2_test5 = 10;
// roots:  {-3,1}  {1,-3}  {1, 2}  {2, 1}



int main() 
{
	cout << "Non-linear equation\n";
	cout << "-------------------\n";

	int iter;

	Function f = f_mytest2;
	Function df = df_mytest2;
	Segment segment = s_mytest2;

	vector<Segment> segments = locate_roots(f, segment, 0.03);
	if (segments.size() == 0) {
		segments.push_back(segment);
		cout << "The segments of roots cannot be located - the initial segment is taken:\n";
	}
	else {
		cout << "Located segments of roots:\n";
	}
	for (Segment& s : segments) {
		cout << "[" << s.a << "," << s.b << "]" << "\n";
	}

	cout << "\nBi-section method:\n";
	for (Segment& s : segments) {
		cout << solve_bisection(f, s, eps, iter) << "\n";
		cout << "iterations: " << iter << "\n";
	}

	cout << "\nNewton method:\n";
	for (Segment& s : segments) {
		double x0 = horde(f, s);
		cout << solve_newton(f, df, s, x0, eps, iter) << "\n";
		cout << "iterations: " << iter << "\n";
	}

	cout << "------------------------------\n\n";


	cout << "System of non-linear equations\n";
	cout << "------------------------------\n";

	Function2 f1 = f1_test4;
	Function2 f2 = f2_test4;
	Function2* j = 0; //j_test4;
	double L1 = L1_test4;
	double L2 = L2_test4;


	Matrix X0 = { {-2}, {1} };
	Matrix root = solve_newton(f1, f2, j, L1, L2, X0, eps, iter);
	cout << root;
	cout << "iterations: " << iter << "\n";


	cout << "------------------------------\n\n";

	cout << "Convergence area analysis\n";
	cout << "-------------------------\n";

	convergence_area(f1, f2, j, L1, L2, 10, eps, "convergence4.csv");

	return 0;
}