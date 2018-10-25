#include "MyMatrix.h"

int main(int argc, char **argv)
{
	mario::MyMatrix<double> a, b(3, 3, 1);
	//std::cout << b.get(0, 1);
	//a = b;
	//std::cout << a.rows() << a.cols() << a.get(1, 1);
	mario::MyMatrix<double> c(a);
	c.set(0, 0, 1);
	c.set(1, 1, 0);
	c.set(2, 2, 5);
	c.set(0, 2, 7);
	//std::cout << c.rows() << c.cols() << c.get(1, 1) << std::endl;
	std::cout << c << c * c << std::endl << c.inv_LU() * c << c.getDet() << std::endl;

	mario::MyMatrix<double> d = mario::matrixRandom_Double(5, 5, 0., 1.);

	std::cout << d << d.inv_LU() * d;

	d = mario::diag(5, 2.1);

	std::cout << d;

	return 0;
}