#include "Linear_Regression.h"

#define N 10000
#define M 4
#define NOISE		20

#define ITENUM		500
#define BATCH			200
#define LEARNINGRATE	0.05
#define DECAY 0.9
#define DECAYNUM		80

int main(int argc, char **argv)
{
	mario::MyMatrix<double> X = mario::matrixRandom_Double(N, M, -10, 10);

	mario::MyMatrix<double> t = mario::matrixRandom_Double(M, 1, -100, 100);
	mario::MyMatrix<double> noise = mario::matrixRandom_Double(N, 1, -NOISE, NOISE);

	mario::MyMatrix<double> Y = X * t + noise;
	mario::MyMatrix<double> theta = mario::matrixRandom_Double(M, 1, -20, 20);


	clock_t start = clock();
	if (!LeastSquares(X, Y, theta))
	{
		std::cerr << "Error in main£ºLeastSquares return false.\n";
	}
	std::cout << "Time of leats squares£º" << clock() - start << "ms." << std::endl;
	std::cout << "Error of least squares£º\n" << theta - t;

	
	srand((unsigned)time(NULL));
	theta = mario::matrixRandom_Double(M, 1, -50, 50);
	start = clock();
	if (!GradientDescent(X, Y, theta, ITENUM, BGD, BATCH, LEARNINGRATE, DECAY, DECAYNUM))
	{
		std::cerr << "\nError in main£ºGradientDescent return false.\n";
	}
	std::cout << "Time of batch gradient descent£º" << clock() - start << "ms." << std::endl;
	std::cout << "Error of batch gradient descent£º\n" << theta - t;

	srand((unsigned)time(NULL));
	theta = mario::matrixRandom_Double(M, 1, -50, 50);
	start = clock();
	if (!GradientDescent(X, Y, theta, ITENUM, SGD, BATCH, LEARNINGRATE, DECAY, DECAYNUM))
	{
		std::cerr << "Error in main£ºGradientDescent return false.\n";
	}
	std::cout << "Time of stochastic gradient descent£º" << clock() - start << "ms." << std::endl;
	std::cout << "Error of stochastic gradient descent£º\n" << theta - t;

	return 0;
}
