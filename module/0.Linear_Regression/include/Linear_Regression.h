#ifndef LINEAR_REGRESSION_H_
#define LINEAR_REGRESSION_H_

#include "MyMatrix.h"


template<class T>
bool LeastSquares(const mario::MyMatrix<T> &X, const mario::MyMatrix<T> &Y, mario::MyMatrix<T> &theta)
{
	if (X.rows() != Y.rows() || X.cols() != theta.rows() || theta.cols() != Y.cols())
	{
		std::cerr << "Error in bool LeastSquares(const mario::MyMatrix<T>&, const mario::MyMatrix<T>&, mario::MyMatrix<T>&)£º\
														the size of array is not match.\n";
		return false;
	}

	mario::MyMatrix<T> XTX = X.t() * X;

	if (XTX.getDet() <= mario::EPSINON && XTX.getDet() >= -1.*mario::EPSINON)
	{
		std::cerr << "Error in LeastSquares(const mario::MyMatrix<T>&, const mario::MyMatrix<T>&, mario::MyMatrix<T>&)£º\
					 														the input is linearly dependent.\n";
		return false;
	}

	//		(X.T * X).inv * (X.T * Y)
	mario::MyMatrix<T> XTX_inv = XTX.inv_LU();
	mario::MyMatrix<T> XTY = X.t() * Y;
	theta = XTX_inv * XTY;

	return true;
}


#define BGD			0
#define SGD			1
template<class T>
bool GradientDescent(const mario::MyMatrix<T> &X, const mario::MyMatrix<T> &Y, mario::MyMatrix<T> &theta
	, int iter, int mode = BGD, int batch = 1,double learningRate = 0.1, double decay = 1., int decayNum = 100)
{
	if (X.cols() != theta.rows() || X.rows() != Y.rows() || theta.cols() != Y.cols())
	{
		std::cerr << "Error in bool GradientDescent(const mario::MyMatrix<T>&, const mario::MyMatrix<T>&, mario::MyMatrix<T>&\
								, int, int, double, double, int, int)£ºthe size of array is not match.\n";
		return false;
	}
	if (batch > X.rows() || batch <=0)
	{
		std::cerr << "Error in bool GradientDescent(const mario::MyMatrix<T>&, const mario::MyMatrix<T>&, mario::MyMatrix<T>&\
					 								, int, int, double, double, int, int)£ºthe batch must < number of samples and must > 0.\n";
		return false;
	}

	// GD = X.T * X * theta- X.T * Y
	mario::MyMatrix<T> deltaTheta(theta.rows(), theta.cols());
	int num = 0;
	switch (mode)
	{
	case BGD:
		for (int i = 0; i < iter; ++i)
		{
			if (num >= decayNum)
			{
				learningRate *= decay;
				num = 0;
			}

			//X.T * X * theta - X.T * Y
			mario::MyMatrix<T> y = X * theta;
			mario::MyMatrix<T> deltaTheta = (X.t() * y - X.t() * Y) / X.rows();
			theta -= deltaTheta * learningRate;
			++num;
		}
		break;

	case SGD:
		for (int i = 0; i < iter; ++i)
		{
			
			for (int j = 0; j < X.rows() / batch; ++j)
			{
				if (num >= decayNum)
				{
					learningRate *= decay;
					num = 0;
				}

				int startIndex = i*batch%X.rows();
				int endIndex = (startIndex + batch) < X.rows() ? (startIndex + batch) : X.rows();
				mario::MyMatrix<T> x = X.block(startIndex, 0, endIndex - startIndex, X.cols());
				mario::MyMatrix<T> y_ = Y.block(startIndex, 0, endIndex - startIndex, Y.cols());

				//X.T * X * theta - X.T * Y_
				mario::MyMatrix<T> y = x * theta;
				deltaTheta = (x.t() * y - x.t() * y_) / x.rows();
				theta -= deltaTheta * learningRate;
				++num;
			}
		}
		break;

	default:
		std::cerr << "Error in bool GradientDescent(const mario::MyMatrix<T>&, const mario::MyMatrix<T>&, mario::MyMatrix<T>&\
					 								, int, int, double, double, int, int)£ºthe mode is wrong.\n";
		return false;
	}

	return true;
}


#endif