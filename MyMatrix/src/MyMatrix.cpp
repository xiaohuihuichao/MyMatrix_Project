/*
#include "MyMatrix.h"

namespace mario
{
	template<class T>
	MyMatrix<T>::MyMatrix(const int &row, const int &col, const T &data)
	{
		if (row <= 0 || row <= 0)
		{
			std::cerr << "Error in MyMatrix<T>::MyMatrix(const int row, const int col, const T &data)：rows and columns must be > 0.\n";
		}

		iRow = row;
		iCol = col;
		pData = new T[iRow * iCol];
		for (int i = 0; i < iRow*iCol; ++i)
		{
			*(pData + i) = data;
		}
	}


	template<class T>
	MyMatrix<T>::MyMatrix(const MyMatrix<T> &m)
	{
		pData = new T[m.iRow * m.iCol];
		iRow = m.iRow;
		iCol = m.iCol;
		memcpy((void*)pData, (void*)m.pData, iRow*iCol*sizeof(T));
	}


	template<class T>
	T MyMatrix<T>::get(const int &row, const int &col) const
	{
		if (row < 0 || row >= iRow
			|| col < 0 || col >= iCol)
		{
			std::cerr << "Error in MyMatrix<T>::get(int row, int col) const：the index is out of data.(row, col)=(" << row << ", " << col << ")\n";
		}

		return pData[row*iCol + col];
	}

	template<class T>
	bool MyMatrix<T>::set(const int &row, const int &col, const T&data)
	{
		if (row < 0 || row >= iRow
			|| col < 0 || col >= iCol)
		{
			std::cerr << "the index is out of data.(row, col)=(" << row << ", " << col << ")\n";
			return false;
		}

		*(pData + row*iCol + col) = data;
		return true;
	}

	template<class T>
	bool MyMatrix<T>::swapRow(const int &i, const int &j)
	{
		if (i < 0 || j < 0
			|| i >= iRow || j >= iCol)
		{
			std::cerr << "Error in bool MyMatrix<T>::swap(int, int)：the index is out of boundary.\n";
			return false;
		}

		T temp = 0;
		T *pI = pData + i*iCol;
		T *pJ = pData + j*iCol;
		for (int k = 0; k < iCol; ++k)
		{
			temp = *(pI + k);
			*(pI + k) = *(pJ + k);
			*(pJ + k) = temp;
		}

		return true;
	}


	template<class T>
	T MyMatrix<T>::getDet() const
	{
		if (!isSquareMatrix())
		{
			std::cerr << "Error in double MyMatrix<T>::getDet() const：the matrix is not a square matrix.\n";
			return 0;
		}

		MyMatrix<T> copyMatrix(*this);
		//记录行变换次数
		int iter = 0;

		//行列式结果
		T det = 1;
		for (int i = 0; i < iRow; ++i)
		{
			if (0 == copyMatrix(i, i))
			{
				for (int j = i; j < iRow; ++j)
				{
					if (copyMatrix(j, i) != 0)
					{
						if (!copyMatrix.swapRow(i, j))
						{
							return 0;
						}
						++iter;
					}
				}
			}

			for (int k = i + 1; k < iRow; ++k)
			{
				double scale = -1 * copyMatrix(k, i) / copyMatrix(i, i);
				double data = 0;
				for (int u = 0; u < iRow; ++u)
				{
					data = copyMatrix(k, u) + copyMatrix(i, u)*scale;
					copyMatrix.set(k, u, data);
				}
			}
		}

		for (int i = 0; i < iRow; ++i)
		{
			det *= copyMatrix(i, i);
		}

		if (1 == iter % 2)
		{
			det *= -1;
		}

		return det;

	}


	template<class T>
	MyMatrix<double> MyMatrix<T>::inv_LU() const
	{
		if (!isSquareMatrix())
		{
			std::cerr << "Error in MyMatrix<double> MyMatrix<T>::inv)LU() const：this is not a square matrix.(row, col)=(" << iRow << ", " << iCol << ")\n";
		}

		double det = getDet();
		if (det >= -EPSINON && det <= EPSINON)
		{
			std::cerr << "Error in MyMatrix<double> MyMatrix<T>::inv_LU() const：the determinant is equal to 0.\n";
		}

		int n = iRow;
		//建立l、l_inverse、u、u_inverse矩阵
		MyMatrix<double> l(n, n);
		MyMatrix<double> l_inverse(n, n);
		MyMatrix<double> u(n, n);
		MyMatrix<double> u_inverse(n, n);

		//计算矩阵对角线
		for (int i = 0; i < n; ++i)
		{
			l.set(i, i, 1);
			l_inverse.set(i, i, 1);
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = i; j < n; ++j)
			{
				double s = 0;
				for (int k = 0; k < i; ++k)
				{
					s += l(i, k) * u(k, j);	//按行计算u值
				}
				u.set(i, j, get(i, j) - s);
			}

			for (int j = i + 1; j < n; ++j)
			{
				double s = 0;
				for (int k = 0; k < i; ++k)
				{
					s += l(j, k) * u(k, i);
				}
				l.set(j, i, (get(j, i) - s) / u(i, i));	//按列计算l值
			}
		}

		for (int i = 1; i < n; ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				double s = 0;
				for (int k = 0; k < i; ++k)
				{
					s += l(i, k) * l_inverse(k, j);
				}
				l_inverse.set(i, j, -s);
			}
		}

		for (int i = 0; i < n; ++i)
		{
			u_inverse.set(i, i, 1 / u(i, i));
		}
		for (int i = 1; i < n; ++i)
		{
			for (int j = i - 1; j >= 0; --j)
			{
				double s = 0;
				for (int k = j + 1; k <= i; ++k)
				{
					s += u(j, k) * u_inverse(k, i);
				}
				u_inverse.set(j, i, -s / u(j, j));
			}
		}

		MyMatrix<double> result(u_inverse * l_inverse);

		for (int i = 0; i < result.rows(); ++i)
		{
			for (int j = 0; j < result.cols(); ++j)
			{
				if (_isnan(result.get(i, j)))
				{
					std::cerr << "Error in MyMatrix<double> MyMatrix<T>::inv_LU() const：the result is 1.#INF or -1.#IND or 1.#INF000 or -1.#INF000.\n";
				}
			}
		}

		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::t() const
	{
		MyMatrix<T> result(iCol, iRow);
		T data;
		for (int i = 0; i < iCol; ++i)
		{
			for (int j = 0; j < iRow; ++j)
			{
				data = get(j, i);
				result.set(i, j, data);
			}
		}
		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::mul(const MyMatrix<T> &m) const
	{
		if (iRow != m.iRow || iCol != m.iCol)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::mul(const MyMatrix<T>&) const：the index of two matrix is not equal.\n";
		}
		if (0 == iRow || 0 == iCol)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::mul(const MyMatrix<T>&) const：the row or col of matrix is equal to 0.\n";
		}

		T e;
		MyMatrix result(iRow, iCol);
		for (int i = 0; i < iRow*iCol; ++i)
		{
			e = *(m.pData + i);
			*(result.pData + i) = *(pData + i) * e;
		}

		return result;
	}


	template<class T>
	T mario::MyMatrix<T>::operator()(const int &row, const int &col) const
	{
		if (row < 0 || row >= iRow
			|| col < 0 || col >= iCol)
		{
			std::cerr << "Error in MyMatrix<T>::operator()(int row, int col) const：the index is out of data.(row, col)=(" << row << ", " << col << ")\n";
		}
		return get(row, col);
	}


	//运算符=重载
	template<class T>
	MyMatrix<T>& MyMatrix<T>::operator=(const MyMatrix<T> &m)
	{
		//can not copy itself
		if (this == &m)
		{
			return *this;
		}
		if (NULL != pData
			&& (iRow != m.iRow || iCol != m.iCol))
		{
			std::cerr << "Error in MyMatrix<T>& MyMatrix<T>::operator=(const MyMatrix<T>&)：the size of two matrix is not equal.\n";
		}

		pData = new T[m.iRow*m.iCol];
		iRow = m.iRow;
		iCol = m.iCol;
		memcpy((void*)pData, (void*)m.pData, iRow*iCol*sizeof(T));
		return *this;
	}


	template<class T>
	bool MyMatrix<T>::operator==(const MyMatrix<T> &m) const
	{
		if (iRow != m.iRow || iCol != m.iCol)
		{
			return false;
		}

		for (int i = 0; i < iRow*iCol; ++i)
		{
			if (*(pData + i) != *(m.pData + i))
			{
				return false;
			}
		}

		return true;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator+(const MyMatrix<T> &m) const
	{
		if (iRow != m.iRow || iCol != m.iCol)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::operator+(const MyMatrix<T>&) const：the size of two matrix is not equal, so can not +.\n";
		}

		T e;
		MyMatrix<T> result(iRow, iCol);
		for (int i = 0; i < iRow*iCol; ++i)
		{
			e = *(m.pData + i);
			*(result.pData + i) = *(pData + i) + e;
		}

		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator+(const T &e) const
	{
		MyMatrix<T> result(iRow, iCol);
		for (int i = 0; i < iRow*iCol; ++i)
		{
			*(result.pData + i) = *(pData + i) + e;
		}

		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator-(const MyMatrix<T> &m) const
	{
		if (iRow != m.iRow || iCol != m.iCol)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::operator-(const MyMatrix<T>&) const：the size of two matrix is not equal, so can not -.\n";
		}

		T e;
		MyMatrix<T> result(iRow, iCol);
		for (int i = 0; i < iRow*iCol; ++i)
		{
			e = *(m.pData + i);
			*(result.pData + i) = *(pData + i) - e;
		}

		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator-(const T& e) const
	{
		MyMatrix<T> result(iRow, iCol);
		for (int i = 0; i < iRow*iCol; ++i)
		{
			*(result.pData + i) = *(pData + i) - e;
		}

		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator*(const MyMatrix<T> &m) const
	{
		if (iCol != m.cols())
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::operator*(const MyMatrix<T>&) const：the size of two matrix is not match.\n";
		}

		MyMatrix<T> result(iRow, m.iCol);
		int index1 = 0;
		int index2 = 0;
		int indexResult = 0;
		for (int i = 0; i < result.iRow; ++i)
		{
			for (int j = 0; j < result.iCol; ++j)
			{
				indexResult = i*result.iCol + j;
				*(result.pData + indexResult) = 0;
				for (int k = 0; k < iCol; ++k)
				{
					index1 = i*iCol + k;
					index2 = k*m.iCol + j;
					*(result.pData + indexResult) += *(pData + index1) * (*(m.pData + index2));
				}
			}
		}
		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator*(const T &e) const
	{
		MyMatrix<T> result(iRow, iCol);

		for (int i = 0; i < iRow*iCol; ++i)
		{
			*(result.pData + i) = *(pData + i) * e;
		}

		return result;
	}


	template<class T>
	MyMatrix<double> MyMatrix<T>::operator/(const MyMatrix<T> &m) const
	{
		if (!m.isSquareMatrix() || m.iRow != iCol)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::operator/(const MyMatrix<T>&)：this is not a square matrix, or the size of two matrix is not match.\n";
		}
		MyMatrix<double> m_inv = m.inv_LU();
		return (*this) * m_inv;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::operator/(const T &e) const
	{
		if (e >= -EPSINON && e <= EPSINON)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::operator/(const T&) const：the element is equal to 0.\n";
		}

		MyMatrix<T> result(iRow, iCol);

		for (int i = 0; i < iRow*iCol; ++i)
		{
			*(result.pData + i) = *(pData + i) / e;
		}

		return result;
	}


	template<class T>
	MyMatrix<T> MyMatrix<T>::block(const int &row, const int &col, const int &height, const int &width) const
	{
		if (row < 0 || col < 0 || height <= 0 || width <= 0)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::block(int, int, int, int) const：the row, col must >= 0, and height and width must > 0.\n";
		}
		if (row + height > iRow || col + width > iCol)
		{
			std::cerr << "Error in MyMatrix<T> MyMatrix<T>::block(int, int, int, int) const：the block is out of boundary.\n";
		}

		MyMatrix<T> result(height, width);
		for (int i = 0; i < height; ++i)
		{
			for (int j = 0; j < width; ++j)
			{
				*(result.pData + i*width + j) = *(pData + (row + i)*iCol + col + j);
			}
		}

		return result;
	}


	template<class T>
	std::ostream& operator<<(std::ostream& out, const MyMatrix<T> &m)
	{
		out << "[\n";
		for (int i = 0; i < m.rows(); ++i)
		{
			out << "  ";
			for (int j = 0; j < m.cols(); ++j)
			{
				out << m.get(i, j) << ",\t";
			}
			out << "\n";
		}
		out << "]\n";
		return out;
	}

	template<class T>
	MyMatrix<T> diag(const int &n, const T &data)
	{
		if (n <= 0)
		{
			std::cerr << "Error in MyMatrix<T> diag(const int&, const int&)：the size of matrix must > 0.\n";
		}
		
		MyMatrix<T> m(n, n, 0);
		for (int i = 0; i < n; ++i)
		{
			m.set(i, i, data);
		}
		return m;
	}


	MyMatrix<int> matrixRandom(const int &row, const int &col, const int &start, const int &end)
	{
		if (row <= 0 || col <= 0)
		{
			std::cerr << "Error in MyMatrix<T> matrixRandom(const int&, const int&, const int&, const int&)：the size of matrix must > 0.\n";
		}
		if (start >= end)
		{
			std::cerr << "Error in MyMatrix<T> matrixRandom(const int&, const int&, const int&, const int&)： \
													the end number must > start number.\n";
		}

		MyMatrix<int> m(row, col);
		int data = 0;
		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
			{
				srand((unsigned)time(NULL));
				data = rand() % (end - start) + end + 1;
				m.set(i, j, data);
			}
		}
		return m;
	}

	
	MyMatrix<double> matrixRandom(const int &row, const int &col, const double &start, const double &end)
	{
		if (row <= 0 || col <= 0)
		{
			std::cerr << "Error in MyMatrix<T> matrixRandom(const int&, const int&, const int&, const int&)：the size of matrix must > 0.\n";
		}
		if (start >= end)
		{
			std::cerr << "Error in MyMatrix<T> matrixRandom(const int&, const int&, const int&, const int&)： \
						 													the end number must > start number.\n";
		}

		MyMatrix<double> m(row, col);
		double data = 0;
		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
			{
				srand((unsigned)time(NULL));
				data = rand() / double(RAND_MAX) *(end - start) + start;
				m.set(i, j, data);
			}
		}
		return m;
	}

}//end of mario
*/