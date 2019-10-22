#include <cmath>
#include <vector>
#include "NewMatrixClass.h"

// double *cpy_vector(double *tmp, double *x, int size)
// {
//
// 	for (int i = 0; i < size; i++)
// 		tmp[i] = x[i];
// 	return (tmp);
// }

bool	degenerate_matrix(const VectorMatrix <double> & R)
{
	for (int i = 0; i < R.size; i++)
		if (fabs(R.value[i][i]) < EPS)
			return (true);
	return (false);
}

// reverse Hause
void find_x(const VectorMatrix <double> & A, std::vector <double> b_, std::vector <double> &x)
{
	for (int i = A.count - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < A.count; j++)
			b_[i] = b_[i] - A.value[i][j] * b_[j];
		b_[i] = b_[i] / A.value[i][i];
	}
	x = b_;
}

// find matrix T
void VectorMatrix::QR_find_x(VectorMatrix <double> & A)
{
	T	c;
	T	s;
	T	tmp;

	for (int k = 0; k < count; k++)
	{
		for (int i = k + 1; i < count; i++)
		{
			if (fabs(A.value[k][i]) > 10e-8)
			{
				c = A.value[k][k] / sqrt(A.value[k][k] * A.value[k][k] + A.value[i][k] * A.value[i][k]);
				s = A.value[i][k] / sqrt(A.value[k][k] * A.value[k][k] + A.value[i][k] * A.value[i][k]);
				for (int j = 0; j <= i; j++) // change T-matrix
				{
					tmp = value[k][j];
					value[k][j] = value[k][j] * c + value[i][j] * s;
					value[i][j] = c * value[i][j] - s * tmp;
				}
				for (int j = k; j < size; j++) // change A-matrix
				{
					tmp = A.value[k][j];
					A.value[k][j] = c * A.value[k][j] + s * A.value[i][j];
					A.value[i][j] = c * A.value[i][j] - s * tmp;
				}
			}
		}
	}
	std::cout << "QR method matrix R:" << std::endl;
	A.print();

	std::cout << "QR method matrix T:" << std::endl;
	print();
}
