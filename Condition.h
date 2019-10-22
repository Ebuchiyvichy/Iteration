//#pragma once
#include "QR.h"


//multiply matrix * vector
double *multi_vect(std::vector<double> I, const VectorMatrix  <double> & T)
{
	double *tmp;

	tmp = new double[T.count];
	for (int i = 0; i < T.count; i++)
		tmp[i] = 0;
	for (int j = 0; j < T.count; j++)
		for (int i = 0; i < T.count; i++)
			tmp[j] += T.value[j][i] * I[i];
	return (tmp);
}

//find inverse matrix
void Matrix::inverse_matrix(const VectorMatrix  <double>& R, VectorMatrix <double>& T)
{
/*
	//with Hausse
	double *x;
	double *y;
	x = new double [T.count];
	y = new double [T.count];
	Matrix E(R.count, 1);

	for (int i = 0; i < T.count; i++)
	{
	for (int j = 0; j < T.count; j++)
	{
	if (j == i)
	y[j] = 1;
	else
	y[j] = 0;
	}
	x = find_x(R, multi_vect(y, T));
	for (int j = 0; j < T.count; j++)
	{
	value[mymap[i] - value + j] = x[j];
	}
	}
	Tranc();
*/

	//with inverse triangle matrix
	  VectorMatrix <double> <double> E(R.count);
	E.onebyone();
	for (int k = R.count - 1; k >= 0; k--)
	{
		E.value[k][k] = E.value[k][k] / R.value[k][k];
		for (int j = k - 1; j >= 0; j--)
			{
				for (int i = R.count - 1; i > j; i--)
					E.value[j][k] -= R.value[j][i] * E.value[i][k];
				E.value[j][k] /= R.value[j][j];
			}
	}
/*	for (int k = 0; k < R.count; k++)
	{
		value[mymap[k] - value + k] = value[mymap[k] - value + k] / R.value[R.mymap[k] - R.value + k];
		for (int j = k - 1; j >= 0; j--)
		{
			for (int i = R.count - 1; i > j; i--)
				value[mymap[j] - value + k] -= R.value[R.mymap[j] - R.value + i] * value[mymap[i] - value + k];
			value[mymap[j] - value + k] /= R.value[R.mymap[j] - R.value + j];
		}
	}
	*/
	for (int i = 0; i < R.count; i++)
	{
		for (int j = 0; j < R.count; j++)
		{
			value[i][j] = 0;
			for (int k = 0; k < R.count; k++)
				value[i][j] = E.value[i][k] * T.value[k][j];
		}
	}
}
//
// //find cube norm
// double	cube_norm(const Matrix& A)
// {
// 	double	norm = 0;
// 	double	sum = 0;
//
// 	for (int i = 0; i < A.count; i++)
// 	{
// 		for (int j = 0; j < A.count; j++)
// 			sum += fabs(A.value[A.mymap[i] - A.value + j]);
// 		if (sum > norm)
// 			norm = sum;
// 		sum = 0;
// 	}
// 	return (norm);
// }
//
// //find octahedral norm
// double	octah_norm(const Matrix& A)
// {
// 	T	norm = 0;
// 	T	sum = 0;
//
// 	for (int i = 0; i < A.count; i++)
// 	{
// 		for (int j = 0; j < A.count; j++)
// 			sum += fabs(A.value[A.mymap[j] - A.value + i]);
// 		if (sum > norm)
// 			norm = sum;
// 		sum = 0;
// 	}
// 	return (norm);
// }
//
// //find sferical norm
// T	sfer_norm(const Matrix& A)
// {
// 	T	sum = 0;
//
// 	for (int i = 0; i < A.count; i++)
// 	{
// 		for (int j = 0; j < A.count; j++)
// 			sum += fabs(A.value[A.mymap[i] - A.value + j]);
// 	}
// 	return (sqrt(sum));
// }
