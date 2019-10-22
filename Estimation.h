//#pragma once
#include "Condition.h"
#include <vector>

//find octahedral norme of vector
double	octah_vect_norm(std::vector <double> x, int size)
{
	double	norm = 0;

	for (int i = 0; i < size; i++)
		norm += fabs(x[i]);
	return (norm);
}

//find cube norme of vector
double	cube_vect_norm(std::vector <double> x, int size)
{
	double	norm = 0;

	for (int i = 0; i < size; i++)
		if (norm < x[i])
			norm = x[i];
	return (norm);
}

//find sferical norm of vector
double	sfer_vect_norm(std::vector <double> x, int DIM)
{
	double	norm = 0;

	for (int i = 0; i < DIM; i++)
		norm += x[i] * x[i];
	return (sqrt(norm));
}

//find absolutely difference of two vectors
double	*abs_diff_vector(std::vector <double> a, std::vector <double> b, int DIM)
{
	double	*c = new double[DIM];

	for (int i = 0; i < DIM; i++)
		c[i] = fabs(a[i] - b[i]);
	return (c);
}

// //find difference of two vectors
// double	*diff_vector(double *a, double *b, int DIM)
// {
// 	double	*c = new double[DIM];
//
// 	for (int i = 0; i < DIM; i++)
// 		c[i] = a[i] - b[i];
// 	return (c);
// }

//find conditional number's estimate

//for cube norme
double	cube_estimate_number(std::vector<double> b, const VectorMatrix <double> &R, const VectorMatrix <double> &T)
{
	int DIM = T.count;
  std::vector<double> b_(DIM);
  std::vector<double> tmp(DIM);
	std::vector<double> x(DIM);
  std::vector<double> x_(DIM);
	std::vector<double>dx(DIM);
	double	my_dx = 0;
	std::vector<double> db(DIM) ;

	b_ = multi_vect(b, T);
	x = find_x(R, b_);
	for (int i = 0; i < 2 * T.count; i++)
	{
		tmp = b;
		if (i < T.count)
			tmp[i] += pert;
		else
			tmp[i - T.count] -= pert;
		b_ = multi_vect(tmp, T);
		x_ = find_x(R, b_);
		dx = abs_diff_vector(x, x_, T.count);
		if (my_dx < cube_vect_norm(dx, T.count))
			my_dx = cube_vect_norm(dx, T.count);
	}
	db = abs_diff_vector(tmp, b, T.count);
	return ((my_dx / cube_vect_norm(x, T.count)) / (cube_vect_norm(db, T.count) / cube_vect_norm(b, T.count)));
}

//for octahedral norme
double	octah_estimate_number(std::vector<double> b, const VectorMatrix <double> &R, const VectorMatrix<double> &T)
{
	double	pert = 0.1;
	int DIM = T.count;
	std::vector<double> b_(DIM);
  std::vector<double> tmp (DIM);
	std::vector<double> x(DIM);
  std::vector<double> x_(DIM);
	std::vector<double> dx(DIM);
	double	my_dx = 0;
	std::vector<double> db(DIM);

	b_ = multi_vect(b, T);
	x = find_x(R, b_);
	for (int i = 0; i < 2 * DIM; i++)
	{
		tmp = b;
		if (i < DIM)
			tmp[i] += pert;
		else
			tmp[i - DIM] -= pert;
		b_ = multi_vect(tmp, T);
		x_ = find_x(R, b_);
		dx = abs_diff_vector(x, x_, DIM);
		if (my_dx < octah_vect_norm(dx, DIM))
			my_dx = octah_vect_norm(dx, DIM);
	}
	db = abs_diff_vector(tmp, b, DIM);
	return ((my_dx / octah_vect_norm(x, DIM)) / (octah_vect_norm(db, DIM) / octah_vect_norm(b, DIM)));
}

//for sferical norm
double	sfer_estimate_number(std::vector<double> b, const   VectorMatrix <double>& R, const   VectorMatrix <double>& T)
{
	double	pert = 0.1;
	int DIM = T.count;
  std::vector<double> b_(DIM);
  std::vector<double> tmp (DIM);
	std::vector<double> x(DIM);
  std::vector<double> x_(DIM);
	std::vector<double> dx(DIM);
	double	my_dx = 0;
	std::vector<double> db(DIM);

	b_ = multi_vect(b, T);
	x = find_x(R, b_);
	for (int i = 0; i < 2 * DIM; i++)
	{
		tmp = tmp;
		if (i < DIM)
			tmp[i] += pert;
		else
			tmp[i - DIM] -= pert;
		b_ = multi_vect(tmp, T);
		x_ = find_x(R, b_);
		dx = abs_diff_vector(x, x_, DIM);
		if (my_dx < sfer_vect_norm(dx, DIM))
			my_dx = sfer_vect_norm(dx, DIM);
	}
	db = abs_diff_vector(tmp, b, DIM);
	return ((my_dx / sfer_vect_norm(x, DIM)) / (sfer_vect_norm(db, DIM) / sfer_vect_norm(b, DIM)));
}
