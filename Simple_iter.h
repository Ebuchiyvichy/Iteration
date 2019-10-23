#pragma once
#include "MatrixClass.h"
//���������� ����� ������� � � ������ ������� ��������
mytype	norm_C_cube(const VectorMatrix <mytype> &A, mytype tau)
{
	mytype	norm = 0;
	mytype	max = 0;
	mytype	tmp;

	std::cout << "Matrix C in simple iteration:" << std::endl;
	for (int i = 0; i < A.count; i++)
	{
		norm = 0;
		for (int j = 0; j < A.count; j++)
		{
			if (i != j)
			{
				tmp = A.value[i][j] * tau;
				norm += fabs(tmp);
				std::cout << -tmp << '\t';
			}
			else
			{
				tmp = A.value[i][j] * tau - 1;
				norm += fabs(tmp);
				std::cout << -tmp << '\t';
			}
		}
		if (norm > max)
			max = norm;
		std::cout << std::endl;
	}
	std::cout << "Cube norm of matrix C = " << max << std::endl;
	return (max);
}
// //�������������� ����� ������� � � ������ ������� ��������
mytype	norm_C_oct(const VectorMatrix <mytype> &A, mytype tau)
{
	mytype	norm = 0;
	mytype	max = 0;
	mytype	tmp;

	std::cout << "Matrix C in simple iteration:" << std::endl;
	for (int i = 0; i < A.count; i++)
	{
		norm = 0;
		for (int j = 0; j < A.count; j++)
		{
			if (i != j)
			{
				tmp = A.value[j][i] * tau;
				norm += fabs(tmp);
				std::cout << tmp << '\t';
			}
			else
			{
				tmp = A.value[j][i] * tau - 1;
				norm += fabs(tmp);
				std::cout << tmp << '\t';
			}
		}
		if (norm > max)
			max = norm;
		std::cout << std::endl;
	}
	std::cout << "Octahedral norm of matrix C = " << max << std::endl;
	return (max);
}
//�������� ��������
std::vector<mytype> diff_vector(std::vector<mytype> a, std::vector<mytype> b, int DIM)
{
	std::vector<mytype> c(DIM);

	for (int i = 0; i < DIM; i++)
		c[i] = a[i] - b[i];
	return (c);
}

//���������� ����� �������
mytype	cube_vect_norm(std::vector<mytype> x, int size)
{
	mytype	norm = 0;

	for (int i = 0; i < size; i++)
		if (norm < fabs(x[i]))
			norm = fabs(x[i]);
	return (norm);
}
//�������������� ����� �������
mytype	octah_vect_norm(std::vector<mytype> x, int size)
{
	mytype	norm = 0;

	for (int i = 0; i < size; i++)
		norm += fabs(x[i]);
	return (norm);
}

// //�������� �������� � � b
std::vector<mytype> sum_vect(std::vector<mytype> a, std::vector<mytype> b, std::vector<mytype> c, int size)
{
	for (int i = 0; i < size; i++)
		c[i] = a[i] + b[i];
	return (c);
}


std::vector<mytype> cpy_vector(std::vector<mytype> x, int size, mytype tau)
{
	std::vector<mytype> tmp(size);

	for (int i = 0; i < size; i++)
		tmp[i] = x[i] * tau;
	return (tmp);
}

std::vector<mytype> multi_vect_matr_C(const VectorMatrix <mytype> &A, std::vector<mytype> x, mytype tau)
{
	std::vector<mytype> tmp(A.count);

	for (int j = 0; j < A.count; j++)
	{
		tmp[j] = 0;
		for (int i = 0; i < A.count; i++)
		{
			if (i != j)
				tmp[j] += A.value[j][i] * x[i] * tau;
			else
				tmp[j] += (A.value[j][i] * tau - 1) * x[i];

		}
	}
	return (tmp);
}

std::vector<mytype> diff_vect(std::vector<mytype> a, std::vector <mytype> b, std::vector<mytype> c, int size, mytype tau)
{
	for (int i = 0; i < size; i++)
		c[i] = b[i] * tau - a[i];
	return (c);
}

std::vector<mytype> simp_iter(const VectorMatrix <mytype> &A, int flag)
{
	std::vector<mytype> x(A.count);
	std::vector<mytype> x0(A.count);
	mytype	norm;
	std::vector<mytype> tmp(A.count);
	mytype	tau = 0.03;

	x = cpy_vector(A.rvalue, A.count, tau);
	if (flag == 0)
	{
		norm = norm_C_cube(A, tau);
		do
		{
			x0 = x;
			tmp = multi_vect_matr_C(A, x0, tau);
			x = diff_vect(tmp, A.rvalue, x, A.count, tau);
		} while (cube_vect_norm(diff_vector(x0, x, A.count), A.count) > ((1 - norm) / norm) * EPS);
	}
	else if (flag == 1)
	{
		tau *= EPS;
		norm = norm_C_oct(A, tau);
		do
		{
			x0 = x;
			tmp = multi_vect_matr_C(A, x0, tau);
			x = diff_vect(tmp, A.rvalue, x, A.count, tau);
		} while (octah_vect_norm(diff_vector(x0, x, A.count), A.count) > ((1 - norm) / norm) * EPS);
	}
	return (x);
}
