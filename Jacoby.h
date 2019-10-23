
#include "Simple_iter.h"

mytype	norm_oct(const VectorMatrix <mytype> &A)
{
	mytype	norm = 0;
	mytype	max = 0;

	for (int i = 0; i < A.count; i++)
	{
		norm = 0;
		for (int j = 0; j < A.count; j++)
				norm += fabs(A.value[j][i]);
		if (norm > max)
			max = norm;
	}
	return (max);
}

mytype	norm_cube(const VectorMatrix <mytype> &A)
{
	mytype	norm = 0;
	mytype	max = 0;

	for (int i = 0; i < A.count; i++)
	{
		norm = 0;
		for (int j = 0; j < A.count; j++)
				norm += fabs(A.value[i][j]);
		if (norm > max)
			max = norm;
	}
	return (max);
}

std::vector<mytype>	multi_vect(std::vector<mytype> I, const VectorMatrix<mytype> &T)
{
	std::vector<mytype> tmp(T.count);

	for (int i = 0; i < T.count; i++)
		tmp[i] = 0;
	for (int j = 0; j < T.count; j++)
		for (int i = 0; i < T.count; i++)
			tmp[j] += T.value[j][i] * I[i];
	return (tmp);
}

// should we check A[i][i] != 0 ????
std::vector<mytype>	jacoby(const VectorMatrix<mytype> &A, int flag)
{
	std::vector<mytype> x0(A.rvalue);
	std::vector<mytype>	x(A.rvalue);
	std::vector<mytype>	b(A.rvalue);
	std::vector<mytype> tmp(A.count);
	mytype	norm;
	VectorMatrix<mytype>	C(A.count);

	for (int i = 0; i < A.count; i++)
		b[i] /= A.value[i][i];
	for (int i = 0; i < A.count; i++)
		for (int j = 0; j < A.count; j++)
		{
			if (i != j)
				C.value[i][j] = -A.value[i][j] / A.value[i][i];
		}
	std::cout << "Matrix C in Jacoby method:" << std::endl;
	C.print();
	if (flag == 0)
	{
		norm = norm_cube(C);
		std::cout << "Matrix C cube norme = " << norm << std::endl;
		do
		{
			x0 = x;
			tmp = multi_vect(x0, C);
			x = sum_vect(tmp, b, x, A.count);
		} while (cube_vect_norm(diff_vector(x0, x, A.count), A.count) > ((1 - norm) / norm) * EPS);
	}
	else if (flag == 1)
	{
		norm = norm_oct(C);
		std::cout << "Matrix C octahedral norme = " << norm << std::endl;
		do
		{
			x0 = x;
			tmp = multi_vect(x0, C);
			x = sum_vect(tmp, b, x, A.count);
		} while (octah_vect_norm(diff_vector(x0, x, A.count), A.count) > ((1 - norm) / norm) * EPS);
	}
	return (x);
}
