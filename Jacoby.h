
#include "Simple_iter.h"

// mytype	*cpy_vector(std::vector<mytype> x, int size)
// {
// 	mytype	*tmp = new mytype[size];
//
// 	for (int i = 0; i < size; i++)
// 		tmp[i] = x[i];
// 	return (tmp);
// }

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

mytype	*multi_vect(mytype* I, const VectorMatrix<mytype> &T)
{
	mytype *tmp;

	tmp = new mytype[T.count];
	for (int i = 0; i < T.count; i++)
		tmp[i] = 0;
	for (int j = 0; j < T.count; j++)
		for (int i = 0; i < T.count; i++)
			tmp[j] += T.value[j][i] * I[i];
	return (tmp);
}

std::vector<mytype> multi_vect(std::vector<mytype> I, const VectorMatrix<mytype> &T)
{
	std::vector<mytype> tmp(T.count, 0);

	for (int j = 0; j < T.count; j++)
		for (int i = 0; i < T.count; i++)
			tmp[j] += T.value[j][i] * I[i];
	return (tmp);
}
//should we check A[i][i] != 0 ????
// mytype	*jacoby(const VectorMatrix<mytype> &A, int flag)
// {
// 	mytype	*x0 = new mytype[A.count];
// 	mytype	*x = new mytype[A.count];
// 	mytype	*b, *tmp;
// 	mytype	norm;
// 	VectorMatrix<mytype>	C(A.count);
//
// 	b = cpy_vector(A.rvalue, A.count);
// 	for (int i = 0; i < A.count; i++)
// 		b[i] /= A.value[i][i];
// 	x = cpy_vector(x, b, A.count);
// 	for (int i = 0; i < A.count; i++)
// 		for (int j = 0; j < A.count; j++)
// 		{
// 			if (i != j)
// 				C.value[i][j] = -A.value[i][j] / A.value[i][i];
// 		}
// 	std::cout << "Matrix C in Jacoby method:" << std::endl;
// 	C.print();
// 	if (flag == 0)
// 	{
// 		norm = norm_cube(C);
// 		std::cout << "Matrix C cube norme = " << norm << std::endl;
// 		do
// 		{
// 			x0 = cpy_vector(x0, x, A.count);
// 			tmp = multi_vect(x0, C);
// 			x = sum_vect(tmp, b, x, A.count);
// 			delete[] tmp;
// 		} while (cube_vect_norm(diff_vector(x0, x, A.count), A.count) > ((1 - norm) / norm) * EPS);
// 	}
// 	else if (flag == 1)
// 	{
// 		norm = norm_oct(C);
// 		std::cout << "Matrix C octahedral norme = " << norm << std::endl;
// 		do
// 		{
// 			x0 = cpy_vector(x0, x, A.count);
// 			tmp = multi_vect(x0, C);
// 			x = sum_vect(tmp, b, x, A.count);
// 			delete[] tmp;
// 		} while (octah_vect_norm(diff_vector(x0, x, A.count), A.count) > ((1 - norm) / norm) * EPS);
// 	}
// 	delete[] x0;
// 	delete[] b;
// 	return (x);
// }
