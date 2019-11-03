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

std::vector<mytype>	jacoby(const VectorMatrix<mytype> &A, int flag)
{
	std::vector<mytype>		x0(A.rvalue);
	std::vector<mytype>		x(A.rvalue);
	std::vector<mytype>		b(A.rvalue);
	std::vector<mytype>		tmp(A.count);
	mytype					norm;
	VectorMatrix<mytype>	C(A.count);
	int						iter = 0;

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
			iter++;
			x0 = x;
			tmp = multi_vect(x0, C);
			x = sum_vect(tmp, b, x, A.count);
			std::cout << "Residual vector with cube norme at " << iter << " step = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
		} while (cube_vect_norm(diff_vector(x0, x, A.count), A.count) > (fabs(1 - norm) / norm) * EPS);
	//	} while (cube_vect_norm(diff_vector(x0, x, A.count), A.count) > EPS);
	//	} while (cube_vect_norm(diff_vector(x0, x, A.count), A.count) > 10 * EPS);
	//	} while (cube_vect_norm(diff_vector(multi_vect(x, A), A.rvalue, A.count), A.count) > EPS);
	}
	else if (flag == 1)
	{
		norm = norm_oct(C);
		std::cout << "Matrix C octahedral norme = " << norm << std::endl;
		do
		{
			iter++;
			x0 = x;
			tmp = multi_vect(x0, C);
			x = sum_vect(tmp, b, x, A.count);
			std::cout << "Residual vector with octahedral norme at " << iter << " step = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
		} while (octah_vect_norm(diff_vector(x0, x, A.count), A.count) > (fabs(1 - norm) / norm) * EPS);
	}
	std::cout << "Number of iteration: " << iter << std::endl;
	return (x);
}
