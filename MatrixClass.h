#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
typedef double mytype;

mytype	EPS = 10e-3;

template <typename T> class VectorMatrix
{
public:
	int count;
	std::vector<T> rvalue;
	std::vector<std::vector <T>> value;

public:
	VectorMatrix(int size)
	{
		count = size;
		for (int i = 0; i != count; i++)
		{
			value.push_back(std::vector<T>(count));
		}
	}

	~VectorMatrix()
	{
	}

	void VectorMatrixC(int size)
	{
		count = size;
		for (int i = 0; i != 3; i++)
		{
			value.push_back(std::vector<mytype>(count));
		}
		for (int i = 0; i != count; i++)
		{
			for (int j = 0; j != 3; j++)
			{
				if (j == 0) {
					if (i == 0)
						value[j][i] = 0;
					else
						value[j][i] = 1;
				}
				else if (j == 1)
						value[j][i] = 4;
				else if (j == 2) {
					if (i == 201)
						value[j][i] = 0;
					else
						value[j][i] = 1;
				}
			}
		}
		std::cout << "I'm really works fine\n";
		for (int i = 0; i != count; i++)
		{
			//rvalue.push_back();
			if (i == 0)
				rvalue.push_back(6);
			else if (i == 201)
				rvalue.push_back(9 - 3 * (202 % 2));
			else
				rvalue.push_back(10 - 2 * (i % 2));
		}
}

	void printC()
	{
		for (int j = 0; j != 3; j++)
		{
			for (int i = 0; i != count; i++)
			{
				std::cout << value[j][i] << " ";
			}
			std::cout << std::endl;
		}
	for (int i = 0; i != 202; i++)
		std::cout << rvalue[i] << " ";
	}

	void init()
	{
		std::ifstream file;
		file.open("matrix.txt");

		for (int i = 0; i != count; i++)
		{
			for (int j = 0; j <= count; j++)
			{
				T numb;
				file >> numb;
				if (j != count)
					value[i][j] = numb;
				else
					rvalue.push_back(numb);
			}
		}
		file.close();
	}

	void cpy(const VectorMatrix <T> &A)
	{
		for (int i = 0; i != A.count; i++)
		{
			for (int j = 0; j != A.count; j++)
			{
				value[i][j] = A.value[i][j];
			}
		}
	}

	void onebyone()
	{
		std::cout << "I'm right!";
		for (int i = 0; i != count; i++)
		{
			for (int j = 0; j != count; j++)
			{
				if (i == j)
					value[i][j] = 1;
				else
					value[i][j] = 0;
			}
		}
	}
	void trunc()
	{
		T temp;
		for (int i = 0; i != count; i++)
		{
			for (int j = i + 1; j != count; j++)
			{
				temp = value[i][j];
				value[i][j] = value[j][i];
				value[j][i] = temp;
			}
		}
	}
	void print()
	{
		for (int n = 0; n != count; n++)
		{
			std::cout << std::endl;;
			for (int j = 0; j != count; j++)
			{
				std::cout << std::setw(8) << value[n][j] << "\t";
			}
		}
		std::cout << std::endl;;
	}
	// friend VectorMatrix VectorMatrix::operator *(VectorMatrix A, VectorMatrix B)
	// {
	//   VectorMatrix <T> C
	//   for (int i = 0; i != A.count; i++)
	//   {
	//     T sum = 0;
	//     for (int j = 0; j != A.count; i++)
	//     {
	//       for (inr k = 0; k != A.count; k++)
	//       {
	//         sum += A.matrix[i][j] * B.matrix[k][j];
	//       }
	//       C.matrix[][] = sum;
	//     }
	//   }
	// }
	friend	mytype	norm_C_cube(const VectorMatrix <mytype> &A);
	friend	mytype	norm_C_oct(const VectorMatrix <mytype> &A);
	friend	std::vector<mytype>	diff_vector(std::vector<mytype>a, std::vector<mytype>b, int DIM);
	friend	mytype	cube_vect_norm(std::vector<mytype>x, int size);
	friend	mytype	octah_vect_norm(std::vector<mytype>x, int size);
	friend	std::vector<mytype> cpy_vector(std::vector<mytype>x, int size);
	friend	std::vector<mytype> multi_vect_matr_C(const VectorMatrix <mytype> &A, std::vector<mytype>x);
	friend	std::vector<mytype> sum_vect(std::vector<mytype>a, std::vector<mytype>b, std::vector<mytype>c, int size);
	friend	std::vector<mytype> simp_iter(const VectorMatrix <mytype> &A, std::vector<mytype>b, int flag);
	friend void printC(const VectorMatrix& C);
	friend VectorMatrix<mytype> countbigmatrix(const VectorMatrix<mytype> & A, mytype omega);
	friend std::vector<mytype> BigRelax(const VectorMatrix<mytype> &A);
	friend std::vector<mytype> BigZeldel(const VectorMatrix<mytype> & A);
	friend std::vector<mytype> SmallZeydel(const VectorMatrix<mytype> & A);
	friend VectorMatrix<mytype> countmatrix(const VectorMatrix<mytype>& A, double omega);
	friend std::vector<mytype> SmallRelax(const VectorMatrix<mytype> & A);
	friend mytype norm_C_big_oct(const VectorMatrix <mytype> &C);
	friend mytype norm_C_big_cube(const VectorMatrix<mytype> &C);
	friend mytype Ufunc(const VectorMatrix<mytype> &A, bool flag);
	friend mytype Lfunc(const VectorMatrix<mytype> &A, bool flag);
	friend VectorMatrix<mytype> countbigmatrix(const VectorMatrix<mytype> &A, mytype omega);
	friend mytype searchbigomeganorm(const VectorMatrix<mytype> &A);
	friend mytype searchbigomegaiter(VectorMatrix<mytype> &A);
};
//double *abs_diff_vector(T *a, T *b, int DIM)
