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
		std::cout << std::endl;
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
	// friend std::vector<T>& operator -(const std::vector<T> & a, const std::vector<T> & b)
	// {
	// 	std::vector<T> c(a);
	// 	for (int i = 0; i != c.size(); i++)
	// 		c[i]=a[i]-b[i];
	// 	return c;
	// }
	void randombigmatrix()
	{
		count = 202;
	  std::vector<std::vector <mytype>> value;

	  for (int i = 0; i != 3; i++)
	  {
	    for (int j = 0; j != 202; j++)
	    {
	      if (i > 3)
	        value[i].push_back(1);
	      switch (i)
	      {
	        case 0: {
						if (j == 0)
							value[i][j] = 0;
						else value[i][j] = 1;}
	        case 1: value[i][j] = 4;
	        case 2: value[i][j] = 1;
	      }
	    }
	  }
	  for (int i = 0; i != count; i++)
	  {
	    if (i > 3)
	      rvalue.push_back();
	    switch (i) {
	      case 0: rvalue[i] = 6;
	      case (count - 1): rvalue[i] = 9 - 3 * (count % 2);
			}
		}

		value[2][201] = 0;
	}
	friend std::vector<mytype> SmallRelax(const VectorMatrix<mytype> & A);
	friend VectorMatrix<mytype> countmatrix(const VectorMatrix<mytype>& A, double omega);
	friend	mytype	norm_C_cube(const VectorMatrix <mytype> &A);
	friend	mytype	norm_C_oct(const VectorMatrix <mytype> &A);
	friend std::vector<mytype> diff_vector(std::vector<mytype> a, std::vector<mytype> b, int DIM);
	friend	mytype	cube_vect_norm(std::vector<mytype> x, int size);
	friend	mytype	octah_vect_norm(std::vector<mytype> x, int size);
	friend	std::vector<mytype> cpy_vector(std::vector<mytype> x, int size);
	friend	std::vector<mytype>	multi_vect_matr_C(const VectorMatrix <mytype> &A, std::vector<mytype> x);
	friend	std::vector<mytype> sum_vect(std::vector<mytype> a, std::vector<mytype> b, std::vector<mytype> c, int size);
	friend	std::vector<mytype> simp_iter(const VectorMatrix <mytype> &A, std::vector<mytype> b, int flag);
	friend std::vector<mytype> multi_vect(std::vector<double> I, const VectorMatrix  <double> & F);
};
//double *abs_diff_vector(T *a, T *b, int DIM);
