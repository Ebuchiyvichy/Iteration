#include "Jacoby.h"

void printC(const VectorMatrix<mytype> &C)
{
	for (int j = 0; j != 3; j++)
	{
		for (int i = 0; i != C.count; i++)
			std::cout << C.value[j][i] << " ";
		std::cout << std::endl;
	}
}

//mytype norm_C_big_cube(const VectorMatrix<mytype> &C)
//{
//	mytype	norm = 0;
//	mytype	max = 0;
//	for (int i = 0; i != C.count; i++)
//	{
//		norm = 0;
//		for (int j = 0; j != 3; j++)
//			norm += fabs(C.value[j][i]);
//		if (norm > max)
//			max = norm;
//	}
//	return(max);
//}

//mytype norm_C_big_oct(const std::vector<std::vector<mytype>> &C)
//{
//	mytype	norm = 0;
//	mytype	max = 0;
//
//	for (int i = 0; i != C.count; i++)
//	{
//		norm = 0;
//		if (i == 0)
//			norm = fabs(C.value[1][i] + C.value[2][i + 1]);
//		else if (i == 201)
//			norm = fabs(C.value[1][i] + C.value[0][i - 1]);
//		else
//			norm = fabs(C.value[1][i] + C.value[0][i - 1] + C.value[2][i + 1]);
//		if (norm > max)
//			max = norm;
//	}
//	return(max);
//}

 //std::vector<std::vector<mytype>> countbigmatrix(const VectorMatrix<mytype> & A, mytype omega)
 //{
	// std::vector<std::vector<mytype>>	C;

 //  for (int k = 0; k != A.count; k++)
 //  {
 //    std::vector<mytype> x(A.count);
 //    std::vector<mytype> p(A.count);
 //    x[k] = 1;
 //    p[k] = 1;
 //    for (int i = 0; i != A.count; i++)
 //    {
 //      switch (i) {
 //        case 0: x[i] = (-1)*omega*p[i+1] *A.value[i][2]/A.value[i][1] +(1-omega) * x[i];
 //        case 201: x[i] = (-1)*omega*x[i-1] *A.value[i][0]/A.value[i][1]+(1-omega) * x[i];
 //        default:
 //            x[i] = ((-1)*p[i+1]*A.value[i][2] - x[i-1]*A.value[i][0])/A.value[i][1] + (1-omega) * x[i];
 //          }
 //    }
 //    for (int i = 0; i != A.count; i++)
 //    {
 //      if (x[i] != 0)
 //      {
 //        for (int l = 0; l != 3; l++)
 //          C.value[i][l]=x[i];
 //        break ;
 //      }
 //    }
 //  }
 //  return C;
 //}

// std::vector<mytype> BigRelax(const VectorMatrix<mytype> &A )
// {
//   std::vector<mytype> x(A.rvalue);
//   std::vector<mytype> p(A.rvalue);
//   double omega = 0.01;
//   //start iteration process
//   do
//   {
//     for (int i = 0; i != A.count; i++)
//       p[i] = x[i];
//     for (int i = 0; i != A.count; i++)
//     {
//       switch (i) {
//         case 0: x[i] = ((-1)*omega*p[i+1] *A.value[i][2]+omega*A.rvalue[i])/A.value[i][1] +(1-omega) * p[i];
//         case 201: x[i] = ((-1)*omega*x[i-1] *A.value[i][0]+ omega*A.rvalue[i])/A.value[i][1] +(1-omega) * p[i];
//         default:
//             x[i] = ((-1)*omega*p[i+1]*A.value[i][2] - omega*x[i-1]*A.value[i][0]+ omega * A.rvalue[i])/A.value[i][1] + (1-omega) * p[i];
//           }
//     }
//   }
//   while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_cube(C)) / norm_cube(C)) * EPS);
//   return x;
// }


 //std::vector<mytype> BigZeydel(const VectorMatrix<mytype> & A)
 //{
 //  std::vector<mytype> x(A.rvalue);
 //  std::vector<mytype> p(A.rvalue);
 //  //start iteration process
 //  do
 //  {
 //    for (int i = 0; i != A.count; i++)
 //      p[i] = x[i];
 //    for (int i = 0; i != A.count; i++)
 //    {
 //      switch (i)
	//   {
 //        case 0: 
	//		 x[i] = ((-1)*p[i+1] *A.value[i][2]+A.rvalue[i])/A.value[i][1];
 //        case 201:
	//		 x[i] = ((-1)*x[i-1] *A.value[i][0]+A.rvalue[i])/A.value[i][1];
 //        default:
 //            x[i] = ((-1)*p[i-1]*A.value[i][2] - x[i-1]*A.value[i][0]+ A.rvalue[i])/A.value[i][1];
 //      }
 //    }
 //  }
 // // while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_C_big_cube(C)) / norm_cube(C)) * EPS);
 //  while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > EPS);
 //  return x;
 //}

std::vector<mytype> SmallZeydel(const VectorMatrix<mytype> & A, int flag)
{
	std::vector <mytype>	x(A.rvalue);
	std::vector <mytype>	p(A.rvalue);
	VectorMatrix <mytype>	C(A.count);
	int						iter = 0;
	//start iteration process
	if (flag == 0)
	{
		do
		{
			for (int i = 0; i != A.count; i++)
				p[i] = x[i];
			for (int i = 0; i != A.count; i++)
			{
				mytype	sum = 0;
				for (int j = 0; j < i; j++)
				{
					sum += (A.value[i][j] * x[j]);
					C.value[i][j] = (-1)*A.value[i][j] / A.value[i][i];
				}
				C.value[i][i] = 0;
				for (int j = i + 1; j != A.count; j++)
				{
					sum += (A.value[i][j] * p[j]);
					C.value[i][j] = (-1)*A.value[i][j] / A.value[i][i];
				}
				x[i] = (A.rvalue[i] - sum) / A.value[i][i];
			}
			iter++;
			std::cout << "Residual vector with cube norme at " << iter << " step = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
		} while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > (fabs(1 - norm_cube(C)) / norm_cube(C)) * EPS);
	//	} while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > EPS);
	//	} while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > 10 * EPS);
	//	} while (cube_vect_norm(diff_vector(multi_vect(x, A), A.rvalue, A.count), A.count) > EPS);
		std::cout << "Cube norme of matrix C = " << norm_cube(C) << std::endl;
	}
	else if (flag == 1)
	{
		do
		{
			for (int i = 0; i != A.count; i++)
				p[i] = x[i];
			for (int i = 0; i != A.count; i++)
			{
				mytype	sum = 0;
				for (int j = 0; j < i; j++)
				{
					sum += (A.value[i][j] * x[j]);
					C.value[i][j] = (-1)*A.value[i][j] / A.value[i][i];
				}
				C.value[i][i] = 0;
				for (int j = i + 1; j != A.count; j++)
				{
					sum += (A.value[i][j] * p[j]);
					C.value[i][j] = (-1)*A.value[i][j] / A.value[i][i];
				}
				x[i] = (A.rvalue[i] - sum) / A.value[i][i];
			}
			iter++;
			std::cout << "Residual vector with octahedral norme at " << iter << " step = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
		} while (octah_vect_norm(diff_vector(x, p, A.count), A.count) > (fabs(1 - norm_oct(C)) / norm_oct(C)) * EPS);
		std::cout << "Octahedral norme of matrix C = " << norm_oct(C) << std::endl;
	}
	std::cout << "Number of iteration: " << iter << std::endl;
	return x;
}

VectorMatrix<mytype> countmatrix(const VectorMatrix<mytype> &A, mytype omega)
{
	VectorMatrix<mytype> C(A.count);
	for (int k = 0; k != A.count; k++)
	{
		std::vector<mytype> x(A.count);
		std::vector<mytype> p(A.count);
		x[k] = 1;
		p[k] = 1;
		// for (int i = 0; i != A.count; i++)
		//   p[i] = x[i];
		for (int i = 0; i != A.count; i++)
		{
			mytype	sum = 0;
			for (int j = 0; j < i; j++)
				sum += (omega * A.value[i][j] * x[j]) / A.value[i][i];

			for (int j = i + 1; j != A.count; j++)
				sum += (omega * A.value[i][j] * p[j]) / A.value[i][i];
			x[i] = (1 - omega) * p[i] - sum;
		}
		for (int i = 0; i != A.count; i++)
		{
			C.value[i][k] = x[i];
		}
	}
	return (C);
}

std::vector<mytype> SmallRelax(const VectorMatrix<mytype> &A, int flag)
{
	mytype	omega = 1;
	std::vector<mytype> x(A.rvalue);
	std::vector<mytype> p(A.rvalue);
	int					iter = 0;
	//VectorMatrix<mytype> C(countmatrix(A, omega));
	std::wcout << "Omega = " << omega << std::endl;
	if (flag == 0)
	{
		do
		{
			for (int i = 0; i != A.count; i++)
				p[i] = x[i];
			for (int i = 0; i != A.count; i++)
			{
				x[i] = (A.rvalue[i] / A.value[i][i] - p[i]) * omega + p[i];
				mytype	sum = 0;
				for (int j = 0; j < i; j++)
					sum -= (omega * A.value[i][j] * x[j]) / A.value[i][i];
				for (int j = i + 1; j != A.count; j++)
					sum -= (omega * A.value[i][j] * p[j]) / A.value[i][i];
				x[i] += sum;
			}
			iter++;
			std::cout << "Residual vector with cube norme at " << iter << " step = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
		} while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > (fabs(1 - norm_cube(countmatrix(A, omega)) / norm_cube(countmatrix(A, omega))) * EPS));
	//	} while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > EPS);
	//	} while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > 10 * EPS);
	//	} while (cube_vect_norm(diff_vector(multi_vect(x, A), A.rvalue, A.count), A.count) > EPS);
		std::cout << "Cube norme of matrix C = " << norm_cube(countmatrix(A, omega)) << std::endl;
	}
	else if (flag == 1)
	{
		do
		{
			for (int i = 0; i != A.count; i++)
				p[i] = x[i];
			for (int i = 0; i != A.count; i++)
			{
				x[i] = (A.rvalue[i] / A.value[i][i] - p[i]) * omega + p[i];
				mytype	sum = 0;
				for (int j = 0; j < i; j++)
					sum -= (omega * A.value[i][j] * x[j]) / A.value[i][i];
				for (int j = i + 1; j != A.count; j++)
					sum -= (omega * A.value[i][j] * p[j]) / A.value[i][i];
				x[i] += sum;
			}
			iter++;
			std::cout << "Residual vector with octahedral norme at " << iter << " step = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
		} while (octah_vect_norm(diff_vector(x, p, A.count), A.count) > (fabs(1 - norm_oct(countmatrix(A, omega)) / norm_oct(countmatrix(A, omega))) * EPS));
		std::cout << "Octahedral norme of matrix C = " << norm_oct(countmatrix(A,omega)) << std::endl;
	}
	std::cout << "Number of iteration: " << iter << std::endl;
	return x;
}