#include "Jacoby.h"

mytype Lfunc(const VectorMatrix<mytype> &A, bool flag)
{
	mytype	norm = 0;
	mytype	max = 0;
	//true if small matrix, false if else
	if (flag == true)
	{
		for (int i = 0; i < A.count; i++)
		{
			norm = 0;
			for (int j = i+1; j < A.count; j++)
				norm += fabs(A.value[i][j]);
			if (norm > max)
				max = norm;
		}
	}
	else
	{
		for (int i = 0; i != A.count; i++)
		{
			if (fabs(A.value[2][i]) > max)
				max = fabs(A.value[2][i]);
		}
	}
	return(max);
}
mytype Ufunc(const VectorMatrix<mytype> &A, bool flag)
{
	mytype	norm = 0;
	mytype	max = 0;
	//true if small matrix, false if else
	if (flag == true)
	{
		for (int i = 0; i < A.count; i++)
		{
			norm = 0;
			for (int j = 0; j < i; j++)
				norm += fabs(A.value[i][j]);
			if (norm > max)
				max = norm;
		}
	}
	else
	{
		for(int i = 0; i != A.count; i++)
		{
			if (fabs(A.value[0][i]) > max)
				max = fabs(A.value[0][i]);
		}
	}
	return(max);
}

mytype norm_C_big_cube(const VectorMatrix<mytype> &C)
{
	mytype	norm = 0;
	mytype	max = 0;
	for (int i = 0; i != C.count; i++)
	{
		norm = 0;
		for (int j = 0; j != 3; j++)
			norm += fabs(C.value[j][i]);
		if (norm > max)
			max = norm;
	}
	return(max);
}

mytype norm_C_big_oct(const VectorMatrix<mytype> &C)
{
	mytype	norm = 0;
	mytype	max = 0;

	for (int i = 0; i != C.count; i++)
	{
		norm = 0;
		if (i == 0)
			norm = fabs(C.value[1][i] + C.value[2][i + 1]);
		else if (i == 201)
			norm = fabs(C.value[1][i] + C.value[0][i - 1]);
		else
			norm = fabs(C.value[1][i] + C.value[0][i - 1] + C.value[2][i + 1]);
		if (norm > max)
			max = norm;
	}
	return(max);
}

std::vector<mytype> BigRelax(const VectorMatrix<mytype> &A )
{
  std::vector<mytype> x(A.rvalue);
  std::vector<mytype> p(A.rvalue);
	std::cout << "cool!\n";
  mytype omega = 1.5;
  //start iteration process
  do
  {
    for (int i = 0; i != A.count; i++)
      p[i] = x[i];
    for (int i = 0; i != A.count; i++)
    {
      if (i == 0)
          x[i] = ((-1)*omega*p[i+1] *A.value[2][i]+omega*A.rvalue[i])/A.value[1][i] +(1-omega) * p[i];
      else if ( i == 201)
          x[i] = ((-1)*omega*x[i-1] *A.value[0][i]+ omega*A.rvalue[i])/A.value[1][i] +(1-omega) * p[i];
      else
          x[i] = ((-1)*omega*p[i+1]*A.value[2][i] - omega*x[i-1]*A.value[0][i]+ omega * A.rvalue[i])/A.value[1][i] + (1-omega) * p[i];
      }
    }
  while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > EPS);
  //while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_cube(countbigmatrix(A, omega))) /  norm_cube(countbigmatrix(A, omega)) * EPS));
  return x;
}


 std::vector<mytype> BigZeydel(const VectorMatrix<mytype> & A)
 {
  std::vector<mytype> x(A.rvalue);
  std::vector<mytype> p(A.rvalue);
	VectorMatrix<mytype> C(0);
	C.VectorMatrixC(202);
	for (int i = 0; i != A.count; i++)
	{
		C.value[1][i] = 0;
	}
  //start iteration process
  do
  {
		for (int i = 0; i != A.count; i++)
			p[i] = x[i];
		for (int i = 0; i != A.count; i++)
		{
			if (i == 0)
			{
				x[i] = ((-1) * p[i+1] * A.value[2][i] + A.rvalue[i])/A.value[1][i];
				C.value[2][i] = A.value[2][i]/A.value[1][i];
			}
			else if (i == 201)
			{
				x[i] = ((-1) * x[i-1] * A.value[0][i] + A.rvalue[i])/A.value[1][i];
				C.value[0][i] = A.value[0][i]/A.value[1][i];
			}
			else
			{
				x[i] = ((-1) * x[i-1] * A.value[0][i] - p[i+1] * A.value[2][i] + A.rvalue[i])/A.value[1][i];
				C.value[0][i] = A.value[0][i]/A.value[1][i];
				C.value[2][i] = A.value[2][i]/A.value[1][i];
			}
		}
  }
 //while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_C_big_cube(C)) / norm_C_big_cube(C)) * EPS);
  while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > EPS);
  return x;
 }

mytype searchbigomeganorm(const VectorMatrix<mytype> &A)
{
	mytype omega = 0.01;
	mytype omegaopt = omega;
	double h = 0.01;
	mytype min = norm_cube(countmatrix(A, omega));
	while (omega < 2){
		if (norm_cube(countmatrix(A, omega)) < min )
		{
			std::cout << min;
			min = norm_cube(countmatrix(A, omega));
			omegaopt = omega;
	 	}
		omega += h;
	}
	std::cout << "Norm matrix C: " << min << std::endl;
	std::cout << "Omega opt: " << omegaopt << std::endl;
	return(omegaopt);
}

mytype searchbigomegaiter(VectorMatrix<mytype> &A)
{
	mytype omega = 1;
	double h = 0.5;
	int miniter;
	mytype omegaopt;
	while(omega < 2){
		std::cout << "workwork\n";
		int iter = 0;
		std::vector<mytype> x(A.rvalue);
	  std::vector<mytype> p(A.rvalue);
		std::cout << "workwork\n";
	  //start iteration process
//	  do {
//			std::cout << "workwork\n";
//			iter++;
//	    for (int i = 0; i != A.count; i++)
//	      p[i] = x[i];
//	    for (int i = 0; i != A.count; i++)
//	    {
//				mytype	sum = 0;
//				if (i == 0)
//					sum += (omega * (A.value[2][i] * p[i+1] - A.rvalue[i]))/A.value[1][i];
//			else if ( i == 201)
//					sum += (omega * (A.value[0][i] * x[i-1] - A.rvalue[i]))/A.value[1][i];
//				else
//					sum += (omega * (A.value[2][i] * p[i+1] + A.value[0][i] * x[i-1] - A.rvalue[i]))/A.value[1][i];
//				x[i] = (1 - omega) * p[i] - sum;
//				}
//	    }
//		//while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_cube(countbigmatrix(A, omega))) /  norm_cube(countbigmatrix(A, omega)) * EPS));
//      while (cube_vect_norm(diff_vector(multi_vect(x, A), A.rvalue, A.count), A.count) > EPS);
        do
        {
            for (int i = 0; i != A.count; i++)
                p[i] = x[i];
            for (int i = 0; i != A.count; i++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum -=(omega * A.value[i][j] * x[j])/A.value[i][i];
                for (int j = i + 1; j != A.count; j++)
                    sum -= (omega * A.value[i][j] * p[j])/A.value[i][i];
                sum += ((1-omega) * p[i]);
                x[i]= ((omega*A.rvalue[i])/A.value[i][i]+sum);
            }
            iter++;
        } while (cube_vect_norm(diff_vector(x, p, A.count), A.count) > (fabs(1 - norm_cube(countmatrix(A, omega)) / norm_cube(countmatrix(A, omega))) * EPS));
		if (omega == 0.1)
			miniter = iter;
		if (iter < miniter)
		{
			miniter = iter;
			omegaopt = omega;
		}
		omega += h;
	}
	std::cout << "Count of iteration: " << miniter <<"\n";
	std::cout << "Omega opt: " << omegaopt << "\n";
	return (omegaopt);
}

VectorMatrix<mytype> countbigmatrix(const VectorMatrix<mytype> &A, mytype omega)
{
	VectorMatrix<mytype> C(202);
	for (int k = 0; k != A.count; k++)
	{
		std::vector<mytype> x(A.count);
		std::vector<mytype> p(A.count);
		x[k] = 1;
		p[k] = 1;
		for (int i = 0; i != A.count; i++)
		{
			mytype	sum = 0;
			if (i == 0)
				sum += (omega * (A.value[2][i] * p[i+1] - A.rvalue[i]))/A.value[1][i];
		else if ( i == 201)
				sum += (omega * (A.value[0][i] * x[i-1] - A.rvalue[i]))/A.value[1][i];
			else
				sum += (omega * (A.value[2][i] * p[i+1] + A.value[0][i] * x[i-1] - A.rvalue[i]))/A.value[1][i];
			x[i] = (1 - omega) * p[i] - sum;
			}

		for (int i = 0; i != A.count; i++)
		{
			C.value[i][k] = x[i];
		}
	}

	return (C);
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
		for (int i = 0; i != A.count; i++)
		{
			mytype	sum = 0;
			for (int j = 0; j < i; j++)
				sum += (omega * A.value[i][j] * x[j]) / A.value[i][i];

			for (int j = i + 1; j != A.count; j++)
				sum += (omega * A.value[i][j] * p[j] ) / A.value[i][i];
			x[i] = (1 - omega) * p[i] - sum + omega * A.rvalue[i]/A.value[i][i];
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

	for (int t = 0; t != A.count; t++) {
    x[t] = 1;
    p[t] = 1;
    }
	//VectorMatrix<mytype> C(countmatrix(A, omega));
	std::cout << "Omega = " << omega << std::endl;
	if (flag == 0)
	{
		do
		{
			for (int i = 0; i != A.count; i++)
				p[i] = x[i];
            for (int i = 0; i != A.count; i++)
            {
                mytype sum = 0;
                for (int j = 0; j < i; j++)
                    sum -=(omega * A.value[i][j] * x[j])/A.value[i][i];
                for (int j = i + 1; j != A.count; j++)
                    sum -= (omega * A.value[i][j] * p[j])/A.value[i][i];
                sum += ((1-omega) * p[i]);
                x[i]= ((omega*A.rvalue[i])/A.value[i][i]+sum);
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
	countmatrix(A,omega).print();
	std::cout << "Number of iteration: " << iter << std::endl;
	return x;
}
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
		C.print();
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
