#include "Jacoby.h"

VectorMatrix<mytype>::VectorMatrix<typename> randombigmatrix()
{
  count = 204;
  std::vector<std::vector <mytype>> value;

  for (int i = 0; i != 3; i++)
  {
    for (int j = 0; j != count; j++)
  {
    if (i > 3)    
  }
}
 }

std::vector<mytype> SmallZeydel(const VectorMatrix<mytype> & A)
{
  std::vector<mytype> x(A.rvalue);
  std::vector<mytype> p(A.rvalue);
  VectorMatrix<mytype> C(A.count);
  std::vector<mytype> c(A.count);

  //start iteration process
  do
  {
    for (int i = 0; i != A.count; i++)
      p[i] = x[i];
    for (int i = 0; i != A.count; i++)
    {
      double sum = 0;
      for (int j = 0; j < i; j++)
      {
        sum +=(A.value[i][j]*x[j]);
        C.value[i][j] = (-1)*A.value[i][j]/A.value[i][i];
      }
      C.value[i][i] = 0;
      for (int j = i + 1; j != A.count; j++)
      {
        sum += (A.value[i][j] * p[j]);
        C.value[i][j] = (-1)*A.value[i][j]/A.value[i][i];
      }
      x[i]= (A.rvalue[i]-sum) / A.value[i][i];
    }
  }
  while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_cube(C)) / norm_cube(C)) * EPS);
  return x;
}


VectorMatrix<mytype> countmatrix(const VectorMatrix<mytype>& A, double omega)
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
      for (int j = 0; j!= A.count; j++)
        p[j] = x[j];
      x[i] = A.rvalue[i];
      for (int j = 0; j < A.count; j++)
      {
        if (i != j)
          x[i] = x[i] - A.value[i][j]*x[j];
      }
      x[i] /= A.value[i][i];
      x[i] = omega * x[i] + (1 - omega)*p[i];
    }
    for (int j = 0; j != A.count; j++)
    {
      C.value[j][k] = x[j];
    }
  }
  C.print();
  return C;
}

 std::vector<mytype> SmallRelax(const VectorMatrix<mytype> & A)
 {
   double omega = 2;
   std::vector<mytype> x(A.rvalue);
   std::vector<mytype> p(A.rvalue);
   VectorMatrix<mytype> C(countmatrix(A, omega));
   do
   {
     int n =0;
     while (norm_cube(countmatrix(A, omega)) > 1)
     {
       n++;
       omega -= 0.01;
       std::cout << omega << std::endl;
       if (n == 20)
       {
         std::cout << norm_cube(countmatrix(A, omega));
         break ;
       }
     }
     for (int i = 0; i != A.count; i++)
       p[i] = x[i];
     for (int i = 0; i != A.count; i++)
     {
      double sum = 0;
      for (int j = 0; j < i; j++)
        sum -=(omega*A.value[i][j]*x[j])/A.value[i][i];
      sum += ((1-omega) * p[i]);
      x[i]= ((omega*A.rvalue[i])/A.value[i][i]+sum);
     }
   } while (cube_vect_norm(diff_vector(x,p, A.count), A.count) > ((1 - norm_cube(countmatrix(C, omega)) / norm_cube(countmatrix(C, omega))) * EPS));
   return x;
 }
