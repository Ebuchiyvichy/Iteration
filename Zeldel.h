#include "Jacoby.h"

void printC(const VectorMatrix<mytype> &C)
{
  for (int i = 0; i != C.count; i++)
  {
    for (int j = 0; j != 3; j++)
      std::cout << C.value[j][i] << "   ";
    std::cout << std::endl;
    }
  }

// VectorMatrix<mytype> countbigmatrix(const VectorMatrix<mytype> & A, mytype omega)
// {
//   for (int k = 0; k != A.count; k++)
//   {
//     std::vector<mytype> x(A.count);
//     std::vector<mytype> p(A.count);
//     x[k] = 1;
//     p[k] = 1;
//     for (int i = 0; i != A.count; i++)
//     {
//       switch (i) {
//         case 0: x[i] = (-1)*omega*p[i+1] *A.value[i][2]/A.value[i][1] +(1-omega) * x[i];
//         case 201: x[i] = (-1)*omega*x[i-1] *A.value[i][0]/A.value[i][1]+(1-omega) * x[i];
//         default:
//             x[i] = ((-1)*p[i+1]*A.value[i][2] - x[i-1]*A.value[i][0])/A.value[i][1] + (1-omega) * x[i];
//           }
//     }
//     for (int i = 0; i != A.count; i++)
//     {
//       if (x[i] != 0)
//       {
//         for (int l = 0; l != 3; l++)
//           C.value[i][l]=x[i];
//         break ;
//       }
//     }
//   }
//   return C;
// }

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

//
// std::vector<mytype> BigZeldel(const VectorMatrix<mytype> & A)
// {
//   std::vector<mytype> x(A.rvalue);
//   std::vector<mytype> p(A.rvalue);
//   //start iteration process
//   do
//   {
//     for (int i = 0; i != A.count; i++)
//       p[i] = x[i];
//     for (int i = 0; i != A.count; i++)
//     {
//       switch (i) {
//         case 0: x[i] = ((-1)*p[i+1] *A.value[i][2]+A.rvalue[i])/A.value[i][1];
//         case 201: x[i] = ((-1)*x[i-1] *A.value[i][0]+A.rvalue[i])/A.value[i][1];
//         default:
//             x[i] = ((-1)*p[i-1]*A.value[i][2] - x[i-1]*A.value[i][0]+ A.rvalue[i])/A.value[i][1];
//           }
//     }
//   }
//   while (cube_vect_norm(diff_vector(x,p,A.count), A.count) > ((1 - norm_cube(C)) / norm_cube(C)) * EPS);
//   return x;
// }

std::vector<mytype> SmallZeydel(const VectorMatrix<mytype> & A)
{
  std::vector<mytype> x(A.rvalue);
  std::vector<mytype> p(A.rvalue);
  VectorMatrix<mytype> C(A.count);
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
    // for (int i = 0; i != A.count; i++)
    //   p[i] = x[i];
    for (int i = 0; i != A.count; i++)
    {
      double sum = 0;
      for (int j = 0; j < i; j++)
        sum +=(omega * A.value[i][j] * x[j])/A.value[i][i];

      for (int j = i + 1; j != A.count; j++)
        sum += (omega * A.value[i][j] * p[j])/A.value[i][i];
      x[i]= (1 - omega) * p[i] - sum;
    }
    for (int i = 0; i != A.count; i++)
    {
      C.value[i][k] = x[i];
    }
  }
  return (C);
}

 std::vector<mytype> SmallRelax(const VectorMatrix<mytype> & A)
 {
   double omega = 0;
   std::vector<mytype> x(A.rvalue);
   std::vector<mytype> p(A.rvalue);
   //VectorMatrix<mytype> C(countmatrix(A, omega));
   do
   {
     int n = 0;
     while (norm_cube(countmatrix(A, omega)) >= 1)
     {
       std::cout << norm_cube(countmatrix(A, omega)) << "\n";
       n++;
       std::cout << omega << std::endl;
     }
     std::cout << norm_cube(countmatrix(A, omega)) << "\n";
     std::cout << "Meow\n";

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
   } while (cube_vect_norm(diff_vector(x,p, A.count), A.count) > ((1 - norm_cube(countmatrix(A, omega)) / norm_cube(countmatrix(A, omega))) * EPS));
   return x;
 }
