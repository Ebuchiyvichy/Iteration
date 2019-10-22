
#include "Estimation.h"
  void delim(VectorMatrix<double>  &A, int k)
{
 int i = k;
 double temp = A.value[i][k];
 for (int j = k; j != A.count; j++)
 {
   A.value[i][j] = A.value[i][j] / temp;
 }
 A.rvalue[i] = A.rvalue[i] / temp;
}

  void vych(VectorMatrix <double> &A, int k)
{
 for (int i = (k + 1); i != A.count; i++)
 {
   //columns
   double temp = A.value[i][k];
   for (int j = 0; j != A.count; j++)
     A.value[i][j] = A.value[i][j] - temp * A.value[k][j];
   A.rvalue[i] -= temp * A.rvalue[k];
 }
}

  void GaussRight(VectorMatrix <double>  &A)
{
 // for k-string our matrix
 for (int k = 0; k != A.count; k++)
 {
   //search max in column
  double max = A.value[k][k];
   int maxstring = k;
   for (int i = k; i != A.count; i++)
   {
     if (A.value[i][k] > max)
     {
       max = A.value[i][k];
       maxstring = i;
     }
   }
   //change rvalue
   double temp = A.rvalue[k];
   A.rvalue[k] = A.rvalue[maxstring];
   A.rvalue[maxstring] = temp;
   // change strings
   std::swap(A.value[k], A.value[maxstring]);
   delim(A, k);
   vych(A, k);
 }
 std::cout << "Matrix A after straight run of Gauss:" << std::endl;
 A.print();
}

  void GaussLeft(VectorMatrix<double> &A)
{
 for (int i = A.count - 1; i >= 0; i--)
 {
   for (int j = i + 1; j < A.count; j++)
     A.rvalue[i] = A.rvalue[i] - A.value[i][j] * A.rvalue[j];
   A.rvalue[i] = A.rvalue[i] / A.value[i][i];
 }
 std::cout << "Vector in Gauss x:" << std::endl;
 for (int i = 0; i != A.count; i++)
   std::cout << std::setw(8) << A.rvalue[i] << std::endl;
 std::cout << std::endl;
}
