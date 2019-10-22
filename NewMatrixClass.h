
#include <iostream>
#include <iomanip>
#include <fstream>

double	EPS = 10e-8;

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
    std::cout << "OK";
  }
   ~VectorMatrix()
   {
   }
   void init()
   {
     std::ifstream file;
     file.open("matrixexample.txt");

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
       for (int j = i+1; j != count; j++)
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
       std::cout << "\n";
       for (int j = 0; j != count; j++)
       {
         std::cout << value[n][j] << "  ";
       }
     }
     std::cout << "\n";
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

   friend  void delim(VectorMatrix <T> &A, int k);
   friend void vych(VectorMatrix <T> &A, int k);
   friend void GaussRight(VectorMatrix<T> &A);
   friend void GaussLeft(VectorMatrix<T> &A);
   void	QR_find_x(VectorMatrix<T>& A);
 	 void	inverse_matrix(const VectorMatrix<T>& R,VectorMarix<T> & G);
 	 friend  double *find_x(const VectorMatrix<T>& A, double *b_);
   friend bool		degenerate_matrix(const VectorMatrix<T>& R);
	 double *multi_vect(std::vector<double> I, const VectorMatrix  <double> & F);
};
double *abs_diff_vector(T *a, T *b, int DIM);
