<<<<<<< HEAD

#include "Zeldel.h"

int	main()
{
//	std::ifstream file;
//	file.open("matrixexample.txt");
//	int	size;
//	file >> size;
	VectorMatrix <mytype>	A(4);
	std::vector<mytype> x;
=======
#include "Gauss.cpp"
int	main()
{
>>>>>>> a76e75838e30bcd6d8429795141cb69842d5a2ed

	VectorMatrix <double>	A(4);
	A.init();
	VectorMatrix <double>	R(4);
	R.init();
	//VectorMatrix	A_cpy;	//copy of A matrix for inverse Gauss
	VectorMatrix <double>	T(R.count);
	T.onebyone();
	std::vector <double> b_;
	std::vector <double> x;

	std::cout << "matrix A:" << std::endl;
	A.print();
<<<<<<< HEAD
	std::cout << "Zeldel for small one method:" << std::endl;
	// std::cout << "Simple iteration method:" << std::endl;
	x = SmallZeydel(A);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
  std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(x, multi_vect(A.rvalue, A), A.count), A.count) << std::endl;
	 std::cout << std::endl;
	//  std::cout << "Relax for small one method:" << std::endl;
	//  x = SmallRelax(A);
 	// std::cout << "Vector X with cube norme:" << std::endl;
 	// for (int i = 0; i < A.count; i++)
  //  std::cout << std::setw(8) << x[i] << std::endl;
	//x = simp_iter(A, 1);
	//std::cout << "Vector X with octahedral norme:" << std::endl;
	//for (int i = 0; i < A.count; i++)
	//	std::cout << std::setw(8) << x[i] << std::endl;
	//delete[] x;
	// std::cout << "Jacoby method:" << std::endl;
	// x = jacoby(A, 0);
	// std::cout << "Vector X with cube norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// 	std::cout << std::setw(8) << x[i] << std::endl;
	/*std::cout << "Jacoby method:" << std::endl;
	x = jacoby(A, 1);
	std::cout << "Vector X with octahedral norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	delete[] x;*/

//	system("pause");
=======
	std::cout << "vector b:" << std::endl;
	for (int j = 0; j < A.count; j++)
		std::cout << std::setw(8) << A.rvalue[j] << std::endl;
	std::cout << std::endl;

	//QR
	T.QR_find_x(R);
	if (degenerate_matrix(R) == false)
	{
	 	b_.push_back(multi_vect(R.rvalue, T));
	 	find_x(R, b_, x);
	 	std::cout << "Vector x in QR:" << std::endl;
	 	for (int j = 0; j < A.count; j++)
	  {
	 		std::cout << std::setw(8) << x[j] << std::endl;
	 	}
	 	std::cout << std::endl;

		//find inverse matrix with QR
  	VectorMatrix <double> A_(A.count);
	 	A_.onebyone();
	 	A_.inverse_matrix(R, T);
	 	std::cout << "Inverse matrix with QR method A_:" << std::endl;
	 	A_.print();

		// //find conditional number with different norm
		 std::cout << "Estimation of condition number with cube norme number >= " << cube_estimate_number(A.rvalue, R, T) << std::endl;
		// std::cout << "Conditional number with cube norme = " << cube_norm(A_) * cube_norm(A) << std::endl;
		// std::cout << std::endl;
		// std::cout << "Estimation of condition number with octahedral norme number >= " << octah_estimate_number(A.rvalue, R, T) << std::endl;
		// std::cout << "Conditional number with octahedral norme = " << octah_norm(A_) * octah_norm(A) << std::endl;
		// std::cout << std::endl;
		// std::cout << "Estimation of condition number with octahedral norme number >= " << sfer_estimate_number(A.rvalue, R, T) << std::endl;
		// std::cout << "Conditional number with sferical norme = " << sfer_norm(A_) * sfer_norm(A) << std::endl;
		// std::cout << std::endl;
		//
		// //find residual vector
		// std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.size), A.size) << std::endl;
		// std::cout << std::endl;
		// std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.size), A.size) << std::endl;
		// std::cout << std::endl;
		// std::cout << "Residual vector with sferical norme = " << sfer_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.size), A.size) << std::endl;
		// std::cout << std::endl;

	//}
	//else
	//	std::cout << "Wrong Matrix" << std::endl;

	//delete variable of QR
	// if (x)
	// 	delete[] x;
	// if (b_)
	// 	delete[] b_;

	//Gauss
	GaussRight(A);
	if (degenerate_matrix(A) == false)
	{
		GaussLeft(A);
		//find inverse matrix with Gauss
		VectorMatrix <double> E (A.size);
		E.onebyone();
	//	reversematrixgauss(A_cpy, E);
	//	std::cout << "Inverse matrix with Gauss method A_:" << std::endl;
		E.print();
	}
	else
		std::cout << "Wrong Matrix" << std::endl;

	// system("pause");
>>>>>>> a76e75838e30bcd6d8429795141cb69842d5a2ed
	return (0);
}
