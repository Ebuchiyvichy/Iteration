
#include "Zeldel.h"

int	main()
{
	VectorMatrix <mytype>	A(202);
	for (int i = 1; i != A.count-1; i++)
	{
		A.value[i][i+1] = 1;
		A.value[i][i] = 4;
		A.value[i][i-1] = 1;
	}
	A.value[0][0] = 4;
	A.value[0][1] = 1;
	A.value[201][200] = 1;
	A.value[201][201] = 4;
	for (int i = 0; i != A.count; i++)
	{
		//rvalue.push_back();
		if (i == 0)
			A.rvalue.push_back(6);
		else if (i == 201)
			A.rvalue.push_back(9 - 3 * (202 % 2));
		else
			A.rvalue.push_back(10 - 2 * (i % 2));
	}
   // A.VectorMatrixC(202);
   //std::cout << "\n";
	std::vector<mytype> x;
	// mytype w = searchbigomegaiter(A);
	//smallzeidel
	//std::cout << "Zeldel for small one method:" << std::endl;
	//x = SmallRelax(A, 0);
	//std::cout << "Vector X with cube norme:" << std::endl;
	//for (int i = 0; i < 202; i++)
	//std::cout << i << "  " <<std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
  // x = SmallZeydel(A, 1);
	// std::cout << "Vector X with octahedral norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// 	std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	//

//	bigZeidel
//std::cout << "Zeldel for big one method:" << std::endl;
    x = SmallRelax(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
	std::cout << std::setw(8) << x[i] << std::endl;

	//relax
 	//std::cout << "Relax for big one method:" << std::endl;

	// std::cout << "Relax for small one method:" << std::endl;
	// x = SmallRelax(B, 0);
	// std::cout << "Vector X with cube norme:" << std::endl;
	// for (int i = 0; i < B.count; i++)
	// std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
  // x = SmallRelax(A, 1);
	// std::cout << "Vector X with octahedral norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// 	std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	//
	// //simple_iter
	// std::cout << "Simple iteration method:" << std::endl;
	// x = simp_iter(A, 0);
	// std::cout << "Vector X with cube norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// 	std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// x = simp_iter(A, 1);
	// std::cout << "Vector X with octahedral norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// 	std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	//
	// //jacoby
	// std::cout << "Jacoby method:" << std::endl;
	// x = jacoby(A, 0);
	// std::cout << "Vector X with cube norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	// std::cout << "Jacoby method:" << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// x = jacoby(A, 1);
	// std::cout << "Vector X with octahedral norme:" << std::endl;
	// for (int i = 0; i < A.count; i++)
	// std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;

	//system("pause");

	return (0);
}
