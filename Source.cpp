
#include "Zeldel.h"

int	main()
{

	VectorMatrix <mytype>	A(4);
	A.init();
	A.print();
	std::vector<mytype> x;
	//smallzeidel
	std::cout << "Zeldel for small one method:" << std::endl;
	x = SmallZeydel(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
	std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
/*	x = SmallZeydel(A, 1);
	std::cout << "Vector X with octahedral norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;*/

	
	//bigZeidel
/*	std::cout << "Zeldel for big one method:" << std::endl;
	x = BigZeydel(A);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
	std::cout << std::setw(8) << x[i] << std::endl;*/

	//relax
	std::cout << std::endl;
	std::cout << "Relax for small one method:" << std::endl;
	x = SmallRelax(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
	std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
/*	x = SmallRelax(A, 1);
	std::cout << "Vector X with octahedral norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;*/

	//simple_iter
/*	std::cout << "Simple iteration method:" << std::endl;
	x = simp_iter(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	x = simp_iter(A, 1);
	std::cout << "Vector X with octahedral norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;*/

	//jacoby
/*	std::cout << "Jacoby method:" << std::endl;
	x = jacoby(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
	std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.count), A.count) << std::endl;
	std::cout << "Jacoby method:" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	x = jacoby(A, 1);
	std::cout << "Vector X with octahedral norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
	std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;*/
	
	system("pause");

	return (0);
}