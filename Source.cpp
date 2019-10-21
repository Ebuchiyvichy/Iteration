#include "Jacoby.h"

int	main()
{
//	std::ifstream file;
//	file.open("matrixexample.txt");
//	int	size;
//	file >> size;
	VectorMatrix <mytype>	A(4);
	mytype					*x;

	A.init();
	std::cout << "Matrix A:" << std::endl;
	A.print();
	std::cout << "Simple iteration method:" << std::endl;
	x = simp_iter(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(x, multi_vect(A.rvalue, A), A.count), A.count) << std::endl;
	std::cout << std::endl;
	delete[] x;
	//x = simp_iter(A, 1);
	//std::cout << "Vector X with octahedral norme:" << std::endl;
	//for (int i = 0; i < A.count; i++)
	//	std::cout << std::setw(8) << x[i] << std::endl;
	//delete[] x;
	std::cout << "Jacoby method:" << std::endl;
	x = jacoby(A, 0);
	std::cout << "Vector X with cube norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	/*std::cout << "Jacoby method:" << std::endl;
	x = jacoby(A, 1);
	std::cout << "Vector X with octahedral norme:" << std::endl;
	for (int i = 0; i < A.count; i++)
		std::cout << std::setw(8) << x[i] << std::endl;
	delete[] x;*/

	system("pause");
	return (0);
}