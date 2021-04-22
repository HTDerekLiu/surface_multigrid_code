#include <Eigen/Core>
#include <Eigen/Sparse>

// Sparse Version LBS
void lumped_mass_matrix(const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::SparseMatrix<double> & M);