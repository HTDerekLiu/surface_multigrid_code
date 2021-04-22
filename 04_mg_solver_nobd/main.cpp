#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/boundary_loop.h>
#include <igl/massmatrix.h>
#include <igl/writeOBJ.h>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>
#include <stdlib.h>

#include <mg_data.h>
#include <mg_precompute.h>
#include <normalize_unit_area.h>
#include <min_quad_with_fixed_mg.h>

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

	MatrixXd V;
	MatrixXi F;
	VectorXi b;
	{
		// load mesh
		igl::read_triangle_mesh("../../meshes/hilbert_cube.obj", V, F);

		// load constrained points
		MatrixXd V_known;
		MatrixXi F_known;
		igl::read_triangle_mesh("../../meshes/hilbert_cube_known.obj", V_known, F_known);

		// get constrained vertex index (TODO: make it faster by using spatial hash)
		b.resize(V_known.rows());
		for (int ii=0; ii<V_known.rows(); ii++)
		{
			double minD = std::numeric_limits<double>::max();
			int minIdx = -1;
			RowVector3d v = V_known.row(ii);
			for (int vIdx=0; vIdx<V.rows(); vIdx++)
			{
				double distance = (V.row(vIdx) - v).norm();
				if (distance < minD)
				{
					minD = distance;
					minIdx = vIdx;
				}
			}
			b(ii) = minIdx;
		}

		// normalize mesh area
		normalize_unit_area(V,F);
		cout << "original mesh: |V| " << V.rows() << ", |F|: " << F.rows() << endl;
	}

	// construct multigrid hierarchy
	int min_coarsest_nV = 500;
	float coarsening_ratio = 0.25;
	int decimation_type = 1;
	vector<mg_data> mg;
	mg_precompute(V,F,coarsening_ratio, min_coarsest_nV, decimation_type, mg);

	// construct a toy Poisson problem 
	// solve A*z = B, s.t. z(b) = bval
	// where b are the boundary vertices (thus this demo is for meshes with boundaries)
	SparseMatrix<double> A;
	igl::cotmatrix(V,F,A);
	A = -A;

	// fixed boundary value
	VectorXd bval(b.size());
	bval.setZero();
	
	// get RHS
	SparseMatrix<double> M;
	igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);

	VectorXd ones(V.rows());
	ones.setOnes();
	VectorXd B = M * ones;
	for (int ii = 0; ii<b.size(); ii++)
		B(b(ii)) = bval(ii);

	// initial guess z0
	srand((unsigned int) time(0));
	VectorXd z0(V.rows());
	z0.setRandom();

	// get some solver data
	VectorXd z = z0;
	min_quad_with_fixed_mg_data solverData;
	SimplicialLDLT<SparseMatrix<double>> coarseSolver;
	min_quad_with_fixed_mg_precompute(A,b,solverData, mg, coarseSolver);

	// multigrid solve
	vector<double> rHis;
	// min_quad_with_fixed_mg_solve(solverData, B, bval, z0, coarseSolver, mg, z, rHis);
	min_quad_with_fixed_mg_solve(solverData, B, bval, z0, coarseSolver, 1e-10, mg, z, rHis);

	// Plot the results
	igl::opengl::glfw::Viewer viewer;
	{
		viewer.data().set_mesh(V, F);

		// set mesh color	
		viewer.data().set_data(z);

		// set background color
		Vector4f backColor;
		backColor << 208/255., 237/255., 227/255., 1.;
		viewer.core().background_color = backColor;

		// set edge color
		Vector4f edgeColor;
		edgeColor << .2,.2,.2, .8;
		viewer.data().line_color = edgeColor;
		viewer.launch();
	}
}
