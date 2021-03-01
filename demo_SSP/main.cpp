#include <igl/read_triangle_mesh.h>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>

#include <single_collapse_data.h>
#include <query_fine_to_coarse.h>
#include <get_prolong.h>

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

	// load mesh
	MatrixXd VO,V;
	MatrixXi FO,F;
	{
		igl::read_triangle_mesh("../../meshes/ogre_sim.obj", VO, FO);
		cout << "original mesh: |V| " << VO.rows() << ", |F|: " << FO.rows() << endl;
	}

	// decimate 
	SparseMatrix<double> P;
	int tarF = 1000; // target number of faces
	int dec_type = 0; // decimation type (0:qslim, 1:midpoint, 2:vertex removal)
	get_prolong(VO,FO,tarF,dec_type,V,F,P);

	MatrixXd pt = P * V;

	// Plot the results
	igl::opengl::glfw::Viewer viewer;
	{
		viewer.data().set_mesh(V, F);
			viewer.data().add_points(pt, Eigen::RowVector3d(0, 0, 0));
		viewer.data().point_size = 10;
		viewer.launch();
	}
	
}
