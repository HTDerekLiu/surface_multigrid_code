#include <igl/read_triangle_mesh.h>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>

#include <mg_data.h>
#include <mg_precompute.h>

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

	// load mesh
	MatrixXd VO,V;
	MatrixXi FO,F;
	{
		igl::read_triangle_mesh("../../meshes/bunny.obj", VO, FO);
		cout << "original mesh: |V| " << VO.rows() << ", |F|: " << FO.rows() << endl;
	}

	// construct multigrid hierarchy
	int min_coarsest_nV = 500;
	float coarsening_ratio = 0.25;
	int decimation_type = 1;
	vector<mg_data> mg;
	mg_precompute(VO,FO,coarsening_ratio, min_coarsest_nV, decimation_type, mg);

	//Viewer that shows all functions: zexact, znoisy, zl, zh
	MatrixXd pt;
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(mg[0].V, mg[0].F);
	viewer.callback_key_down =
	[&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
	{
			switch(key) {
			case '0':
				viewer.data().clear();
				viewer.data().set_mesh(mg[0].V,mg[0].F);
				break;
			case '1':
				pt = mg[1].P * mg[1].V;
				viewer.data().clear();
				viewer.data().set_mesh(mg[1].V,mg[1].F);
				viewer.data().add_points(pt, Eigen::RowVector3d(0, 0, 0));
				viewer.data().point_size = 7;
				break;
			case '2':
				pt = mg[2].P * mg[2].V;
				viewer.data().clear();
				viewer.data().set_mesh(mg[2].V,mg[2].F);
				viewer.data().add_points(pt, Eigen::RowVector3d(0, 0, 0));
				viewer.data().point_size = 10;
				break;
			default:
				return false;
			}
			return true;
	};
	viewer.launch();
	
}
