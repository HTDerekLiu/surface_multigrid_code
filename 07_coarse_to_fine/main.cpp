#include <igl/read_triangle_mesh.h>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>

#include <SSP_decimate.h>
#include <single_collapse_data.h>
#include <query_coarse_to_fine.h>

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

	// decimate the input mesh using SSP
	SparseMatrix<double> P;
	int tarF = 1000; // target number of faces
	int dec_type = 0; // decimation type (0:qslim, 1:midpoint, 2:vertex removal)
	VectorXi IM, FIM;
	vector<single_collapse_data> decInfo;
	vector<vector<int>> decIM;
	VectorXi IMF; 
	SSP_decimate(VO,FO,tarF, dec_type, V,F,IMF, IM, decInfo, decIM, FIM);

	// get barycentric coordinates on the coarse mesh for querying (I set it to the coarse vertices in this case)
	MatrixXd BC(V.rows(),3); BC.setZero();
	MatrixXi BF(V.rows(),3); BF.setZero();
	VectorXi FIdx(V.rows()); FIdx.setZero();
	for (int fIdx=0; fIdx<F.rows(); fIdx++)
	{
		for (int ii = 0; ii<F.cols(); ii++)
		{
			int vIdx = F(fIdx,ii);
			if (BC.row(vIdx).sum() == 0.0)
			{
				BC(vIdx,ii) = 1;
				BF.row(vIdx) = F.row(fIdx);
				FIdx(vIdx) = fIdx;	
			}
		}
	}

	// query coarse to fine 
	query_coarse_to_fine(decInfo, IM, decIM, IMF, BC, BF, FIdx);

	// compute point location on the fine mesh
	MatrixXd pt(BC.rows(),3); pt.setZero();
	for (int ii = 0; ii<BC.rows(); ii++)
	{
		pt.row(ii) = BC(ii,0) * VO.row(BF(ii,0)) + BC(ii,1) * VO.row(BF(ii,1)) + BC(ii,2) * VO.row(BF(ii,2));
	}

	// visualize the prolongation operator
	igl::opengl::glfw::Viewer viewer;
	{
		viewer.data().set_mesh(VO, FO);
		// change mesh color
    const Eigen::RowVector3d blue(149.0/255, 217.0/255, 244.0/255);
    viewer.data().set_colors(blue);
		// change background color
    Vector4f backColor;
    backColor << 208/255., 237/255., 227/255., 1.;
    viewer.core().background_color = backColor;
		// plot points
		viewer.data().add_points(pt, Eigen::RowVector3d(0, 0, 0));
		viewer.data().point_size = 10;
		viewer.launch();
	}
	
}
