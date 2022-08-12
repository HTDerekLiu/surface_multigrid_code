#include <igl/read_triangle_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/writeOBJ.h>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>
#include <string> 
#include <stdexcept>

#include <SSP_random_qslim.h>
#include <single_collapse_data.h>
#include <query_coarse_to_fine.h>

int find_row_with_elements(
	const Eigen::MatrixXi & F,
	const std::vector<int> & adjV)
{
	for (int r=0; r<F.rows(); r++)
	{
		int count = 0;
		for (int c=0; c<F.cols(); c++)
		{
			int vIdx = F(r,c);
			for (int ii=0; ii<adjV.size(); ii++)
			{
				if (vIdx == adjV[ii])
				{
					count += 1;
					break;
				}
			}
		}
		if (count == adjV.size())
			return r;
	}
	return -1;
}

void loop_upsample_connecitivy( 
	const int & n_verts, 
	const Eigen::MatrixXi & F,
	Eigen::SparseMatrix<double> & S,
	Eigen::MatrixXi & NF)
{
	// This function is the modified version of igl::upsample. This upsampling stacks faces differently from igl::upsample so that it is consistent with the neural subdivision code
	using namespace std;
	using namespace Eigen;

	MatrixXi FF, FFi;
	igl::triangle_triangle_adjacency(F,FF,FFi);

	// Compute the number and positions of the vertices to insert (on edges)
	Eigen::MatrixXi NI = Eigen::MatrixXi::Constant(FF.rows(),FF.cols(),-1);
	Eigen::MatrixXi NIdoubles = Eigen::MatrixXi::Zero(FF.rows(), FF.cols());
	int counter = 0;

	// Compute the number and positions of the vertices to insert (on edges)
	for(int i=0;i<FF.rows();++i)
	{
		for(int j=0;j<3;++j)
		{
			if(NI(i,j) == -1)
			{
				NI(i,j) = counter;
				NIdoubles(i,j) = 0;
				if (FF(i,j) != -1) {
					//If it is not a boundary
					NI(FF(i,j), FFi(i,j)) = counter;
					NIdoubles(i,j) = 1;
				}
				++counter;
			}
		}
	}

	const int& n_odd = n_verts;
	const int& n_even = counter;
	const int n_newverts = n_odd + n_even;

	// Construct the subdivision operator S such that NV = S * V
	vector<Triplet<double>> tripletList;

	// Fill the odd vertices position
	for (int i=0; i<n_odd; ++i)
		tripletList.emplace_back(i, i, 1.);
	// Fill even vertices position
	for(int i=0;i<FF.rows();++i)
	{
		for(int j=0;j<3;++j)
		{
			if(NIdoubles(i,j)==0) {
				tripletList.emplace_back(NI(i,j) + n_odd, F(i,j), 1./2.);
				tripletList.emplace_back(NI(i,j) + n_odd, F(i,(j+1)%3), 1./2.);
			}
		}
	}
	S.resize(n_newverts, n_verts);
  	S.setFromTriplets(tripletList.begin(), tripletList.end());

	// Build the new face list (see https://github.com/HTDerekLiu/neuralSubdiv/blob/master/utils_matlab/midPointUpsample.m to understand the notations) 
	NF.resize(F.rows()*4,3);
	int fIdx = 0;
	// put (f0, i2, i1)
	for(int i=0; i<F.rows();++i)
	{
		int f0 = F(i,0);
		int i2 = NI(i,0) + n_odd;
		int i1 = NI(i,2) + n_odd;
		NF.row(fIdx) << f0, i2, i1;
		fIdx += 1;
	}
	// put (f1, i0, i2)
	for(int i=0; i<F.rows();++i)
	{
		int f1 = F(i,1);
		int i0 = NI(i,1) + n_odd;
		int i2 = NI(i,0) + n_odd;
		NF.row(fIdx) << f1, i0, i2;
		fIdx += 1;
	}
	// put (f2, i1, i0)
	for(int i=0; i<F.rows();++i)
	{
		int f2 = F(i,2);
		int i1 = NI(i,2) + n_odd;
		int i0 = NI(i,1) + n_odd;
		NF.row(fIdx) << f2, i1, i0;
		fIdx += 1;
	}
	// put (f2, i1, i0)
	for(int i=0; i<F.rows();++i)
	{
		int i0 = NI(i,1) + n_odd;
		int i1 = NI(i,2) + n_odd;
		int i2 = NI(i,0) + n_odd;
		NF.row(fIdx) << i0, i1, i2;
		fIdx += 1;
	}
}

void loop_upsample_connecitivy(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF,
  const int number_of_subdivs)
{
	// This function is the modified version of igl::upsample. This upsampling stacks faces differently from igl::upsample so that it is consistent with the neural subdivision code
	using namespace std;
	using namespace Eigen;

	NV = V;
	NF = F;
	for(int i=0; i<number_of_subdivs; ++i)
	{
		MatrixXi tempF = NF;
		SparseMatrix<double> S;
		loop_upsample_connecitivy(NV.rows(), tempF, S, NF);
		NV = (S*NV).eval(); // This .eval is super important
	}
}

void loop_upsample_barycentric(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const int num_subdiv_iters,
	Eigen::MatrixXd & BC,
	Eigen::MatrixXi & BF,
	Eigen::VectorXi & FIdx,
	Eigen::MatrixXi & NF)
{
	using namespace std;
	using namespace Eigen;
	
	// construct mid point subdivision operator S such that V_new = S*V
	NF = F; // new face list
	Eigen::SparseMatrix<double> S; // mid point subdivision operator
	{
		for(int ii=0; ii<num_subdiv_iters; ++ii)
		{
			MatrixXi tempF = NF;
			if (ii == 0) // first iteration
			{
				loop_upsample_connecitivy(V.rows(), tempF, S, NF);
			}
			else // further iterations
			{
				Eigen::SparseMatrix<double> SS;
				loop_upsample_connecitivy(S.rows(), tempF, SS, NF);
				S = SS * S;
			}
		}
	}

	// extract barycentric coordinates from S
	Eigen::SparseMatrix<double,Eigen::RowMajor> SR = S;
	int nV = S.rows();
	BC.resize(nV,3); BC.setZero();
	BF.resize(nV,3); BF.setZero();
	FIdx.resize(nV); FIdx.setZero();
	for (int vIdx=0; vIdx<SR.outerSize(); ++vIdx)
	{
		vector<int> adjV; adjV.reserve(3); // adjacent vertex indices
		vector<double> b; b.reserve(3); // barycentric coordinates

		for (SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(SR,vIdx); it; ++it)
		{
			adjV.push_back(it.col());
			b.push_back(it.value());
		}

		// find which face has this vertex (this is a slow implementation)
		int fIdx = find_row_with_elements(F,adjV);
		assert(fIdx >= 0); // ensure we find a row with all the adjV
		FIdx(vIdx) = fIdx;
		BF.row(vIdx) = F.row(fIdx);

		// put the information back to barycentric coordinates
		for (int ii=0; ii<3; ii++)
		{
			for (int jj=0; jj<adjV.size(); jj++)
			{
				if ( F(fIdx,ii) == adjV[jj])
				{
					BC(vIdx,ii) = b[jj];
					break;
				}
			}
		}
	}
}

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

	// parsing arguments ./random_subdiv_remesh_bin [mesh_path] [target_faces] [number_subdivision] [random_seed]
	string mesh_path;
	int tarF, num_subdivs, random_seed;
	if (argc == 5)
	{
		mesh_path = argv[1]; // path to mesh
		tarF = stoi(argv[2]); // number of faces on the coarse mesh
		num_subdivs = stoi(argv[3]); // number of subdivisions for upsampling
		random_seed = stoi(argv[4]); // random seed
	}
	else
	{
		throw std::invalid_argument( "invalid input arguments. It should be ./random_subdiv_remesh_bin [mesh_path] [target_faces] [number_subdivision] [random_seed]" );
	}

	// load mesh
	MatrixXd VO;
	MatrixXi FO;
	igl::read_triangle_mesh(mesh_path, VO, FO);
	cout << "original mesh: |V| " << VO.rows() << ", |F|: " << FO.rows() << endl;

	// decimate the input mesh using SSP
	MatrixXd V; // coarse vertices
	MatrixXi F; // coarse faces
	SparseMatrix<double> P;
	VectorXi IM, FIM;
	vector<single_collapse_data> decInfo;
	vector<vector<int>> decIM;
	VectorXi IMF; 
	// SSP_decimate(VO,FO,tarF, dec_type, V,F,IMF, IM, decInfo, decIM, FIM);
	srand(random_seed);
	SSP_random_qslim(VO,FO,tarF,V,F,IMF,IM,decInfo,decIM,FIM);

	// loop upsample the mesh 
	MatrixXd BC; // barycentric coordinates 
	MatrixXi BF; // barycentric faces 
	VectorXi FIdx; // indices of barycentric faces
	MatrixXi SF; // subdivided faces
	loop_upsample_barycentric(V,F,num_subdivs,BC,BF,FIdx,SF);
	query_coarse_to_fine(decInfo, IM, decIM, IMF, BC, BF, FIdx);

	// compute subdivided vertex locations 
	MatrixXd SV; // subdivided vertices
	SV.resize(BC.rows(),3); SV.setZero();
	for (int ii = 0; ii<BC.rows(); ii++)
	{
		SV.row(ii) = BC(ii,0) * VO.row(BF(ii,0)) + BC(ii,1) * VO.row(BF(ii,1)) + BC(ii,2) * VO.row(BF(ii,2));
	}

	// split the subdivided meshes into levels
	for (int iter=0; iter<=num_subdivs; iter++)
	{
		MatrixXd NV;
		MatrixXi NF;
		loop_upsample_connecitivy(V,F,NV,NF,iter);
		int nV = NV.rows();
		NV = SV.block(0,0,nV,3);
		std::string output_name = "../output_s" + std::to_string(iter) + ".obj";
		igl::writeOBJ(output_name, NV,NF);
	}
}
