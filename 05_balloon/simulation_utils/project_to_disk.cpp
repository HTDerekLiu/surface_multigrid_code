#include "project_to_disk.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>

using namespace Eigen;

void project_to_disk(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    Eigen::VectorXi b;
    igl::boundary_loop(F,b);
    Eigen::MatrixXd bc;
    igl::map_vertices_to_circle(V,b,bc);
    
    Eigen::MatrixXd UV;
    igl::harmonic(F,b,bc,1,UV);
    UV.col(1) *= -1.;
    U.resizeLike(V);
    // U << UV.col(0),0.02*Eigen::VectorXd::Random(V.rows()),UV.col(1);
    U << UV.col(0),Eigen::VectorXd::Zero(V.rows()),UV.col(1);
}