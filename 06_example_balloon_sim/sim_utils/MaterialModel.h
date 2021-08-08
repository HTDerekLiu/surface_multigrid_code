#ifndef MATERIALMODEL_H
#define MATERIALMODEL_H

#include <Eigen/Core>

class MeshConnectivity;

template <class SFF>
class MaterialModel
{
public:
    virtual double stretchingEnergy(
        const MeshConnectivity &mesh,
        const Eigen::MatrixXd &curPos,
        double thickness,
        const Eigen::Matrix2d &abar,
        int face,
        Eigen::Matrix<double, 1, 9> *derivative, // F(face, i)
        Eigen::Matrix<double, 9, 9> *hessian) const = 0;

    virtual double bendingEnergy(
        const MeshConnectivity &mesh,
        const Eigen::MatrixXd &curPos,
        const Eigen::VectorXd &extraDOFs,
        double thickness,
        const Eigen::Matrix2d &abar, const Eigen::Matrix2d &bbar,
        int face,
        Eigen::Matrix<double, 1, 18 + 3*SFF::numExtraDOFs> *derivative, // F(face, i), then the three vertices opposite F(face,i), then the extra DOFs on oppositeEdge(face,i)
        Eigen::Matrix<double, 18 + 3*SFF::numExtraDOFs, 18 + 3*SFF::numExtraDOFs> *hessian) const = 0;
};

#endif