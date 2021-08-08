#include <Eigen/Core>
#include <vector>
#include <iostream>

#include <igl/min_quad_with_fixed.h>
#include <igl/writeDMAT.h>

#include "MaterialModel.h"
#include "MeshConnectivity.h"
#include "ElasticShell.h"

#include <profc.h>

template <class SFF>
double implicit_euler_balloon(const MeshConnectivity & mesh,
    const Eigen::SparseMatrix<double> & M,
    Eigen::MatrixXd & curPos,
    Eigen::VectorXd & qdot,
    const Eigen::VectorXd & fExt,
    const Eigen::VectorXi & bi,
    const Eigen::VectorXd & curEdgeDOFs,
    const MaterialModel<SFF> & mat,
    const double & dt,
    const Eigen::VectorXd & thicknesses,
    const std::vector<Eigen::Matrix2d> & abars,
    const std::vector<Eigen::Matrix2d> & bbars)
{
    Eigen::VectorXd qdot0 = qdot;
    Eigen::MatrixXd curPos0 = curPos;

    Eigen::VectorXd dx; // update of qdot i.e. dqot

    Eigen::VectorXd bc(bi.rows()); // place holder
    bc.setZero();
    Eigen::VectorXd Beq;
    Eigen::SparseMatrix<double> Aeq;

    //compute newton step direction
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    for (int i = 0; i < 10; i++) {

        Eigen::VectorXd G;
        std::vector<Eigen::Triplet<double>> hessian;

        double energy = ElasticShell<SFF>::elasticEnergy(mesh, curPos, curEdgeDOFs, mat, thicknesses, abars, bbars, ElasticShell<SFF>::EnergyTerm::ET_STRETCHING, &G, &hessian);
        Eigen::SparseMatrix<double> K(qdot.rows(), qdot.rows());
        K.setFromTriplets(hessian.begin(), hessian.end());

        Eigen::SparseMatrix<double> tmp_H;
        Eigen::VectorXd tmp_g;

        // dynamic solve
        // Eigen::SparseMatrix<double> I(K.rows(), K.rows());
        // I.setIdentity();
        tmp_H = M + dt * dt * K;
        tmp_g = M * (qdot - qdot0) + dt * G + dt * fExt;

        {
            PROFC_NODE("direct solve time");
            // min_quad_with_fixed: dirichlet boundary condition
            igl::min_quad_with_fixed_data<double> mqwf;
            igl::min_quad_with_fixed_precompute(tmp_H, bi, Aeq, true, mqwf);
            igl::min_quad_with_fixed_solve(mqwf, tmp_g, bc, Beq, dx);
            std::cout << "euler converge: " << tmp_g.transpose() * dx << std::endl;
        }
        // //check for convergence
        // if(tmp_g.transpose() * dx > -1e-3) {
        //     std::cout << "break: " << i << std::endl;
        //     break;
        // }

        // solver.compute(tmp_H); // precompute
        // if (solver.info() != Eigen::Success) {
        //     std::cout<<"newtons_method: decomposition failed\n";
        //     return std::numeric_limits<double>::infinity();
        // }
        // dx = -solver.solve(tmp_g); // solve
        // if (solver.info() != Eigen::Success) {
        //     std::cout<<"newtons_method: solve failed\n";
        //     return std::numeric_limits<double>::infinity();
        // }


        //perform backtracking on the newton search direction
        //Guarantee that step is a descent step and that it will make sufficient decreate
        double alpha = 1;
        double p = 0.5;
        double c = 1e-8;


        auto f = [&](const Eigen::VectorXd & tmp_qdot) -> double {
            double E_total = 0;
            double Ek = 0.5 * (tmp_qdot - qdot0).transpose() * M * (tmp_qdot - qdot0);
            Eigen::MatrixXd newPos = curPos0;
            for (int j = 0; j < newPos.rows(); j++)
            {
                newPos.row(j) += dt * tmp_qdot.segment<3>(3 * j);
                E_total += newPos.row(j) * fExt.segment<3>(3 * j);
            }
            double V = ElasticShell<SFF>::elasticEnergy(mesh, newPos, curEdgeDOFs, mat, thicknesses, abars, bbars, ElasticShell<SFF>::EnergyTerm::ET_STRETCHING, NULL, NULL);
            E_total = E_total + Ek + V;

            return E_total;
        };

        double f0 = f(qdot);

        double s = f0 + c * tmp_g.transpose() * dx; // sufficient decrease

        while (alpha > 1e-8) {
   
            if (f(qdot + alpha * dx) <= s) {
                qdot += alpha * dx;
                break;
            }
            alpha *= p;
        }
        std::cout << "alpha: " << alpha << std::endl;

        // update curPos
        curPos = curPos0;
        for (int j = 0; j < curPos.rows(); j++) {
            curPos.row(j) += dt * qdot.segment<3>(3 * j);
        }

        // // a simple collision handler
        // for (int j = 0; j < curPos.rows(); j++) {
        //     if (curPos(j,1) < 0) {
        //         curPos(j,1) = 0;
        //         const double cr = 10;
        //         qdot(3*j+1) = -qdot(3*j+1) / cr;
        //     }
        // }

    }

    return 0;

}