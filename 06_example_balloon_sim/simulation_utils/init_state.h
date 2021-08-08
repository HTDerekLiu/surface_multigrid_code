#include <Eigen/Core>

void init_state(const Eigen::MatrixXd & V,
                Eigen::VectorXd & q,
                Eigen::VectorXd & qdot) {

    q.resize(V.rows()*V.cols());
    qdot.resize(V.rows()*V.cols());

    Eigen::MatrixXd Vt = V.transpose();
    q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows()*Vt.cols());
    qdot.setZero();

}
