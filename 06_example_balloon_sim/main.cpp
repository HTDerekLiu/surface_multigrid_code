#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>

// simulation
#include <MeshConnectivity.h>
#include <ElasticShell.h>
#include <MidedgeAngleTanFormulation.h>
#include <MidedgeAngleSinFormulation.h>
#include <MidedgeAverageFormulation.h>
#include <StVKMaterial.h>
#include <TensionFieldStVKMaterial.h>
#include <NeoHookeanMaterial.h>
#include <init_state.h>
#include <lumped_mass_matrix.h>
#include <implicit_euler_mg_balloon.h>
#include <implicit_euler_balloon.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <vector>
#include <stdlib.h>

#include <mg_data.h>
#include <mg_precompute_block.h>
#include <normalize_unit_area.h>
#include <min_quad_with_fixed_mg.h>

std::vector<mg_data> mg;

int numSteps;
double thickness;
double mg_tolerance = 2e-1;
double poisson;
double dt = 0.001;
int frame_num = 1000;
int matid;
int sffid;
bool useMG = true;

Eigen::MatrixXd curPos, curPosT, origV;
Eigen::VectorXd q, qdot, fExt;
Eigen::MatrixXi F, origF;
MeshConnectivity mesh;

Eigen::SparseMatrix<double> M, Mv;
Eigen::VectorXd Mvd;

using namespace std;

int counter = 0;
double young = 6e6;

void lameParameters(double & young, double &alpha, double &beta)
{
    alpha = young * poisson / (1.0 - poisson * poisson);
    beta = young / 2.0 / (1.0 + poisson);
}

template <class SFF>
void runSimulation(igl::opengl::glfw::Viewer &viewer, 
    const MeshConnectivity &mesh, 
    Eigen::MatrixXd &curPos, 
    Eigen::VectorXd &fExt,
    const Eigen::VectorXi &bi,
    const Eigen::VectorXd &thicknesses,
    double lameAlpha,
    double lameBeta,
    int matid)
{
    // initialize default edge DOFs (edge director angles)
    Eigen::VectorXd edgeDOFs;
    SFF::initializeExtraDOFs(edgeDOFs, mesh, origV);
    // initialize first fundamental forms to those of input mesh
    std::vector<Eigen::Matrix2d> abar;
    ElasticShell<SFF>::firstFundamentalForms(mesh, origV, abar);

    // // initialize second fundamental forms to those of input mesh
    std::vector<Eigen::Matrix2d> bbar;
    ElasticShell<SFF>::secondFundamentalForms(mesh, origV, edgeDOFs, bbar);

    MaterialModel<SFF> *mat;
    switch (matid)
    {
    case 0:
        mat = new NeoHookeanMaterial<SFF>(lameAlpha, lameBeta);
        break;
    case 1:
        mat = new StVKMaterial<SFF>(lameAlpha, lameBeta);
        break;
    case 2:
        mat = new TensionFieldStVKMaterial<SFF>(lameAlpha, lameBeta);
        break;
    default:
        assert(false);
    }

    // newton steps
    double reg = 1e-6;
    for (int j = 0; j <= numSteps; j++)
    {
        std::cout << "iter: " << j << std::endl;

        // compute external force
        Eigen::MatrixXd Nv;
        igl::per_vertex_normals(curPos,mesh.faces(),Nv);
        const double avg = igl::avg_edge_length(curPos,mesh.faces());
        igl::massmatrix(curPos,mesh.faces(),igl::MASSMATRIX_TYPE_DEFAULT,Mv);
        Mvd = Mv.diagonal();
        for (int vi = 0; vi < Nv.rows(); vi++) {
            Eigen::RowVector3d fe = Nv.row(vi) * Mvd(vi);
            fExt.segment<3>(3 * vi) = -fe.transpose() * 1000000;
        }

        // time stepping
        if (useMG)
        {
            PROFC_NODE("implicit euler time");
            implicit_euler_mg_balloon(mesh, M, curPos, qdot, fExt, bi, edgeDOFs, *mat, dt, thicknesses, abar, bbar, mg, mg_tolerance);
        }
        else
        {
            PROFC_NODE("implicit euler time");
            implicit_euler_balloon(mesh, M, curPos, qdot, fExt, bi, edgeDOFs, *mat, dt, thicknesses, abar, bbar);
        }

        // // save the meshes
        // {
        //     string name = "output";
        //     name.append(6-to_string(counter).length(), '0');
        //     if (useMG)
        //         name = "../mg_results/" + name + to_string(counter) + ".obj";
        //     else
        //         name = "../direct_results/" + name + to_string(counter) + ".obj";
        //     igl::writeOBJ(name, curPos, mesh.faces());
        //     counter++;
        // }
    }

    viewer.data().set_vertices(curPos);
    viewer.data().compute_normals();
    // delete mat;
}

int main(int argc, char *argv[])
{    
    using namespace std;
    using namespace Eigen;

    numSteps = 1;

    // set up material parameters
    thickness = 1e-1; 
    poisson = 0.5;
    matid = 0;
    sffid = 2;

    igl::readOBJ("../../meshes/bunny_15K_init.obj", origV, F);
    // igl::readOBJ("../../meshes/bunny_140K_init.obj", origV, F);
    
    // precompute multigrid hierarchy 
    mg_precompute_block(origV,F,mg);
     
    // set up mesh connectivity
    mesh = MeshConnectivity(F);
    origF = mesh.faces();

    // initial position
    curPos = origV;
    init_state(curPos, q, qdot); // initial position and velocity
    lumped_mass_matrix(origV, origF, M);
    M = 1000 * M;

    // external forces
    Eigen::VectorXd gravity;
    {
        Eigen::Vector3d g = Eigen::Vector3d(0., -900.8, 0.);
        Eigen::VectorXd gCol = g.replicate(M.rows()/3, 1);
        gravity = -M*gCol;
    }
    // fExt = gravity;
    fExt.resizeLike(gravity);
    fExt.setZero();

    // fixed point
    Eigen::VectorXi bi;

    // libigl viewer
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
   
    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        if (ImGui::Button("Reset", ImVec2(-1, 0)))
        {
            curPos = origV;
            viewer.data().set_vertices(curPos);
        }

        if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDouble("Thickness", &thickness);
            ImGui::InputDouble("Poisson's Ration", &poisson);
            ImGui::Combo("Material Model", &matid, "NeoHookean\0StVK\0\0");
            ImGui::Combo("Second Fundamental Form", &sffid, "TanTheta\0SinTheta\0Average\0\0");
            ImGui::Checkbox("use multigrid ", &useMG);
            ImGui::InputDouble("mg tolerance", &mg_tolerance);
        }

        
        if (ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen))
        {            
            ImGui::InputInt("Num Steps", &numSteps);
            if (ImGui::Button("Optimize Some Step", ImVec2(-1,0)))
            {
                Eigen::VectorXd thicknesses(mesh.nFaces());
                thicknesses.setConstant(thickness);
                double lameAlpha, lameBeta;
                lameParameters(young, lameAlpha, lameBeta);

                switch (sffid)
                {
                case 0:
                    runSimulation<MidedgeAngleTanFormulation>(viewer, mesh, curPos, fExt, bi, thicknesses, lameAlpha, lameBeta, matid);
                    break;
                case 1:
                    runSimulation<MidedgeAngleSinFormulation>(viewer, mesh, curPos, fExt, bi, thicknesses, lameAlpha, lameBeta, matid);
                    break;
                case 2:
                    runSimulation<MidedgeAverageFormulation>(viewer, mesh, curPos, fExt, bi, thicknesses, lameAlpha, lameBeta, matid);
                    break;
                default:
                    assert(false);
                }
            }
        }
    };

    viewer.data().set_face_based(false);
    viewer.data().set_mesh(curPos, mesh.faces()); 
    const Eigen::RowVector3d blue(149.0/255, 217.0/255, 244.0/255);
    viewer.data().set_colors(blue);
    Vector4f backColor;
    backColor << 208/255., 237/255., 227/255., 1.;
    viewer.core().background_color = backColor;
    viewer.launch();

}
