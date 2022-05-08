#pragma once
#include "FEMUtility.h"


enum PhysicsType : int
{
    _Mechanics = 1,
    _Thermal = 2,
    _Chemical = 4
};


enum class SolutionScheme
{
    monolithic = 0,
    staggered = 1
};

template <unsigned int dim>
class CoupledFieldFEM
{
public:
    void generate_mesh();
    void initialize();
    void solve();
    void post_print(std::string file);

    Geometry<dim> geom;
    Material mat;
    Mesh2<dim> mesh;
    std::vector<DirichletBC> dirBCs;
    int phytype;    // Physics Type
    SolutionScheme sol_scheme;

    QuadratureType quadrature_type;
    unsigned int Nq;
    Eigen::MatrixXd Qpts;
    Eigen::VectorXd wts;
    Eigen::MatrixXd N;
    Eigen::MatrixXd dNdxi[dim];



private:
    void init_shp_functions();
    void assemble();
    void compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx, Eigen::MatrixXd& xe, Eigen::VectorXd& de);
    void compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::MatrixXd& jac,
        Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac);
    void compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kl);
    void apply_dirBC();
    void create_dirichlet_map();
    void reform_full_sol();
    void raw_sol_print();

private:
    Eigen::MatrixXd Ku;
    Eigen::VectorXd Fu;
    Eigen::VectorXd du;

    Eigen::MatrixXd Kud;
    Eigen::VectorXd Fud;
    Eigen::VectorXd dud;

    unsigned int gdofs;
    unsigned int ucdofs;

    std::vector<unsigned int> dirmap;
};