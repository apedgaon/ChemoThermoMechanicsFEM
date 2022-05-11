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
    void solve2();
    void solve3();
    void post_print(std::string file);

    Geometry<dim> geom;
    Material mat;
    Mesh2<dim> mesh;
    std::vector<DirichletBC> dirBCs;
    std::vector<DirichletBC> dirBCsT;
    std::vector<DirichletBC> dirBCsc;
    int phytype;    // Physics Type
    SolutionScheme sol_scheme;
    unsigned int nsteps;
    double end_time;

    QuadratureType quadrature_type;
    unsigned int Nq;
    Eigen::MatrixXd Qpts;
    Eigen::VectorXd wts;
    Eigen::MatrixXd N;
    Eigen::MatrixXd dNdxi[dim];



private:
    void init_shp_functions();
    void assemble();
    void assemble2();
    void assemble3();
    void compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx, Eigen::MatrixXd& xe, Eigen::VectorXd& de);
    void compute_elem_quantitiesT(connectivity& conn, Eigen::VectorXui& gidxT, Eigen::VectorXd& Te);
    void compute_elem_quantitiesc(connectivity& conn, Eigen::VectorXui& gidxc, Eigen::VectorXd& ce);
    void compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::MatrixXd& jac,
        Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac);
    void compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kl);
    void compute_quad_stiffTT(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& KTTl);
    void compute_quad_stiffuT(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& KuTl);
    void compute_quad_massTT(unsigned iq, double det_jac, Eigen::MatrixXd& MTTl);
    void compute_quad_stiffcc(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kccl);
    void compute_quad_stiffuc(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kucl);
    void compute_quad_stiffcT(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& KcTl, Eigen::MatrixXd& KTcl);
    void apply_dirBC();
    void create_dirichlet_map();
    void reform_full_sol();
    void raw_sol_print();

private:
    // Mechanics Matrices
    Eigen::MatrixXd Kuu;
    Eigen::VectorXd Fu;
    Eigen::VectorXd du;
    Eigen::MatrixXd KuT;
    Eigen::MatrixXd Kuc;

    Eigen::MatrixXd Kuud;
    Eigen::VectorXd Fud;
    Eigen::VectorXd dud;
    Eigen::MatrixXd KuTd;
    Eigen::MatrixXd Kucd;

    unsigned int gdofsu;
    unsigned int ucdofsu;
    std::vector<unsigned int> dirmapu;

    // Thermal Transport Matrices
    Eigen::MatrixXd KTT;
    Eigen::MatrixXd KTc;
    Eigen::VectorXd FT;
    Eigen::VectorXd dT;
    Eigen::MatrixXd MTT;

    Eigen::MatrixXd KTTd;
    Eigen::MatrixXd KTcd;
    Eigen::VectorXd FTd;
    Eigen::VectorXd dTd;
    Eigen::MatrixXd MTTd;

    unsigned int gdofsT;
    unsigned int ucdofsT;
    std::vector<unsigned int> dirmapT;

    // Chemical Transport Matrices
    Eigen::MatrixXd Kcc;
    Eigen::MatrixXd KcT;
    Eigen::VectorXd Fc;
    Eigen::VectorXd dc;
    Eigen::MatrixXd Mcc;

    Eigen::MatrixXd Kccd;
    Eigen::MatrixXd KcTd;
    Eigen::VectorXd Fcd;
    Eigen::VectorXd dcd;
    Eigen::MatrixXd Mccd;

    unsigned int gdofsc;
    unsigned int ucdofsc;
    std::vector<unsigned int> dirmapc;
};