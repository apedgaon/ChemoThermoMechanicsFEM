#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "CoupledFieldFEM.h"

template<unsigned int dim>
void CoupledFieldFEM<dim>::generate_mesh()
{
    mesh.generate(geom);
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::initialize()
{
    // Initialize shape functions
    init_shp_functions();

    // Initialize global vectors
    gdofsu = mesh.Nnd * dim;
    Kuu = Eigen::MatrixXd::Zero(gdofsu, gdofsu);
    Fu = Eigen::VectorXd::Zero(gdofsu);
    du = Eigen::VectorXd::Zero(gdofsu);

    if (phytype & PhysicsType::_Thermal)
    {
        gdofsT = mesh.Nnd;
        KTT = Eigen::MatrixXd::Zero(gdofsT, gdofsT);
        MTT = Eigen::MatrixXd::Zero(gdofsT, gdofsT);
        FT = Eigen::VectorXd::Zero(gdofsT);
        dT = Eigen::VectorXd::Zero(gdofsT);
        KuT = Eigen::MatrixXd::Zero(gdofsu, gdofsT);
    }

    if ((phytype & PhysicsType::_Thermal) && (phytype & PhysicsType::_Chemical))
    {
        gdofsc = mesh.Nnd;
        Kcc = Eigen::MatrixXd::Zero(gdofsc, gdofsc);
        Mcc = Eigen::MatrixXd::Zero(gdofsc, gdofsc);
        Fc = Eigen::VectorXd::Zero(gdofsc);
        dc = Eigen::VectorXd::Zero(gdofsc);
        Kuc = Eigen::MatrixXd::Zero(gdofsu, gdofsc);
        KcT = Eigen::MatrixXd::Zero(gdofsc, gdofsT);
        KTc = Eigen::MatrixXd::Zero(gdofsT, gdofsc);
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::solve()
{
    if (phytype == PhysicsType::_Mechanics)
        solveu();
    else if (phytype == (PhysicsType::_Mechanics | PhysicsType::_Thermal))
        solveuT();
    else if (phytype == (PhysicsType::_Mechanics | PhysicsType::_Thermal | PhysicsType::_Chemical))
        solveuTc();
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::solveu()
{
    // Assemble stiffness and force
    assembleu();

    // Apply Dirichlet Boundary conditions
    apply_dirBC();

    // Solve linear system of equations using Cholesky
    dud = Kuud.ldlt().solve(Fud);

    // Reform full solution vector
    reform_full_sol();
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::solveuT()
{
    // Assemble stiffness and force
    assembleuT();

    // Apply Dirichlet Boundary conditions
    apply_dirBC();

    // Solve System of Equations
    double delta_t = end_time / nsteps;
    Eigen::VectorXd vTd;
    dTd = Eigen::VectorXd::Zero(ucdofsT);
    for (unsigned int istep = 1; istep <= nsteps; ++istep)
    {
        vTd = MTTd.ldlt().solve(FTd - KTTd * dTd);
        dTd += delta_t * vTd;

        dud = Kuud.ldlt().solve(Fud + KuTd * dTd);

        if (istep == 1 || (istep % 100 == 0))
        {
            reform_full_sol();
            post_print(std::string("2p"), istep);
        }
    }

    // Reform full solution vector
    reform_full_sol();
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::solveuTc()
{
    // Assemble stiffness and force
    assembleuTc();

    // Apply Dirichlet Boundary conditions
    apply_dirBC();

    // Solve System of Equations
    double delta_t = end_time / nsteps;
    Eigen::VectorXd vTd;
    Eigen::VectorXd vcd;
    dTd = Eigen::VectorXd::Zero(ucdofsT);
    dcd = Eigen::VectorXd::Zero(ucdofsc);
    for (unsigned int istep = 1; istep <= nsteps; ++istep)
    {
        vcd = Mccd.ldlt().solve(Fcd - Kccd * dcd - KcTd * dTd);
        vTd = MTTd.ldlt().solve(FTd - KTcd * dcd - KTTd * dTd);
 
        dcd += delta_t * vcd;
        dTd += delta_t * vTd;

        dud = Kuud.ldlt().solve(Fud + KuTd * dTd + Kucd * dcd);

        if (istep == 1 || (istep % 100 == 0))
        {
            reform_full_sol();
            post_print(std::string("3p"), istep);
        }
    }

    // Reform full solution vector
    reform_full_sol();
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::post_print(std::string file)
{
    //raw_sol_print();

    // Write vtk file with mesh and solutions
    std::cout << "*** Writing VTK file - " << file.c_str() << ".vtu\n";
    fem_to_vtk_vector(file, mesh.elem_conn, mesh.nodal_coords, du);

    if (phytype & PhysicsType::_Thermal)
    {
        std::string fileT = file + std::string("T");
        std::cout << "*** Writing VTK file - " << fileT.c_str() << ".vtu\n";
        fem_to_vtk_scalar(fileT, mesh.elem_conn, mesh.nodal_coords, dT);
    }

    if (phytype & PhysicsType::_Chemical)
    {
        std::string filec = file + std::string("c");
        std::cout << "*** Writing VTK file - " << filec.c_str() << ".vtu\n";
        fem_to_vtk_scalar(filec, mesh.elem_conn, mesh.nodal_coords, dc);
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::post_print(std::string file, unsigned int id)
{
    // Write vtk file with mesh and solutions
    std::string fileu = file + std::to_string(id);
    std::cout << "*** Writing VTK file - " << fileu.c_str() << ".vtu\n";
    fem_to_vtk_vector(fileu, mesh.elem_conn, mesh.nodal_coords, du);

    if (phytype & PhysicsType::_Thermal)
    {
        std::string fileT = file + std::string("T") + std::to_string(id);
        std::cout << "*** Writing VTK file - " << fileT.c_str() << ".vtu\n";
        fem_to_vtk_scalar(fileT, mesh.elem_conn, mesh.nodal_coords, dT);
    }

    if (phytype & PhysicsType::_Chemical)
    {
        std::string filec = file + std::string("c") + std::to_string(id);
        std::cout << "*** Writing VTK file - " << filec.c_str() << ".vtu\n";
        fem_to_vtk_scalar(filec, mesh.elem_conn, mesh.nodal_coords, dc);
    }
}


template<unsigned int dim>
void CoupledFieldFEM<dim>::init_shp_functions()
{
    // Construct quadrature points
    // And construct and store shape functions and their derivatives at those points

    if (mesh.mesh_type == MeshType::triangular)
    {
        if (quadrature_type == QuadratureType::one_point)
        {
            Nq = 1;
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = 1.0 / 3;
            Qpts(1, 0) = 1.0 / 3;

            wts.resize(Nq);
            wts(0) = 1.0;
        }
        else if (quadrature_type == QuadratureType::two_point)
        {
            Nq = 3;
            Qpts.resize(2, Nq);
            Qpts(0, 0) = 1.0 / 6;
            Qpts(1, 0) = 1.0 / 6;
            Qpts(0, 1) = 2.0 / 3;
            Qpts(1, 1) = 1.0 / 6;
            Qpts(0, 2) = 1.0 / 6;
            Qpts(1, 2) = 2.0 / 3;

            wts.resize(Nq);
            wts(0) = 1.0 / 3;
            wts(1) = 1.0 / 3;
            wts(2) = 1.0 / 3;
        }

        if (mesh.shp_func_type == ShapeFcnType::lagrange_linear)
        {
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_tri(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = LinearLagrange::dNdxi_tri(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = LinearLagrange::dNdxi_tri(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
    }
    else if (mesh.mesh_type == MeshType::quadrilateral)
    {
        if (quadrature_type == QuadratureType::one_point)
        {
            Nq = static_cast<unsigned int>(pow(1, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(1, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(1, 0);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(1, 0) * GsQuad::weights(1, 0);

        }
        else if (quadrature_type == QuadratureType::two_point)
        {
            Nq = static_cast<unsigned int>(pow(2, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 1) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 1) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 3) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 3) = GsQuad::quadrature_pt(2, 1);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(2, 0) * GsQuad::weights(2, 0);
            wts(1) = GsQuad::weights(2, 1) * GsQuad::weights(2, 0);
            wts(2) = GsQuad::weights(2, 1) * GsQuad::weights(2, 1);
            wts(3) = GsQuad::weights(2, 0) * GsQuad::weights(2, 1);
        }
        else if (quadrature_type == QuadratureType::three_point)
        {
            Nq = static_cast<unsigned int>(pow(3, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(3, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(3, 0);
            Qpts(0, 1) = GsQuad::quadrature_pt(3, 1);
            Qpts(1, 1) = GsQuad::quadrature_pt(3, 0);
            Qpts(0, 2) = GsQuad::quadrature_pt(3, 2);
            Qpts(1, 2) = GsQuad::quadrature_pt(3, 0);
            Qpts(0, 3) = GsQuad::quadrature_pt(3, 0);
            Qpts(1, 3) = GsQuad::quadrature_pt(3, 1);
            Qpts(0, 4) = GsQuad::quadrature_pt(3, 1);
            Qpts(1, 4) = GsQuad::quadrature_pt(3, 1);
            Qpts(0, 5) = GsQuad::quadrature_pt(3, 2);
            Qpts(1, 5) = GsQuad::quadrature_pt(3, 1);
            Qpts(0, 6) = GsQuad::quadrature_pt(3, 0);
            Qpts(1, 6) = GsQuad::quadrature_pt(3, 2);
            Qpts(0, 7) = GsQuad::quadrature_pt(3, 1);
            Qpts(1, 7) = GsQuad::quadrature_pt(3, 2);
            Qpts(0, 8) = GsQuad::quadrature_pt(3, 2);
            Qpts(1, 8) = GsQuad::quadrature_pt(3, 2);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(3, 0) * GsQuad::weights(3, 0);
            wts(1) = GsQuad::weights(3, 1) * GsQuad::weights(3, 0);
            wts(2) = GsQuad::weights(3, 2) * GsQuad::weights(3, 0);
            wts(3) = GsQuad::weights(3, 0) * GsQuad::weights(3, 1);
            wts(4) = GsQuad::weights(3, 1) * GsQuad::weights(3, 1);
            wts(5) = GsQuad::weights(3, 2) * GsQuad::weights(3, 1);
            wts(6) = GsQuad::weights(3, 0) * GsQuad::weights(3, 2);
            wts(7) = GsQuad::weights(3, 1) * GsQuad::weights(3, 2);
            wts(8) = GsQuad::weights(3, 2) * GsQuad::weights(3, 2);
        }

        if (mesh.shp_func_type == ShapeFcnType::lagrange_linear)
        {
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = LinearLagrange::dNdxi_quad(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = LinearLagrange::dNdxi_quad(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
        else if (mesh.shp_func_type == ShapeFcnType::lagrange_quadratic)
        {
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = QuadraticLagrange::N_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = QuadraticLagrange::dNdxi_quad(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = QuadraticLagrange::dNdxi_quad(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::assembleu()
{
    // Initialize local arrays and variables
    unsigned int ldofs = mesh.Nen * dim;
    Eigen::MatrixXd Kl;                      // Local Stiffness
    Eigen::VectorXd Fl;                      // Local Force
    Eigen::MatrixXd xe(dim, mesh.Nen);       // Local Coordinates of nodes
    Eigen::VectorXui gidx(ldofs);            // Global index of local dofs
    Eigen::VectorXd de(mesh.Nen * dim);      // Local nodal displacements
    Eigen::MatrixXd jac(dim, dim);           // Quadrature Isoparamteric Jacobian
    Eigen::MatrixXd jac_inv(dim, dim);       // Quadrature Isoparamteric Jacobian Inverse
    Eigen::MatrixXd dNdX(mesh.Nen, dim);     // Quadrature dN / dX matrices
    double det_jac;                          // Determinant of Isoparametric Jacobian

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        // Initialize local stiffness and force
        Kl = Eigen::MatrixXd::Zero(ldofs, ldofs);
        Fl = Eigen::VectorXd::Zero(ldofs);

        // Compute Elemental quantities (elem coords and solutions)
        compute_elem_quantities(conn, gidx, xe, de);

        // Quadrature Loop
        for (unsigned int iq = 0; iq < Nq; ++iq)
        {
            // Compute Quadrature Derivatives
            compute_quad_derivatives(iq, xe, jac, jac_inv, dNdX, det_jac);

            // Compute and add stiffness to local stiffness matrix
            compute_quad_stiff(iq, dNdX, det_jac, Kl);

            // Not impl compute_quad_force()
        }

        // Global Assembly Kl -> K, Fl -> F
        for (unsigned int jdx = 0; jdx < ldofs; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofs; ++idx)
                Kuu(gidx(idx), gidx(jdx)) += Kl(idx, jdx);

            Fu(gidx(jdx)) += Fl(jdx);
        }
    }
}


template<unsigned int dim>
void CoupledFieldFEM<dim>::assembleuT()
{
    // Initialize local arrays and variables
    unsigned int ldofs = mesh.Nen * dim;
    Eigen::MatrixXd Kl;                      // Local Stiffness
    Eigen::VectorXd Fl;                      // Local Force
    Eigen::MatrixXd xe(dim, mesh.Nen);       // Local Coordinates of nodes
    Eigen::VectorXui gidx(ldofs);            // Global index of local dofs
    Eigen::VectorXd de(mesh.Nen * dim);      // Local nodal displacements
    Eigen::MatrixXd jac(dim, dim);           // Quadrature Isoparamteric Jacobian
    Eigen::MatrixXd jac_inv(dim, dim);       // Quadrature Isoparamteric Jacobian Inverse
    Eigen::MatrixXd dNdX(mesh.Nen, dim);     // Quadrature dN / dX matrices
    double det_jac;                          // Determinant of Isoparametric Jacobian

    unsigned int ldofsT = mesh.Nen;
    Eigen::MatrixXd KuTl;
    Eigen::MatrixXd KTTl;
    Eigen::MatrixXd MTTl;
    Eigen::VectorXd FTl;
    Eigen::VectorXui gidxT(ldofsT);
    Eigen::VectorXd Te(ldofsT);

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        // Initialize local stiffness and force
        Kl = Eigen::MatrixXd::Zero(ldofs, ldofs);
        Fl = Eigen::VectorXd::Zero(ldofs);
        KTTl = Eigen::MatrixXd::Zero(ldofsT, ldofsT);
        KuTl = Eigen::MatrixXd::Zero(ldofs, ldofsT);
        MTTl = Eigen::MatrixXd::Zero(ldofsT, ldofsT);
        FTl = Eigen::VectorXd::Zero(ldofsT);

        // Compute Elemental quantities (elem coords and solutions)
        compute_elem_quantities(conn, gidx, xe, de);
        compute_elem_quantitiesT(conn, gidxT, Te);

        // Quadrature Loop
        for (unsigned int iq = 0; iq < Nq; ++iq)
        {
            // Compute Quadrature Derivatives
            compute_quad_derivatives(iq, xe, jac, jac_inv, dNdX, det_jac);

            // Compute and add stiffness to local stiffness matrix
            compute_quad_stiff(iq, dNdX, det_jac, Kl);
            compute_quad_stiffTT(iq, dNdX, det_jac, KTTl);
            compute_quad_stiffuT(iq, dNdX, det_jac, KuTl);

            // Compute and add to local mass matrix
            compute_quad_massTT(iq, det_jac, MTTl);
            

            // Not impl compute_quad_force()
        }

        // Global Assembly Kl -> K, Fl -> F
        for (unsigned int jdx = 0; jdx < ldofs; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofs; ++idx)
                Kuu(gidx(idx), gidx(jdx)) += Kl(idx, jdx);

            Fu(gidx(jdx)) += Fl(jdx);
        }

        for (unsigned int jdx = 0; jdx < ldofsT; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofsT; ++idx)
            {
                KTT(gidxT(idx), gidxT(jdx)) += KTTl(idx, jdx);
                MTT(gidxT(idx), gidxT(jdx)) += MTTl(idx, jdx);
            }

            FT(gidxT(jdx)) += FTl(jdx);

            for (unsigned int idx = 0; idx < ldofs; ++idx)
                KuT(gidx(idx), gidxT(jdx)) += KuTl(idx, jdx);
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::assembleuTc()
{
    // Initialize local arrays and variables
    unsigned int ldofs = mesh.Nen * dim;
    Eigen::MatrixXd Kl;                      // Local Stiffness
    Eigen::VectorXd Fl;                      // Local Force
    Eigen::MatrixXd xe(dim, mesh.Nen);       // Local Coordinates of nodes
    Eigen::VectorXui gidx(ldofs);            // Global index of local dofs
    Eigen::VectorXd de(mesh.Nen * dim);      // Local nodal displacements
    Eigen::MatrixXd jac(dim, dim);           // Quadrature Isoparamteric Jacobian
    Eigen::MatrixXd jac_inv(dim, dim);       // Quadrature Isoparamteric Jacobian Inverse
    Eigen::MatrixXd dNdX(mesh.Nen, dim);     // Quadrature dN / dX matrices
    double det_jac;                          // Determinant of Isoparametric Jacobian

    unsigned int ldofsT = mesh.Nen;
    Eigen::MatrixXd KuTl;
    Eigen::MatrixXd KTTl;
    Eigen::MatrixXd MTTl;
    Eigen::VectorXd FTl;
    Eigen::VectorXui gidxT(ldofsT);
    Eigen::VectorXd Te(ldofsT);

    unsigned int ldofsc = mesh.Nen;
    Eigen::MatrixXd Kucl;
    Eigen::MatrixXd KcTl;
    Eigen::MatrixXd KTcl;
    Eigen::MatrixXd Kccl;
    Eigen::MatrixXd Mccl;
    Eigen::VectorXd Fcl;
    Eigen::VectorXui gidxc(ldofsc);
    Eigen::VectorXd  ce(ldofsc);

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        // Initialize local stiffness and force
        Kl = Eigen::MatrixXd::Zero(ldofs, ldofs);
        Fl = Eigen::VectorXd::Zero(ldofs);
        KTTl = Eigen::MatrixXd::Zero(ldofsT, ldofsT);
        KuTl = Eigen::MatrixXd::Zero(ldofs, ldofsT);
        MTTl = Eigen::MatrixXd::Zero(ldofsT, ldofsT);
        FTl = Eigen::VectorXd::Zero(ldofsT);
        Kccl = Eigen::MatrixXd::Zero(ldofsc, ldofsc);
        Kucl = Eigen::MatrixXd::Zero(ldofs, ldofsc);
        KcTl = Eigen::MatrixXd::Zero(ldofsc, ldofsT);
        KTcl = Eigen::MatrixXd::Zero(ldofsT, ldofsc);
        Mccl = Eigen::MatrixXd::Zero(ldofsc, ldofsc);
        Fcl = Eigen::VectorXd::Zero(ldofsc);

        // Compute Elemental quantities (elem coords and solutions)
        compute_elem_quantities(conn, gidx, xe, de);
        compute_elem_quantitiesT(conn, gidxT, Te);
        compute_elem_quantitiesc(conn, gidxc, ce);

        // Quadrature Loop
        for (unsigned int iq = 0; iq < Nq; ++iq)
        {
            // Compute Quadrature Derivatives
            compute_quad_derivatives(iq, xe, jac, jac_inv, dNdX, det_jac);

            // Compute and add stiffness to local stiffness matrix
            compute_quad_stiff(iq, dNdX, det_jac, Kl);
            compute_quad_stiffTT(iq, dNdX, det_jac, KTTl);
            compute_quad_stiffuT(iq, dNdX, det_jac, KuTl);
            compute_quad_stiffcc(iq, dNdX, det_jac, Kccl);
            compute_quad_stiffuc(iq, dNdX, det_jac, Kucl);
            compute_quad_stiffcT(iq, dNdX, det_jac, KcTl, KTcl);

            // Compute and add to local mass matrix
            compute_quad_massTT(iq, det_jac, MTTl);
            compute_quad_massTT(iq, det_jac, Mccl);


            // Not impl compute_quad_force()
        }

        // Global Assembly Kl -> K, Fl -> F
        for (unsigned int jdx = 0; jdx < ldofs; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofs; ++idx)
                Kuu(gidx(idx), gidx(jdx)) += Kl(idx, jdx);

            Fu(gidx(jdx)) += Fl(jdx);
        }

        for (unsigned int jdx = 0; jdx < ldofsT; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofsT; ++idx)
            {
                KTT(gidxT(idx), gidxT(jdx)) += KTTl(idx, jdx);
                MTT(gidxT(idx), gidxT(jdx)) += MTTl(idx, jdx);
            }

            FT(gidxT(jdx)) += FTl(jdx);

            for (unsigned int idx = 0; idx < ldofs; ++idx)
                KuT(gidx(idx), gidxT(jdx)) += KuTl(idx, jdx);

            for (unsigned int idx = 0; idx < ldofsc; ++idx)
                KcT(gidxc(idx), gidxT(jdx)) += KcTl(idx, jdx);
        }

        for (unsigned int jdx = 0; jdx < ldofsc; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofsc; ++idx)
            {
                Kcc(gidxc(idx), gidxc(jdx)) += Kccl(idx, jdx);
                Mcc(gidxc(idx), gidxc(jdx)) += Mccl(idx, jdx);
            }

            Fc(gidxc(jdx)) += Fcl(jdx);

            for (unsigned int idx = 0; idx < ldofs; ++idx)
                Kuc(gidx(idx), gidxc(jdx)) += Kucl(idx, jdx);

            for (unsigned int idx = 0; idx < ldofsT; ++idx)
                KTc(gidxT(idx), gidxc(jdx)) += KTcl(idx, jdx);
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx,
    Eigen::MatrixXd& xe, Eigen::VectorXd& de)
{
    // Computes elemental nodal coordinates and displacements
    for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            auto loc_idx = dim * ind + idim;
            gidx(loc_idx) = dim * conn[ind] + idim;
            xe(idim, ind) = mesh.nodal_coords[conn[ind]][idim];
            de(loc_idx) = du(gidx(loc_idx));
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_elem_quantitiesT(connectivity& conn, Eigen::VectorXui& gidxT, Eigen::VectorXd& Te)
{
    // Computes global index and local temperatures
    for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
    {
        gidxT(ind) = conn[ind];
        Te(ind) = dT(gidxT(ind));
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_elem_quantitiesc(connectivity& conn, Eigen::VectorXui& gidxc, Eigen::VectorXd& ce)
{
    // Computes global index and local concentration
    for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
    {
        gidxc(ind) = conn[ind];
        ce(ind) = dc(gidxc(ind));
    }
}


template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::MatrixXd& jac,
    Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac)
{
    // Compute isoparametric Jacobian and Inverse
    jac.setZero();
    for (unsigned int jdim = 0; jdim < dim; ++jdim)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
                jac(idim, jdim) += xe(idim, ind) * dNdxi[jdim](ind, iq);
        }
    }

    det_jac = jac.determinant();
    jac_inv = jac.inverse();

    // Compute dNdX from dNdxi
    Eigen::MatrixXd dNdxi_q(mesh.Nen, dim);
    for (unsigned int idim = 0; idim < dim; ++idim)
    {
        for (unsigned int ishp = 0; ishp < mesh.Nen; ++ishp)
            dNdxi_q(ishp, idim) = dNdxi[idim](ishp, iq);
    }

    dNdX = dNdxi_q * jac_inv;
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kl)
{
    unsigned int loc_idx;
    unsigned int loc_kdx;
    double temp;

    // Compute Local Stiffness (-> integral (w_i,j  C_ijkl  u_k,l) dV)
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int kdim = 0; kdim < dim; ++kdim)
        {
            loc_kdx = dim * bshp + kdim;
            for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
            {
                for (unsigned int idim = 0; idim < dim; ++idim)
                {
                    loc_idx = dim * ashp + idim;
                    temp = 0.0;
                    for (unsigned int jdim = 0; jdim < dim; ++jdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                            temp += dNdX(ashp, jdim) * mat.C[idim][jdim][kdim][ldim] * dNdX(bshp, ldim);
                    }

                    Kl(loc_idx, loc_kdx) += temp * wts(iq) * det_jac;
                }
            }
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_stiffTT(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& KTTl)
{
    double temp, temp2;

    // Compute Local Temperature Stiffness (-> integral (w^T_,i  D^T (C_jkmn a_jk a_mn + 1)  T_,i) dV)
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            temp = 0.0;
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                temp2 = 0.0;

                // C_jkmn alpha_jk alpha_mn
                for (unsigned int jdim = 0; jdim < dim; ++jdim)
                {
                    for (unsigned int kdim = 0; kdim < dim; ++kdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                        {
                            for (unsigned int mdim = 0; mdim < dim; ++mdim)
                                temp2 += mat.C[jdim][kdim][ldim][mdim] * mat.alpha_T[jdim][kdim] * mat.alpha_T[ldim][mdim];
                        }
                    }
                }

                // dNdx(A, i) * (C_
                temp += dNdX(ashp, idim) * mat.D_T * (temp2 + 1.0) * dNdX(bshp, idim);
            }

            KTTl(ashp, bshp) += temp * wts(iq) * det_jac;
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_stiffuT(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& KuTl)
{
    unsigned int loc_idx;
    double temp, temp2;

    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                loc_idx = dim * ashp + idim;
                temp = 0.0;
                for (unsigned int jdim = 0; jdim < dim; ++jdim)
                {
                    temp2 = 0.0;
                    for (unsigned int kdim = 0; kdim < dim; ++kdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                            temp2 += mat.C[idim][jdim][kdim][ldim] * mat.alpha_T[kdim][ldim];
                    }

                    temp += dNdX(ashp, jdim) * temp2 * N(bshp, iq);
                }

                KuTl(loc_idx, bshp) += temp * wts(iq) * det_jac;
            }
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_massTT(unsigned iq, double det_jac, Eigen::MatrixXd& MTTl)
{
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
            MTTl(ashp, bshp) += N(ashp, iq) * N(bshp, iq) * wts(iq) * det_jac;
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_stiffcc(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kccl)
{
    double temp, temp2;

    // Compute Local Temperature Stiffness (-> integral (w^T_,i  D^c (C_jkmn a_jk a_mn + 1)  T_,i) dV)
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            temp = 0.0;
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                temp2 = 0.0;

                // C_jkmn alpha_jk alpha_mn
                for (unsigned int jdim = 0; jdim < dim; ++jdim)
                {
                    for (unsigned int kdim = 0; kdim < dim; ++kdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                        {
                            for (unsigned int mdim = 0; mdim < dim; ++mdim)
                                temp2 += mat.C[jdim][kdim][ldim][mdim] * mat.alpha_c[jdim][kdim] * mat.alpha_c[ldim][mdim];
                        }
                    }
                }

                // dNdx(A, i) * (C_
                temp += dNdX(ashp, idim) * mat.D_c * (temp2 + 1.0) * dNdX(bshp, idim);
            }

            Kccl(ashp, bshp) += temp * wts(iq) * det_jac;
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_stiffuc(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kucl)
{
    unsigned int loc_idx;
    double temp, temp2;

    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                loc_idx = dim * ashp + idim;
                temp = 0.0;
                for (unsigned int jdim = 0; jdim < dim; ++jdim)
                {
                    temp2 = 0.0;
                    for (unsigned int kdim = 0; kdim < dim; ++kdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                            temp2 += mat.C[idim][jdim][kdim][ldim] * mat.alpha_c[kdim][ldim];
                    }

                    temp += dNdX(ashp, jdim) * temp2 * N(bshp, iq);
                }

                Kucl(loc_idx, bshp) += temp * wts(iq) * det_jac;
            }
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::compute_quad_stiffcT(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& KcTl,
                                                Eigen::MatrixXd& KTcl)
{
    double temp, temp2;
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            temp = 0.0;
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                temp2 = 0.0;
                for (unsigned int jdim = 0; jdim < dim; ++jdim)
                {
                    for (unsigned int kdim = 0; kdim < dim; ++kdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                        {
                            for (unsigned int mdim = 0; mdim < dim; ++mdim)
                                temp2 += mat.C[jdim][kdim][ldim][mdim] * mat.alpha_c[jdim][kdim] * mat.alpha_T[ldim][mdim];
                        }
                    }
                }

                temp += dNdX(ashp, idim) * temp2 * dNdX(bshp, idim);
            }

            KcTl(ashp, bshp) += mat.D_c * temp * wts(iq) * det_jac;
            KTcl(ashp, bshp) += mat.D_T * temp * wts(iq) * det_jac;
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::create_dirichlet_map()
{
    // Creates a local map of unconstrained dofs
    // dirBCs should be sorted for this to work currently
    unsigned int dirBCsize = static_cast<unsigned int>(dirBCs.size());
    ucdofsu = gdofsu - dirBCsize;       // unconstrained dofs
    dirmapu.reserve(ucdofsu);
    unsigned int offset = 0;
    auto dir_itr = dirBCs.begin();
    for (unsigned int idx = 0; idx < mesh.Nnd; ++idx)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            if (dir_itr != dirBCs.end())
            {
                if (idx == dir_itr->node && idim == dir_itr->dof)
                {
                    ++dir_itr;
                    continue;
                }
            }

            dirmapu.push_back(dim * idx + idim);
        }
    }

    if (phytype & PhysicsType::_Thermal)
    {
        dirBCsize = static_cast<unsigned int>(dirBCsT.size());
        ucdofsT = gdofsT - dirBCsize;       // unconstrained dofs
        dirmapT.reserve(ucdofsT);
        auto dir_itr = dirBCsT.begin();
        for (unsigned int idx = 0; idx < mesh.Nnd; ++idx)
        {
            if (dir_itr != dirBCsT.end())
            {
                if (idx == dir_itr->node)
                {
                    ++dir_itr;
                    continue;
                }
            }

            dirmapT.push_back(idx);
        }
    }

    if (phytype & PhysicsType::_Chemical)
    {
        dirBCsize = static_cast<unsigned int>(dirBCsc.size());
        ucdofsc = gdofsc - dirBCsize;       // unconstrained dofs
        dirmapc.reserve(ucdofsc);
        auto dir_itr = dirBCsc.begin();
        for (unsigned int idx = 0; idx < mesh.Nnd; ++idx)
        {
            if (dir_itr != dirBCsc.end())
            {
                if (idx == dir_itr->node)
                {
                    ++dir_itr;
                    continue;
                }
            }

            dirmapc.push_back(idx);
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::apply_dirBC()
{
    // Applies Dirchlet BC by forming Kd and Fd by removing constrained dofs from K and F

    create_dirichlet_map();

    Kuud = Eigen::MatrixXd::Zero(ucdofsu, ucdofsu);
    Fud = Eigen::VectorXd::Zero(ucdofsu);
    dud.resize(ucdofsu);
    for (unsigned int jdx = 0; jdx < ucdofsu; ++jdx)
    {
        for (unsigned int idx = 0; idx < ucdofsu; ++idx)
            Kuud(idx, jdx) = Kuu(dirmapu[idx], dirmapu[jdx]);

        Fud(jdx) = Fu(dirmapu[jdx]);
    }

    for (auto& dirBC : dirBCs)
    {
        unsigned gidx = dirBC.node * dim + dirBC.dof;
        du(gidx) = dirBC.val;
        for (unsigned int idx = 0; idx < ucdofsu; ++idx)
            Fud(idx) -= Kuu(dirmapu[idx], gidx) * dirBC.val;
    }

    if (phytype & PhysicsType::_Thermal)
    {
        KTTd = Eigen::MatrixXd::Zero(ucdofsT, ucdofsT);
        KuTd = Eigen::MatrixXd::Zero(ucdofsu, ucdofsT);
        MTTd = Eigen::MatrixXd::Zero(ucdofsT, ucdofsT);
        FTd = Eigen::VectorXd::Zero(ucdofsT);
        dTd.resize(ucdofsT);

        for (unsigned int jdx = 0; jdx < ucdofsT; ++jdx)
        {
            for (unsigned int idx = 0; idx < ucdofsT; ++idx)
            {
                KTTd(idx, jdx) = KTT(dirmapT[idx], dirmapT[jdx]);
                MTTd(idx, jdx) = MTT(dirmapT[idx], dirmapT[jdx]);
            }

            FTd(jdx) = FT(dirmapT[jdx]);
        }

        for (unsigned int jdx = 0; jdx < ucdofsT; ++jdx)
        {
            for (unsigned int idx = 0; idx < ucdofsu; ++idx)
                KuTd(idx, jdx) = KuT(dirmapu[idx], dirmapT[jdx]);
        }

        for (auto& dirBC : dirBCsT)
        {
            unsigned gidx = dirBC.node;
            dT(gidx) = dirBC.val;
            for (unsigned int idx = 0; idx < ucdofsT; ++idx)
                FTd(idx) -= KTT(dirmapT[idx], gidx) * dirBC.val;

            for (unsigned int idx = 0; idx < ucdofsu; ++idx)
                Fud(idx) += KuT(dirmapu[idx], gidx) * dirBC.val;
        }
    }

    if ((phytype & PhysicsType::_Chemical) && (phytype & PhysicsType::_Thermal))
    {
        Kccd = Eigen::MatrixXd::Zero(ucdofsc, ucdofsc);
        Kucd = Eigen::MatrixXd::Zero(ucdofsu, ucdofsc);
        Mccd = Eigen::MatrixXd::Zero(ucdofsc, ucdofsc);
        KcTd = Eigen::MatrixXd::Zero(ucdofsc, ucdofsT);
        KTcd = Eigen::MatrixXd::Zero(ucdofsT, ucdofsc);
        Fcd = Eigen::VectorXd::Zero(ucdofsc);
        dcd.resize(ucdofsc);

        for (unsigned int jdx = 0; jdx < ucdofsc; ++jdx)
        {
            for (unsigned int idx = 0; idx < ucdofsc; ++idx)
            {
                Kccd(idx, jdx) = Kcc(dirmapc[idx], dirmapc[jdx]);
                Mccd(idx, jdx) = Mcc(dirmapc[idx], dirmapc[jdx]);
            }

            for (unsigned int idx = 0; idx < ucdofsu; ++idx)
                Kucd(idx, jdx) = Kuc(dirmapu[idx], dirmapc[jdx]);

            for (unsigned int idx = 0; idx < ucdofsT; ++idx)
                KTcd(idx, jdx) = KTc(dirmapT[idx], dirmapc[jdx]);

            Fcd(jdx) = Fc(dirmapc[jdx]);
        }

        for (unsigned int jdx = 0; jdx < ucdofsT; ++jdx)
        {
            for (unsigned int idx = 0; idx < ucdofsc; ++idx)
                KcTd(idx, jdx) = KcT(dirmapc[idx], dirmapT[jdx]);
        }

        for (auto& dirBC : dirBCsc)
        {
            unsigned gidx = dirBC.node;
            dc(gidx) = dirBC.val;
            for (unsigned int idx = 0; idx < ucdofsc; ++idx)
                Fcd(idx) -= Kcc(dirmapc[idx], gidx) * dirBC.val;

            for (unsigned int idx = 0; idx < ucdofsu; ++idx)
                Fud(idx) += Kuc(dirmapu[idx], gidx) * dirBC.val;

            for (unsigned int idx = 0; idx < ucdofsT; ++idx)
                FTd(idx) -= KTc(dirmapT[idx], gidx) * dirBC.val;
        }

        for (auto& dirBC : dirBCsT)
        {
            unsigned gidx = dirBC.node;
            for (unsigned int idx = 0; idx < ucdofsc; ++idx)
                Fcd(idx) -= KcT(dirmapc[idx], gidx) * dirBC.val;
        }
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::reform_full_sol()
{
    // Constructs full solution vector from solutions at unconstrained dofs
    for (unsigned int idx = 0; idx < ucdofsu; ++idx)
        du(dirmapu[idx]) = dud(idx);

    if (phytype & PhysicsType::_Thermal)
    {
        for (unsigned int idx = 0; idx < ucdofsT; ++idx)
            dT(dirmapT[idx]) = dTd(idx);
    }

    if (phytype & PhysicsType::_Chemical)
    {
        for (unsigned int idx = 0; idx < ucdofsc; ++idx)
            dc(dirmapc[idx]) = dcd(idx);
    }
}

template<unsigned int dim>
void CoupledFieldFEM<dim>::raw_sol_print()
{
    // Prints nodal coordinates, element connectivity and solution field
    // Not used currently since fem_to_vtk is implemented natively in C++

    std::ofstream nc("nodal_coords.txt");
    if (nc.is_open())
    {
        for (unsigned int idx = 0; idx < mesh.Nnd - 1; ++idx)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
                nc << std::setprecision(16) << mesh.nodal_coords[idx][idim] << "\t";

            nc << "\n";
        }

        for (unsigned int idim = 0; idim < dim; ++idim)
            nc << std::setprecision(16) << mesh.nodal_coords[mesh.Nnd - 1][idim] << "\t";
        nc.close();
    }

    std::ofstream ec("elem_conn.txt");
    if (ec.is_open())
    {
        for (unsigned int idx = 0; idx < mesh.Nel - 1; ++idx)
        {
            for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
                ec << std::setprecision(16) << mesh.elem_conn[idx][ind] << "\t";

            ec << "\n";
        }

        for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
            ec << std::setprecision(16) << mesh.elem_conn[mesh.Nel - 1][ind] << "\t";
        ec.close();
    }

    std::ofstream sol("sol.txt");
    if (sol.is_open())
    {
        for (unsigned int ind = 0; ind < mesh.Nnd - 1; ++ind)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                unsigned int gidx = dim * ind + idim;
                sol << std::setprecision(16) << du(gidx) << "\t";
            }

            sol << "\n";
        }

        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            unsigned int gidx = dim * (mesh.Nnd - 1) + idim;
            sol << std::setprecision(16) << du(gidx) << "\t";
        }

        sol.close();
    }
}

////////////////////////////////////////////////////////////////////////////////
////                            Instantiation                               ////
////////////////////////////////////////////////////////////////////////////////

template class CoupledFieldFEM<2>;