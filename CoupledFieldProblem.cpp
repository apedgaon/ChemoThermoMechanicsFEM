#include "CoupledFieldFEM.h"

void Part1()
{
    // Initialize fem object
    CoupledFieldFEM<2> fem;

    // Geometry
    fem.geom.length[0] = 2.0;      // cm
    fem.geom.length[1] = 1.0;      // cm

    // Material Parameters
    fem.mat.lambda = 5.769230769230769e-05; // N/cm^2
    fem.mat.mu = 3.846153846153846e-05;     // N/cm^2
    fem.mat.init_elasticity_tensor();

    // Mesh Settings
    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 10;
    fem.mesh.divs[1] = 5;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)// u_x = 0 at x = 0.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
        else if (fabs(x - fem.geom.length[0]) < 1.0e-8)// u_x = 0.1 at x = 2.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.1;
            fem.dirBCs.push_back(dbc);
        }

        if (fabs(y - 0.0) < 1.0e-8)// u_y = 0 at x = 0.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 1;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
    }

    // Quadrature Settings
    fem.quadrature_type = QuadratureType::two_point;

    // Physics and solution settings
    fem.phytype = PhysicsType::_Mechanics;
    fem.nsteps = 1;

    // Initialize, Solve and Post processing
    fem.initialize();
    fem.solve();
    fem.post_print("part1");
}

void Part2()
{
    // Initialize fem object
    CoupledFieldFEM<2> fem;

    // Geometry
    fem.geom.length[0] = 2.0;      // cm
    fem.geom.length[1] = 1.0;      // cm

    // Material Parameters
    fem.mat.lambda = 5.769230769230769e-05; // N/cm^2
    fem.mat.mu = 3.846153846153846e-05;     // N/cm^2
    fem.mat.init_elasticity_tensor();
    fem.mat.D_T = 1e4;  // cm^2/s (= 0.01 W/cm-K x 10^6 cm^3-K/J)
    fem.mat.init_isotropic_thermal_exp(1.0);    // K^-1

    // Mesh Settings
    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 10;
    fem.mesh.divs[1] = 5;
    fem.generate_mesh();

    // Store Dirichlet BC (u and T)
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)// u_x = 0 at x = 0.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
        else if (fabs(x - fem.geom.length[0]) < 1.0e-8)// u_x = 0.1 at x = 2.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.1;
            fem.dirBCs.push_back(dbc);

            dbc.val = 1; // degree C
            fem.dirBCsT.push_back(dbc);// T = 1.0 at x = 2.0
        }

        if (fabs(y - 0.0) < 1.0e-8)// u_y = 0 at y = 0.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 1;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
    }

    // Quadrature Setting
    fem.quadrature_type = QuadratureType::two_point;

    // Physics and Solution Settings
    fem.phytype = PhysicsType::_Mechanics | PhysicsType::_Thermal;
    fem.sol_scheme = SolutionScheme::staggered;
    fem.nsteps = 10000;
    fem.end_time = 0.0005; // sec

    // Initialize, Solve and Post process
    fem.initialize();
    fem.solve2();
    fem.post_print("part2");
}

void Part3()
{
    // Initialize fem object
    CoupledFieldFEM<2> fem;

    // Geometry
    fem.geom.length[0] = 2.0;      // cm
    fem.geom.length[1] = 1.0;      // cm

    // Material paramters
    fem.mat.lambda = 5.769230769230769e-05; // N/cm^2
    fem.mat.mu = 3.846153846153846e-05;     // N/cm^2
    fem.mat.init_elasticity_tensor();
    fem.mat.D_T = 1e4;  // cm^2/s (= 0.01 W/cm-K x 10^6 cm^3-K/J)
    fem.mat.init_isotropic_thermal_exp(1.0);    // K^-1
    fem.mat.D_c = 1e4; // cm^2/s
    fem.mat.init_isotropic_chemical_exp(1.0);    // K^-1

    // Mesh Settings
    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 10;
    fem.mesh.divs[1] = 5;
    fem.generate_mesh();

    // Store Dirichlet BC (u, T and c)
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8) // u_x = 0 at x = 0.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
        else if (fabs(x - fem.geom.length[0]) < 1.0e-8) // u_x = 0.1 at x = 2.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.1;
            fem.dirBCs.push_back(dbc);

            DirichletBC tbc;    // T = 1.0 at x = 2.0
            tbc.node = idx;
            tbc.dof = 0;
            tbc.val = 1.0; // degree C
            fem.dirBCsT.push_back(tbc);
        }

        if (fabs(y - 0.0) < 1.0e-8) // u_y = 0 at y = 0.0
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 1;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
        else if (fabs(y - fem.geom.length[1]) < 1.0e-8) // c = 1.0 at y = 1.0
        {
            DirichletBC cbc;
            cbc.node = idx;
            cbc.dof = 0;
            cbc.val = 1.0;
            fem.dirBCsc.push_back(cbc);
        }
    }

    // Quadrature Setting
    fem.quadrature_type = QuadratureType::two_point;

    // Physics and Solution Settings
    fem.phytype = PhysicsType::_Mechanics | PhysicsType::_Thermal | PhysicsType::_Chemical;
    fem.sol_scheme = SolutionScheme::staggered;
    fem.nsteps = 10000;
    fem.end_time = 0.0005;

    // Initialize, Solve and Post process
    fem.initialize();
    fem.solve3();
    fem.post_print("part3");
}

int main(int argc, char* argv[])
{
    Part1(); // mechanics
    Part2(); // thermal + mechanics
    Part3(); // chemical + thermal + mechanics
    return 0;
}