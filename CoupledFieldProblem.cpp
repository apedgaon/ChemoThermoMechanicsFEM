#include "CoupledFieldFEM.h"

void Part1()
{
    CoupledFieldFEM<2> fem;
    fem.geom.length[0] = 2.0;      // cm
    fem.geom.length[1] = 1.0;      // cm

    fem.mat.lambda = 5.769230769230769e-05; // N/cm^2
    fem.mat.mu = 3.846153846153846e-05;     // N/cm^2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

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
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
        else if (fabs(x - fem.geom.length[0]) < 1.0e-8)
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 0;
            dbc.val = 0.1;
            fem.dirBCs.push_back(dbc);
        }

        if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC dbc;
            dbc.node = idx;
            dbc.dof = 1;
            dbc.val = 0.0;
            fem.dirBCs.push_back(dbc);
        }
    }

    fem.quadrature_type = QuadratureType::two_point;
    fem.phytype = PhysicsType::_Mechanics;

    fem.initialize();
    fem.solve();
    fem.post_print("part1");
}

void Part2()
{
    CoupledFieldFEM<2> fem;
    fem.geom.length[0] = 2.0;      // cm
    fem.geom.length[1] = 1.0;      // cm

    fem.mat.lambda = 5.769230769230769e-05; // N/cm^2
    fem.mat.mu = 3.846153846153846e-05;     // N/cm^2
    fem.mat.init_elasticity_tensor();
    fem.mat.init_isotropic_thermal_exp(1.0);    // K^-1

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 10;
    fem.mesh.divs[1] = 5;
    fem.generate_mesh();

    fem.quadrature_type = QuadratureType::two_point;
    fem.phytype = PhysicsType::_Mechanics | PhysicsType::_Thermal;
    fem.sol_scheme = SolutionScheme::staggered;

    fem.initialize();
}

int main(int argc, char* argv[])
{
    //Part1();
    Part2();
    return 0;
}