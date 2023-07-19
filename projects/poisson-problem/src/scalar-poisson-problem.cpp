// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: Build member variables A (Laplace matrix), M (mass matrix), total area
    this->A = geometry->laplaceMatrix();
    this->M = geometry->massMatrix();
    this->totalArea = geometry->totalArea();
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), 
 * where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    
    SparseMatrix<double> A = this->A;

    double rhoBar = 0;
    for (Eigen::Index i = 0; i < rho.rows(); i++) {
        if (rho[i]) {
            rhoBar += rho[i]*this->M.coeff(i, i);
        }
    }
    rhoBar /= this->totalArea;

    geometrycentral::PositiveDefiniteSolver<double> solver(A);
    Vector<double> rhs1(-this->M * (rho - rhoBar * Vector<double>::Ones(this->A.rows())));
    Vector<double> sol = solver.solve(rhs1);

    return sol;
}
