// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/meshio.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    SparseMatrix<double> A = geometry->laplaceMatrix();
    return M+h*A;
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> A = buildFlowOperator(M, h);

    geometrycentral::PositiveDefiniteSolver<double> solver(A);

    Eigen::MatrixXd vpos(mesh->nVertices(), 3);

    for (Vertex v : mesh->vertices()) {
        Vector3 p = geometry->inputVertexPositions[v];
        vpos.row(v.getIndex()) = Eigen::Vector3d(p.x, p.y, p.z);
    }

    for (int i = 0; i < 3; i++) {
        Vector<double> rhs1(M * vpos.col(i));
        vpos.col(i) = solver.solve(rhs1);
    }

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    for (Vertex v : mesh->vertices()) {
        Eigen::Vector3d p = vpos.row(v.getIndex());
        geometry->inputVertexPositions[v] = {p[0], p[1], p[2]};
    }
}
