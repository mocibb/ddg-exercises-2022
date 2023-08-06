// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    double edgeLength = this->geometry->meanEdgeLength();
    double timeStep = edgeLength*edgeLength;
    int N = mesh->nVertices();
    SparseMatrix<double> M = geometry->massMatrix();
    // Note: core/geometry.cpp has meanEdgeLength() function
    this->A = geometry->laplaceMatrix();
    this->F = M + timeStep * A;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {
    FaceData<Vector3> vf(*mesh, {0, 0, 0});

    for (Face f : mesh->faces()) {
        Vector3 grad{0, 0, 0};
        Vector3 normal = geometry->faceNormal(f);

        for (Halfedge he : f.adjacentHalfedges()) {
            Vector3 ePerp = geometry->inputVertexPositions[he.next().tipVertex()] - geometry->inputVertexPositions[he.next().tailVertex()];
            ePerp = ePerp.rotateAround(normal, M_PI/2);
            grad += ePerp * u[he.vertex().getIndex()];
        }
        vf[f.getIndex()] = -grad.normalizeCutoff();
    }

    return vf;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {
    Vector<double> div = Vector<double>::Zero(mesh->nVertices());

    for (Face f : mesh->faces()) {
        Vector3 Xj = X[f.getIndex()];
        for (Halfedge he : f.adjacentHalfedges()) {
            Vector3 e = geometry->inputVertexPositions[he.tipVertex()] - geometry->inputVertexPositions[he.tailVertex()];
            double val = 0.5*geometry->cotan(he) * dot(e, Xj);
            div[he.tailVertex().getIndex()] += val;
            div[he.tipVertex().getIndex()] += -val;
        }
    }
    return div; // placeholder
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {
    Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
    Vector<double> u = llt.solve(delta);
    FaceData<Vector3> X = computeVectorField(u);
    Vector<double> deltaPhi = computeDivergence(X);

    SparseMatrix<double> A = this->A;
    geometrycentral::PositiveDefiniteSolver<double> solver(A);
    Vector<double> phi = solver.solve(-deltaPhi);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}
