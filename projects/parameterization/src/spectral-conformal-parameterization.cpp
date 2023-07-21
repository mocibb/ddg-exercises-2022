// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {
    int N = mesh->nVertices();

    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    for (Face f : mesh->faces()) {
        for(Halfedge he : f.adjacentHalfedges()) {
            int from = he.tailVertex().getIndex();
            int to = he.tipVertex().getIndex();
            triplets.emplace_back(from, to, std::complex<double>(0, -0.25));
            triplets.emplace_back(to, from, std::complex<double>(0, 0.25));
        }
    }

    SparseMatrix<std::complex<double>> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());

    SparseMatrix<std::complex<double>> Ed = geometry->complexLaplaceMatrix();
    
    return Ed/2 - A;
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {
    SparseMatrix<std::complex<double>> Ec = buildConformalEnergy();
    Vector<std::complex<double>> x = solveInversePowerMethod(Ec);
    VertexData<Vector2> flattening = VertexData<Vector2>(*mesh);

    for (Vertex v : mesh->vertices()) {
      flattening[v] = Vector2::fromComplex(x[v.getIndex()]);
    }

    return flattening;
}
