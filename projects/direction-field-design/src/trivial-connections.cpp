// Implement member functions for TrivialConnections class.
#include "trivial-connections.h"
#include "tree-cotree.h"
#include "harmonic-bases.h"
#include "hodge-decomposition.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TrivialConnections::TrivialConnections(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // Build harmonic bases
    TreeCotree treeCotree(mesh, geometry);
    treeCotree.buildGenerators();
    HodgeDecomposition hodgeDecomposition(mesh, geometry);
    HarmonicBases harmonicBases(mesh, geometry);
    
    this->generators = treeCotree.generators;
    this->bases = harmonicBases.compute(generators, hodgeDecomposition);

    // Build period matrix.
    this->P = this->buildPeriodMatrix();

    // Store DEC operators
    this->A = hodgeDecomposition.A;
    this->hodge1 = hodgeDecomposition.hodge1;
    this->d0 = hodgeDecomposition.d0;
}

/*
 * Builds the period matrix Pij = ∑_{ek ∈ li} (ξj)k, where li is the ith homology generator, ek is a dual edge in li and
 * ξj is the jth harmonic 1-form basis.
 *
 * Input:
 * Returns: A sparse matrix represending the period matrix.
 */
SparseMatrix<double> TrivialConnections::buildPeriodMatrix() const {
    // genus x 2
    int g2 = this->bases.size();
    SparseMatrix<double> P(g2, g2);    

    for (size_t i = 0; i < this->generators.size(); i++) {
        for (size_t j = 0; j < this->bases.size(); j++) {
            double v = 0;
            for (Halfedge he: this->generators[i]) {
                v += (he == he.edge().halfedge() ? 1 : -1) * this->bases[j][he.edge().getIndex()];
            }
            
            P.coeffRef(i, j) = v;
        }
    }
    return P;
}

/*
 * Determine if a mesh satisfies Gauss-Bonnet.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: True if mesh satisfies Gauss-Bonnet, false otherwise.
 */
bool TrivialConnections::satsifyGaussBonnet(const Vector<double>& singularity) const {

    return (abs(singularity.sum() - geometry->eulerCharacteristic()) < 1e-8);
}

/*
 * Compute the dual 0-form potential β by solving the system d𝛿β = -K + 2π * singularity.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: The coexact component 𝛿β.
 */
Vector<double> TrivialConnections::computeCoExactComponent(const Vector<double>& singularity) const {
    Vector<double> u = Vector<double>::Zero(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        auto idx = v.getIndex();
        u[idx] = -geometry->angleDefect(v) + 2*PI*singularity[idx];
    }

    SparseMatrix<double> A = this->A;
    Vector<double> betaTilde = solvePositiveDefinite(A, u);
    return this->hodge1*this->d0*betaTilde;
}


/*
 * Given an initial angle αi in face i, this function computes the new angle αj in the neighboring face j as
 * αj = αi - θij + θji, where θij and θji are the angles between the shared edge e and an arbitrary but fixed reference
 * direction in faces i and j. Repeating this procedure for n consecutive dual edges in a generator gives a sequence of
 * angles α0, . . . , αn with a resulting total angle defect equal to αn - α0. This corresponds to transporting a vector
 * around a generator by unfolding, sliding and refolding it across neighboring faces without any extra in plane
 * rotation.
 *
 * Input: A halfedge lying on the shared edge between face i and j, and the initial angle αi.
 * Returns: The new angle αj.
 */
double TrivialConnections::transportNoRotation(Halfedge he, double alphaI) const {

    Vector3 u = geometry->halfedgeVector(he);

    // Compute two orthonormal tangent vectors for each face.
    Face fi = he.face();
    Face fj = he.twin().face();
    Vector3 e1 = geometry->halfedgeVector(fi.halfedge()).normalize();
    Vector3 e2 = cross(geometry->faceNormal(fi), e1);
    Vector3 f1 = geometry->halfedgeVector(fj.halfedge()).normalize();
    Vector3 f2 = cross(geometry->faceNormal(fj), f1);
    double thetaIJ = atan2(dot(u, e2), dot(u, e1));
    double thetaJI = atan2(dot(u, f2), dot(u, f1));

    return alphaI - thetaIJ + thetaJI;
}

/*
 * Compute the harmonic component γ = ∑_{i = 1, ..., 2g} zi ξi by solving the system Pz = v - ∑𝛿β.
 * v - ∑𝛿β should be normalized to lie between -π and π.
 *
 * Input: The coexact component 𝛿β.
 * Returns: The harmonic component γ.
 */
Vector<double> TrivialConnections::computeHarmonicComponent(const Vector<double>& deltaBeta) const {
    Vector<double> gamma = Vector<double>::Zero(mesh->nEdges());

    if (this->bases.size() > 0) {
        Vector<double> vtilde(this->generators.size());
        for (size_t i = 0; i < this->generators.size(); i++) {
            double sum = 0;
            for (Halfedge he: this->generators[i]) {
                sum += transportNoRotation(he, 0);
                sum -= (he == he.edge().halfedge() ? 1 : -1) * deltaBeta[he.edge().getIndex()];
            }
            vtilde[i] = sum - 2*PI*floor(sum / (2*PI));
        }

        SparseMatrix<double> P = this->P;
        
        Vector<double> z = solveSquare(P, vtilde);
        
        SparseMatrix<double> d1 = geometry->buildExteriorDerivative1Form();
        SparseMatrix<double> hodge1 = geometry->buildHodgeStar1Form();
        SparseMatrix<double> d0T = geometry->buildExteriorDerivative0Form().transpose();
        for (size_t i = 0; i < this->bases.size(); i++) {
            gamma += z[i] * bases[i];
            // check
            Vector<double> dGamma = d1 * bases[i];
            Vector<double> deltaGamma = d0T * hodge1 * bases[i];
            if (dGamma.norm() > 1e-5) {
                std::cout << "dGamma.norm() = " << dGamma.norm() << std::endl;
            }
            if (deltaGamma.norm() > 1e-5) {
                std::cout << "deltaGamma.norm() = " << deltaGamma.norm() << std::endl;
            }
        }
    }
    
    return gamma;
}

/*
 * Compute the dual 1-form connections φ = 𝛿β + γ.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: A vector representing the connections.
 */
Vector<double> TrivialConnections::computeConnections(const Vector<double>& singularity) const {

    if (!this->satsifyGaussBonnet(singularity)) {
        std::cerr << "Singularities do not add up to the Euler characteristic of the mesh" << std::endl;
        return Vector<double>::Zero(mesh->nEdges());
    }

    Vector<double> deltaBeta = computeCoExactComponent(singularity);
    Vector<double> gamma = computeHarmonicComponent(-deltaBeta);
    return deltaBeta + gamma;
}

