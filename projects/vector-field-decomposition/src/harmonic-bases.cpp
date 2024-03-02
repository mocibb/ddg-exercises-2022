// Implement member functions for HarmonicBases class.
#include "harmonic-bases.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HarmonicBases::HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a closed, but not exact, primal 1-form ω.
 *
 * Input: A std::vector of Halfedges representing a homology generator of the mesh.
 * Returns: A vector representing a closed primal 1-form.
 */
Vector<double> HarmonicBases::buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const {
    Vector<double> omega = Vector<double>::Zero(mesh->nEdges());

    // generator from twin --> face
    for (Halfedge he : generator) {
        bool is_oriented = (he == he.edge().halfedge());
        omega[he.edge().getIndex()] = is_oriented ? 1 : -1;
    }
    return omega;
}


/*
 * Compute the harmonic bases [γ1, γ2 ... γn] of the input mesh.
 *
 * Input: A std::vector of homology generators of the mesh (which are in turn represented as std::vectors of halfedges),
 * and a HodgeDecomposition object. Returns:
 */
std::vector<Vector<double>> HarmonicBases::compute(const std::vector<std::vector<Halfedge>>& generators,
                                                   const HodgeDecomposition& hodgeDecomposition) const {

    std::vector<Vector<double>> gammas;
 
    for (size_t i = 0; i < generators.size(); i++) {
        const auto omega = buildClosedPrimalOneForm(generators[i]);
        const auto dAlpha = hodgeDecomposition.computeExactComponent(omega);
        gammas.push_back(omega - dAlpha);
    }
    
    return gammas;
}
