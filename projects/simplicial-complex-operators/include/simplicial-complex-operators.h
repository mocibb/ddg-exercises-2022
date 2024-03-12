#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "mesh_subset.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class SimplicialComplexOperators {

  public:
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    // Using size_t instead of int for (unsigned) adjacency matrices is probably overkill, but doesn't hurt
    SparseMatrix<size_t> A0; // vertex-edge adjacency matrix
    SparseMatrix<size_t> A1; // edge-face adjacency matrix

    // Constructors
    SimplicialComplexOperators(){};
    SimplicialComplexOperators(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
        mesh = inputMesh;
        geometry = inputGeo;
        // this->assignElementIndices();
        A0 = this->buildVertexEdgeAdjacencyMatrix();
        A1 = this->buildFaceEdgeAdjacencyMatrix();
    }

    // Initialize
    void initialize(ManifoldSurfaceMesh* M, VertexPositionGeometry* G) {
        mesh = M;
        geometry = G;
        // this->assignElementIndices();
        A0 = this->buildVertexEdgeAdjacencyMatrix();
        A1 = this->buildFaceEdgeAdjacencyMatrix();
    }

    // Destructor
    ~SimplicialComplexOperators() {
        delete mesh;
        delete geometry;
    }

    inline std::set<size_t> copyFromVector(Vector<size_t> v) const {
        std::set<size_t> ss;
        for (size_t idx = 0; idx < (size_t)v.size(); idx++) {
            if (v[idx]) {
                ss.insert(idx);
            }
        }

        return ss;
    }

    // Operators
    void assignElementIndices();
    SparseMatrix<size_t> buildVertexEdgeAdjacencyMatrix() const;
    SparseMatrix<size_t> buildFaceEdgeAdjacencyMatrix() const;
    Vector<size_t> buildVertexVector(const MeshSubset& subset) const;
    Vector<size_t> buildEdgeVector(const MeshSubset& subset) const;
    Vector<size_t> buildFaceVector(const MeshSubset& subset) const;
    MeshSubset star(const MeshSubset& subset) const;
    MeshSubset closure(const MeshSubset& subset) const;
    MeshSubset link(const MeshSubset& subset) const;
    bool isComplex(const MeshSubset& subset) const;
    int isPureComplex(const MeshSubset& subset) const;
    MeshSubset boundary(const MeshSubset& subset) const;
};
