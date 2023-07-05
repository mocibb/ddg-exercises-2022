// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    Eigen::SparseMatrix<double> H0(mesh.nVertices(), mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tl;
    tl.reserve(mesh.nVertices());

    for (Vertex v : mesh.vertices()) {
        tl.push_back(T(v.getIndex(), v.getIndex(), barycentricDualArea(v)));
    }

    H0.setFromTriplets(tl.begin(), tl.end());
    return H0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    Eigen::SparseMatrix<double> H1(mesh.nEdges(), mesh.nEdges());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tl;
    tl.reserve(mesh.nEdges());

    for (Edge e : mesh.edges()) {
        // 外心和边长之比
        tl.push_back(T(e.getIndex(), e.getIndex(),
                       (halfedgeCotanWeight(e.halfedge()) + halfedgeCotanWeight(e.halfedge().twin()))));
    }

    H1.setFromTriplets(tl.begin(), tl.end());
    return H1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    Eigen::SparseMatrix<double> H2(mesh.nFaces(), mesh.nFaces());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tl;
    tl.reserve(mesh.nFaces());

    for (Face f : mesh.faces()) {
        tl.push_back(T(f.getIndex(), f.getIndex(), 1. / faceArea(f)));
    }

    H2.setFromTriplets(tl.begin(), tl.end());
    return H2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    Eigen::SparseMatrix<double> D0(mesh.nEdges(), mesh.nVertices());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tl;
    tl.reserve(mesh.nEdges() * 2);

    for (Edge e : mesh.edges()) {
        size_t idx = e.getIndex();
        tl.push_back(T(idx, e.firstVertex().getIndex(), -1.0));
        tl.push_back(T(idx, e.secondVertex().getIndex(), 1.0));
    }

    D0.setFromTriplets(tl.begin(), tl.end());
    return D0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    Eigen::SparseMatrix<double> D1(mesh.nFaces(), mesh.nEdges());
    typedef Eigen::Triplet<double> T;
    std::vector<T> tl;
    tl.reserve(mesh.nFaces() * 3);

    bool is_oriented;
    for (Face f : mesh.faces()) {
        size_t idx = f.getIndex();
        Halfedge he = f.halfedge();
        is_oriented = he.tipVertex().getIndex() == he.edge().firstVertex().getIndex();
        tl.push_back(T(idx, he.edge().getIndex(), is_oriented ? 1 : -1));
        he = he.next();
        is_oriented = he.tipVertex().getIndex() == he.edge().firstVertex().getIndex();
        tl.push_back(T(idx, he.edge().getIndex(), is_oriented ? 1 : -1));
        he = he.next();
        is_oriented = he.tipVertex().getIndex() == he.edge().firstVertex().getIndex();
        tl.push_back(T(idx, he.edge().getIndex(), is_oriented ? 1 : -1));
    }

    D1.setFromTriplets(tl.begin(), tl.end());
    return D1;
}

} // namespace surface
} // namespace geometrycentral