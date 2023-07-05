// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <unordered_map>

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    Eigen::SparseMatrix<size_t> A0(mesh->nEdges(), mesh->nVertices());
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tl;
    tl.reserve(mesh->nEdges() * 2);

    size_t idx = 0;
    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
        tl.push_back(T(idx, geometry->vertexIndices[e.firstVertex()], 1));
        tl.push_back(T(idx, geometry->vertexIndices[e.secondVertex()], 1));
    }
    A0.setFromTriplets(tl.begin(), tl.end());

    return A0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    Eigen::SparseMatrix<size_t> A1(mesh->nFaces(), mesh->nEdges());
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tl;
    tl.reserve(mesh->nFaces() * 2);

    size_t idx = 0;
    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
        tl.push_back(T(idx, geometry->edgeIndices[f.halfedge().edge()], 1));
        tl.push_back(T(idx, geometry->edgeIndices[f.halfedge().next().edge()], 1));
        tl.push_back(T(idx, geometry->edgeIndices[f.halfedge().next().next().edge()], 1));
        // tl.push_back(T(idx, geometry->edgeIndices[e.secondVertex()], 1));
    }
    A1.setFromTriplets(tl.begin(), tl.end());

    return A1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    Vector<size_t> vertices = Vector<size_t>::Zero(mesh->nVertices());
    for (std::set<size_t>::iterator it = subset.vertices.begin(); it != subset.vertices.end(); ++it) {
        vertices[*it] = 1;
    }
    return vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> edges = Vector<size_t>::Zero(mesh->nEdges());
    for (std::set<size_t>::iterator it = subset.edges.begin(); it != subset.edges.end(); ++it) {
        edges[*it] = 1;
    }
    return edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> faces = Vector<size_t>::Zero(mesh->nFaces());
    for (std::set<size_t>::iterator it = subset.faces.begin(); it != subset.faces.end(); ++it) {
        faces[*it] = 1;
    }
    return faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;

    for (std::set<size_t>::iterator it = subset.vertices.begin(); it != subset.vertices.end(); ++it) {
        Vertex v = mesh->vertex(*it);

        for (Edge e : v.adjacentEdges()) {
            edges.insert(e.getIndex());
        }

        for (Face f : v.adjacentFaces()) {
            faces.insert(f.getIndex());
        }
    }

    for (std::set<size_t>::iterator it = subset.edges.begin(); it != subset.edges.end(); ++it) {
        Edge e = mesh->edge(*it);

        edges.insert(e.getIndex());

        for (Face f : e.adjacentFaces()) {
            faces.insert(f.getIndex());
        }
    }

    faces.insert(subset.faces.begin(), subset.faces.end());

    MeshSubset stared_subset(subset.vertices, edges, faces);
    return stared_subset;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    std::set<size_t> vertices = subset.vertices;
    std::set<size_t> edges = subset.edges;
    std::set<size_t> faces = subset.faces;

    for (std::set<size_t>::iterator it = faces.begin(); it != faces.end(); ++it) {
        Face f = mesh->face(*it);

        for (Edge e : f.adjacentEdges()) {
            edges.insert(e.getIndex());
        }
    }

    for (std::set<size_t>::iterator it = edges.begin(); it != edges.end(); ++it) {
        Edge e = mesh->edge(*it);

        for (Vertex v : e.adjacentVertices()) {
            vertices.insert(v.getIndex());
        }
    }

    return MeshSubset(vertices, edges, faces);
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset closured_stared = closure(star(subset));
    MeshSubset stared_closured = star(closure(subset));

    closured_stared.deleteSubset(stared_closured);

    return closured_stared;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    return subset.equals(closure(subset));
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    std::set<size_t> vertices;
    std::set<size_t> edges;
    std::set<size_t> faces;
    int degree;
    MeshSubset closured_subset;

    if (!subset.faces.empty()) {
        faces.insert(subset.faces.begin(), subset.faces.end());
        degree = 2;
    } else if (!subset.edges.empty()) {
        edges.insert(subset.edges.begin(), subset.edges.end());
        degree = 1;
    } else {
        vertices.insert(subset.vertices.begin(), subset.vertices.end());
        degree = 0;
    }

    closured_subset = MeshSubset(vertices, edges, faces);
    bool is_pure_complex = subset.equals(closure(closured_subset));

    return !is_pure_complex ? -1 : degree;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    std::set<size_t> faces;
    std::set<size_t> edges;
    std::set<size_t> vertices;

    if (!subset.faces.empty()) {
        std::unordered_map<size_t, int> edges_map;
        for (std::set<size_t>::iterator it = subset.faces.begin(); it != subset.faces.end(); ++it) {
            Face f = mesh->face(*it);

            for (Edge e : f.adjacentEdges()) {
                if (edges_map.count(e.getIndex()) > 0) {
                    edges_map[e.getIndex()]++;
                } else {
                    edges_map[e.getIndex()] = 1;
                }
            }
        }

        for (auto elem : edges_map) {
            if (elem.second == 1) {
                Edge e = mesh->edge(elem.first);
                edges.insert(elem.first);
                vertices.insert(e.firstVertex().getIndex());
                vertices.insert(e.secondVertex().getIndex());
            }
        }
    }

    if (!subset.edges.empty()) {
        std::unordered_map<size_t, int> vertices_map;
        for (std::set<size_t>::iterator it = subset.edges.begin(); it != subset.edges.end(); ++it) {
            Edge e = mesh->edge(*it);

            for (Vertex v : e.adjacentVertices()) {
                if (vertices_map.count(v.getIndex()) > 0) {
                    vertices_map[v.getIndex()]++;
                } else {
                    vertices_map[v.getIndex()] = 1;
                }
            }
        }

        for (auto elem : vertices_map) {
            if (elem.second == 1) {
                vertices.insert(elem.first);
            }
        }
    }

    return MeshSubset(vertices, edges, faces);
}