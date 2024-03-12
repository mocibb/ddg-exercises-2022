// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <unordered_map>
#include <chrono>

// #define USE_SPARSE_MATRIX
#define MEASURE_PERFORMANCE

using namespace geometrycentral;
using namespace geometrycentral::surface;

#ifdef MEASURE_PERFORMANCE
    void report_elapsed_time(const std::chrono::time_point<std::chrono::steady_clock>& start, const std::string& msg) {
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
        std::cout << msg << " = " << elapsed.count() << "ms\n";
    }
#else
    void report_elapsed_time(std::chrono::time_point<clock_t, duration_t> const& start, std::string msg) {

    }
#endif

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
    // 随机设置index
    auto perm = [](int N) {
        std::vector<int> v(N);
        std::iota(std::begin(v), std::end(v), 0);
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle (v.begin(), v.end(), std::default_random_engine(seed));

        return v;
    };

    size_t idx = 0;
    auto vIndex = perm(mesh->nVertices());
    int i = 0;
    for (Vertex v : mesh->vertices()) {
        // geometry->vertexIndices[v] = vIndex[i++];
    }

    auto eIndex = perm(mesh->nEdges());
    i = 0;
    for (Edge e : mesh->edges()) {
        // geometry->edgeIndices[e] = eIndex[i++];
    }

    auto fIndex = perm(mesh->nFaces());
    i = 0;
    for (Face f : mesh->faces()) {
        // geometry->faceIndices[f] = fIndex[i++];
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
    tl.reserve(mesh->nFaces() * 3);

    size_t idx = 0;
    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
        tl.push_back(T(idx, geometry->edgeIndices[f.halfedge().edge()], 1));
        tl.push_back(T(idx, geometry->edgeIndices[f.halfedge().next().edge()], 1));
        tl.push_back(T(idx, geometry->edgeIndices[f.halfedge().next().next().edge()], 1));
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
    for (size_t idx : subset.vertices) {
        vertices[idx] = 1;
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
    for (size_t idx : subset.edges) {
        edges[idx] = 1;
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
    for (size_t idx : subset.faces) {
        faces[idx] = 1;
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
    auto start = std::chrono::steady_clock::now();
#ifdef USE_SPARSE_MATRIX
    MeshSubset stared_subset = subset.deepCopy();
    stared_subset.addEdges(copyFromVector(this->A0*this->buildVertexVector(stared_subset)));
    stared_subset.addFaces(copyFromVector(this->A1*this->buildEdgeVector(stared_subset)));

    report_elapsed_time(start, "star");
    return stared_subset;
#else
    std::set<size_t> vertices = subset.vertices;
    std::set<size_t> edges = subset.edges;
    std::set<size_t> faces = subset.faces;

    for (size_t vid : vertices) {
        const Vertex v = mesh->vertex(vid);

        for (Edge e : v.adjacentEdges()) {
            edges.insert(e.getIndex());
        }
    }

    for (size_t eid : edges) {
        const Edge e = mesh->edge(eid);

        for (Face f : e.adjacentFaces()) {
            faces.insert(f.getIndex());
        }
    }

    report_elapsed_time(start, "star");
    return MeshSubset(vertices, edges, faces);
#endif
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    auto start = std::chrono::steady_clock::now();
#ifdef USE_SPARSE_MATRIX   
    MeshSubset closured_subset = subset.deepCopy();
    closured_subset.addEdges(copyFromVector(this->A1.transpose()*this->buildFaceVector(closured_subset)));
    closured_subset.addVertices(copyFromVector(this->A0.transpose()*this->buildEdgeVector(closured_subset)));

    report_elapsed_time(start, "closure");
    return closured_subset;
#else 
    std::set<size_t> vertices = subset.vertices;
    std::set<size_t> edges = subset.edges;
    std::set<size_t> faces = subset.faces;

    for (size_t fid : faces) {
        const Face f = mesh->face(fid);

        for (Edge e : f.adjacentEdges()) {
            edges.insert(e.getIndex());
        }
    }

    for (size_t eid : edges) {
        const Edge e = mesh->edge(eid);

        for (Vertex v : e.adjacentVertices()) {
            vertices.insert(v.getIndex());
        }
    }

    report_elapsed_time(start, "closure");
    return MeshSubset(vertices, edges, faces);
#endif
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    auto start = std::chrono::steady_clock::now();
    MeshSubset closured_stared = closure(star(subset));
    MeshSubset stared_closured = star(closure(subset));

    closured_stared.deleteSubset(stared_closured);

    report_elapsed_time(start, "link");
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
    auto start = std::chrono::steady_clock::now();
#ifdef USE_SPARSE_MATRIX  
    MeshSubset boundary_subset;

    if (!subset.faces.empty()) {
        Vector<size_t> edgeVector = this->A1.transpose()*this->buildFaceVector(subset);
        for (size_t idx = 0; idx < (size_t)edgeVector.size(); idx++) {
            if (edgeVector[idx] == 1) {
                boundary_subset.addEdge(idx);
                Edge e = mesh->edge(idx);
                boundary_subset.addVertex(e.firstVertex().getIndex());
                boundary_subset.addVertex(e.secondVertex().getIndex());
            }
        }
    }

    if (!subset.edges.empty()) {
        Vector<size_t> vertexVector = this->A0.transpose()*this->buildEdgeVector(subset);
        for (size_t idx = 0; idx < (size_t)vertexVector.size(); idx++) {
            if (vertexVector[idx] == 1) {
                boundary_subset.addVertex(idx);
            }
        }
    }

    report_elapsed_time(start, "boundary");
    return boundary_subset;
#else    
    std::set<size_t> faces;
    std::set<size_t> edges;
    std::set<size_t> vertices;

    if (!subset.faces.empty()) {
        std::unordered_map<size_t, int> edges_map;
        
        for (size_t fid : subset.faces) {
            const Face f = mesh->face(fid);

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
    report_elapsed_time(start, "boundary");
    return MeshSubset(vertices, edges, faces);
#endif
}
