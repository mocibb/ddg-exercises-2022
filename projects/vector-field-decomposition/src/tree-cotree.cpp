// Implement member functions for TreeCotree class.
#include "tree-cotree.h"
#include <queue>

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TreeCotree::TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a primal spanning tree on a mesh without boundary. More specifically, populate the member variable
 * <vertexParent>, which is a std::map that maps each vertex of the input mesh to its parent in the primal spanning
 * tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildPrimalSpanningTree() {
    std::queue<Vertex> queue;
    Vertex root = mesh->vertex(randomInt(0, mesh->nVertices()));

    // 从某个顶点开始
    queue.push(root);
    vertexParent.emplace(root, root);
    std::cout << "root = " << root.getIndex() << std::endl;

    while(!queue.empty()) {
        auto v = queue.front();
        queue.pop();

        for (Halfedge he : v.outgoingHalfedges()) { 
            // do not cross edges in T*
            if (inDualSpanningCotree(he)) {
                continue;
            }       
            Vertex w = he.tipVertex();
            if (!vertexParent.count(w)) {
                vertexParent.emplace(w, v);
                queue.push(w);
            }
        }
    }
}

/*
 * Check whether a halfedge is in the primal spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the primal spanning tree, false otherwise.
 */
bool TreeCotree::inPrimalSpanningTree(Halfedge he) {
    Vertex v = he.tailVertex();
    Vertex w = he.tipVertex();

    if (vertexParent.count(v)) {
        if (vertexParent[v] == w) {
            return true;
        }
    }

    if (vertexParent.count(w)) {
        if (vertexParent[w] == v) {
            return true;
        }
    } 

    return false;
}

/*
 * Build a dual spanning tree on a mesh without boundary. More specificially, populate the member variable <faceParent>,
 * which is a std::map that maps each face of the input mesh to its parent in the dual spanning tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildDualSpanningCoTree() {
    std::queue<Face> queue;
    Face root = mesh->face(randomInt(0, mesh->nFaces()));

    // 从某个顶点开始
    queue.push(root);
    faceParent[root] = root;

    while(!queue.empty()) {
        auto f = queue.front();
        queue.pop();

        for (Halfedge he : f.adjacentHalfedges()) {
            //不考虑带边界情况
            Face g = he.twin().face();

            if (!faceParent.count(g)) {
                faceParent[g] = f;
                queue.push(g);
            }
        }
    }
}

/*
 * Check whether a halfedge is in the dual spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the dual spanning tree, false otherwise.
 */
bool TreeCotree::inDualSpanningCotree(Halfedge he) {
    Face f = he.face();
    Face g = he.twin().face();

    if (faceParent.count(f)) {
        if (faceParent[f] == g) {
            return true;
        }
    }

    if (faceParent.count(g)) {
        if (faceParent[g] == f) {
            return true;
        }
    }    
    return false;
}

/*
 * Returns a halfedge lying on the shared edge between face f and g.
 *
 * Input: Two adjacent faces <f> and <g>.
 * Returns: A halfedge lying on the shared edge between face f and g.
 */
Halfedge TreeCotree::sharedHalfedge(Face f, Face g) const {

    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // Should never get here!
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

/*
 * Compute the homology generators of the mesh.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildGenerators() {
    // 1. 计算对偶Mesh的Spanning Tree
    buildDualSpanningCoTree();
    buildPrimalSpanningTree();


    // 2. Build generators and populate this->generators
    // loop throught each dual edge    
    for (const auto& e : mesh->edges()) {
        Halfedge he = e.halfedge();
        if (!inDualSpanningCotree(he) && !inPrimalSpanningTree(he)) {
            Vertex v1 = he.tipVertex();
            Vertex v2 = he.tailVertex();

            generators.push_back(std::vector<Halfedge>());
            // v2 -> v1
            generators.back().push_back(he);
            
            // v1 -> root
            Vertex v1p = vertexParent.at(v1);
            while(v1p != v1) {
                for (Halfedge he1 : v1.outgoingHalfedges()) {
                    if (he1.tipVertex() == v1p) {
                        generators.back().push_back(he1);
                        break;
                    }
                }
                v1 = v1p;
                v1p = vertexParent.at(v1);
            }
            // root -> v2
            std::vector<Halfedge> backward;
            Vertex v2p = vertexParent.at(v2);
            while(v2p != v2) {
                for (Halfedge he2 : v2.incomingHalfedges()) {
                    if (he2.tailVertex() == v2p) {
                        backward.push_back(he2);
                        break;
                    }
                }
                v2 = v2p;
                v2p = vertexParent.at(v2);
            }
            
            generators.back().insert(generators.back().end(), backward.rbegin(), backward.rend());
        }
    }

    std::cout << "generators.size() = " << generators.size() << std::endl;

}
