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
    Vertex root = *(mesh->vertices().begin());

    // 从某个顶点开始
    queue.push(root);
    vertexParent.emplace(root, root);

    while(!queue.empty()) {
        auto v = queue.front();
        queue.pop();

        for (Halfedge he : v.outgoingHalfedges()) { 
            // do not cross edges in T*
            // if (inDualSpanningCotree(he)) {
            //     continue;
            // }       
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
    Face root = *(mesh->faces().begin());

    // 从某个顶点开始
    queue.push(root);
    faceParent[root] = root;

    while(!queue.empty()) {
        auto f = queue.front();
        queue.pop();       

        for (Halfedge he : f.adjacentHalfedges()) {
            if (inPrimalSpanningTree(he)) {
                continue;
            }
            
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
    buildPrimalSpanningTree();
    buildDualSpanningCoTree();


    // 2. Build generators and populate this->generators
    // loop throught each dual edge    
    // 计算对偶网络的生成元
    for (const auto& e : mesh->edges()) {
        Halfedge he = e.halfedge();
        //if (!inDualSpanningCotree(he) && !inPrimalSpanningTree(he)) {
        if (!inDualSpanningCotree(he) && !inPrimalSpanningTree(he)) {
            Face f1 = he.face();
            Face f2 = he.twin().face();

            generators.push_back(std::vector<Halfedge>());
            std::vector<Halfedge>& generator = generators.back();
            // f2 -> f1 -> root -> f2
            // f2 -> f1
            generator.push_back(he);
            
            // f1 -> root
            Face f1p = faceParent.at(f1);
            while(f1p != f1) {
                generator.push_back(sharedHalfedge(f1p, f1));
                f1 = f1p;
                f1p = faceParent.at(f1);
            }
            
            // root -> f2
            std::vector<Halfedge> backward;
            Face f2p = faceParent.at(f2);
            while(f2p != f2) {
                Halfedge he1 = sharedHalfedge(f2p, f2);
                // 从f2到root时可能提前遇到loop点
                auto it = std::find(generator.begin(), generator.end(), he1);
                if (it != generator.end()) {
                    generator.erase(it, generator.end());
                    break;
                }
                backward.push_back(sharedHalfedge(f2, f2p));
                f2 = f2p;
                f2p = faceParent.at(f2);
            }
            
            generator.insert(generator.end(), backward.rbegin(), backward.rend());
            
        }
    }
}
