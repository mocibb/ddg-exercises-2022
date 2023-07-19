// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total +=  edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {
    Vector3 pB = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pA = inputVertexPositions[he.vertex()];

    Vector3 vecR = pB - pA;
    Vector3 vecL = pC - pA;

    double cotValue = dot(vecR, vecL) / norm(cross(vecR, vecL));
    return cotValue;
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {
    double area = 0;
    for (Face f : v.adjacentFaces()) {
        area += faceArea(f);
    }

    return area / 3;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {
    Halfedge heA = c.halfedge();
    Halfedge heOpp = heA.next();
    Halfedge heB = heOpp.next();

    double lOpp = norm(inputVertexPositions[heOpp.tipVertex()] - inputVertexPositions[heOpp.tailVertex()]);
    double lA = norm(inputVertexPositions[heA.tipVertex()] - inputVertexPositions[heA.tailVertex()]);
    double lB = norm(inputVertexPositions[heB.tipVertex()] - inputVertexPositions[heB.tailVertex()]);

    double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
    q = clamp(q, -1.0, 1.0);

    return std::acos(q);
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
    Face f1 = he.face();
    Face f2 = he.twin().face();

    Halfedge f1h1 = f1.halfedge();
    Halfedge f1h2 = f1h1.next();
    Halfedge f2h1 = f2.halfedge();
    Halfedge f2h2 = f2h1.next();

    Vector3 N1 = cross(inputVertexPositions[f1h1.tipVertex()] - inputVertexPositions[f1h1.tailVertex()],
                       inputVertexPositions[f1h2.tipVertex()] - inputVertexPositions[f1h2.tailVertex()]);
    N1 = unit(N1);
    Vector3 N2 = cross(inputVertexPositions[f2h1.tipVertex()] - inputVertexPositions[f2h1.tailVertex()],
                       inputVertexPositions[f2h2.tipVertex()] - inputVertexPositions[f2h2.tailVertex()]);
    N2 = unit(N2);
    Vector3 pTail = inputVertexPositions[he.vertex()];
    Vector3 pTip = inputVertexPositions[he.next().vertex()];
    Vector3 edgeDir = unit(pTip - pTail);

    return atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {
    Vector3 normal = {0, 0, 0};
    for (Face f : v.adjacentFaces()) {
        Halfedge fh1 = f.halfedge();
        Halfedge fh2 = fh1.next();
        Vector3 N = cross(inputVertexPositions[fh1.tipVertex()] - inputVertexPositions[fh1.tailVertex()],
                          inputVertexPositions[fh2.tipVertex()] - inputVertexPositions[fh2.tailVertex()]);
        N = unit(N);        
        normal += N;
    }
    return unit(normal);
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
    Vector3 normal = {0, 0, 0};

    for (Corner c : v.adjacentCorners()) {
        Face f = c.face();
        Halfedge fh1 = f.halfedge();
        Halfedge fh2 = fh1.next();
        Vector3 N = cross(inputVertexPositions[fh1.tipVertex()] - inputVertexPositions[fh1.tailVertex()],
                          inputVertexPositions[fh2.tipVertex()] - inputVertexPositions[fh2.tailVertex()]);
        N = unit(N);
        normal += angle(c) * N;
    }
    return unit(normal);
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {
    Vector3 normal = {0, 0, 0}; 

    for (Face f : v.adjacentFaces()) {
        Halfedge he = f.halfedge();
        Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
        
        Halfedge he2 = he.next().next();
        Vector3 e_ij1 = inputVertexPositions[he2.tipVertex()] - inputVertexPositions[he.tailVertex()];
        // std::cout << inputVertexPositions[he.tipVertex()] << ", "
        //           << inputVertexPositions[he.tailVertex()] << ", "
        //           << inputVertexPositions[he2.tipVertex()] << std::endl;
        double l_ij = norm(e_ij);
        double l_ij1 = norm(e_ij1);
        normal += cross(e_ij, e_ij1) / (l_ij * l_ij * l_ij1 * l_ij1);

    }
    // std::cout << normal << std::endl;
    return unit(normal);
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {
    Vector3 normal = {0, 0, 0};

    for (Corner c : v.adjacentCorners()) {
        Face f = c.face();
        Halfedge fh1 = f.halfedge();
        Halfedge fh2 = fh1.next();
        Vector3 N = cross(inputVertexPositions[fh1.tipVertex()] - inputVertexPositions[fh1.tailVertex()],
                          inputVertexPositions[fh2.tipVertex()] - inputVertexPositions[fh2.tailVertex()]);
        N = unit(N);
        normal += faceArea(f) * N;
    }
    return unit(normal);
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {
    Vector3 normal = {0, 0, 0};

    for (Corner c : v.adjacentCorners()) {
        Halfedge he = c.halfedge();
        normal +=
            dihedralAngle(he) * unit(inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()]);
    }
    return unit(normal);
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {
    Vector3 normal = {0, 0, 0};

    //for (Corner c : v.adjacentCorners()) {
    for (Halfedge he : v.outgoingHalfedges()) {        
      //  Halfedge he = c.halfedge();

        Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];

        normal += 0.5*(cotan(he)+cotan(he.twin()))*e_ij;
    }
    return unit(normal);
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {
    double sum = 0;
    for (Corner c : v.adjacentCorners()) {
        sum += angle(c);
    }
    // TODO
    return 2 * PI - sum; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {
    double all_sum = 0;
    for (Vertex v : mesh.vertices()) {
        all_sum += angleDefect(v);
    }
    return all_sum;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
    double curvature = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];    
        double phi = dihedralAngle(he);
        curvature += phi * norm(e_ij);
    }

    return 0.5*curvature;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
    double area = 0;
    for (Corner c : v.adjacentCorners()) {
        Halfedge he = c.halfedge();
        Halfedge he2 = he.next().next();
        double l_ij = norm(inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()]);
        double l_ij2 = norm(inputVertexPositions[he2.tipVertex()] - inputVertexPositions[he2.tailVertex()]);    
        area += (cotan(he2)*l_ij2*l_ij2 + cotan(he)*l_ij*l_ij);
    }

    return 1./8*area;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {
    double A = circumcentricDualArea(v);
    double H = scalarMeanCurvature(v) / A;
    double K = angleDefect(v) / A;
    double delta = sqrt(H*H-K);
    // TODO
    return std::make_pair(H-delta, H+delta); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {
    SparseMatrix<double> d0(buildExteriorDerivative0Form());
    SparseMatrix<double> star1(buildHodgeStar1Form());
    return d0.transpose()*star1*d0 + identityMatrix<double>(mesh.nVertices())*1e-8;
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {
    return buildHodgeStar0Form();
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {
    SparseMatrix<std::complex<double>> A = identityMatrix<std::complex<double>>(mesh.nVertices())*1e-8;
    
    // TODO
    return A; // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral
