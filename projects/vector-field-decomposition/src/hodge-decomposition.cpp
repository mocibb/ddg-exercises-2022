// Implement member functions for HodgeDecomposition class.
#include "hodge-decomposition.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // build DEC operators
    // primal 1-form to dual 1-form
    this->hodge1 = geometry->buildHodgeStar1Form();
    // primal 2-form to dual 0-form
    this->hodge2 = geometry->buildHodgeStar2Form();
    // primal 0-form to primal 1-form
    this->d0 = geometry->buildExteriorDerivative0Form();
    // primal 1-form to primal 2-form
    this->d1 = geometry->buildExteriorDerivative1Form();

    // Build operator inverses.
    // Hint: Use the sparseInverseDiagonal() in utils/src/solvers.cpp to invert sparse diagonal matrices.
    // dual 1-form to primal 1-form
    this->hodge1Inv = sparseInverseDiagonal(this->hodge1);
    // dual 2-form to primal 0-form
    this->hodge2Inv = sparseInverseDiagonal(this->hodge2);
    // dual 1-form to dual 2-form
    this->d0T = this->d0.transpose();
    // dual 0-form to dual 1-form
    this->d1T = this->d1.transpose();
    
    // Construct 0-form Laplace matrix.
    // Shift matrix by a small constant (1e-8) to make it positive definite.
    this->A = d0T*hodge1*d0 + identityMatrix<double>(mesh->nVertices())*1e-8;

    // Construct 2-form matrix.
    // 输入为dual 0-form，去掉星是为了保持对称，可以使用PositiveDefiniteSolver
    this->B = d1*hodge1Inv*d1T;
}

/*
 * Compute the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The exact component dα of ω.
 */
Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {
    SparseMatrix<double> A = this->A;
    PositiveDefiniteSolver<double> solver(A);
    Vector<double> alpha = solver.solve(d0T*hodge1*omega);
    return this->d0*alpha;
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {
    SparseMatrix<double> B = this->B;
    PositiveDefiniteSolver<double> solver(B);
    // B = d*d 求出来的是*beta
    Vector<double> beta = solver.solve(d1*omega);
    return hodge1Inv*d1T*beta;
}

/*
 * Compute the harmonic component γ = ω - dα - 𝛿β of ω.
 *
 * Input: A primal 1-form <omega> on the edges of the input mesh, the exact component <dAlpha> of ω, and the coexact
 * component <deltaBeta> of ω.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                                            const Vector<double>& deltaBeta) const {

    
    return omega - dAlpha - deltaBeta;
}
