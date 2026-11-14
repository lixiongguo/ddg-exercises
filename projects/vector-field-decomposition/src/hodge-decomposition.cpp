#include "hodge-decomposition.h"

#include "geometrycentral/numerical/linear_solvers.h"

HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
    mesh = inputMesh;
    geometry = inputGeo;

    this->hodge1 = geometry->buildHodgeStar1Form();
    this->hodge2 = geometry->buildHodgeStar2Form();
    this->d0 = geometry->buildExteriorDerivative0Form();
    this->d1 = geometry->buildExteriorDerivative1Form();
    this->hodge1Inv = sparseInverseDiagonal(this->hodge1);
    this->hodge2Inv = sparseInverseDiagonal(this->hodge2);
    this->d0T = this->d0.transpose();
    this->d1T = this->d1.transpose();

    SparseMatrix<double> hodge1D0 = this->hodge1 * this->d0;
    this->A = this->d0T * hodge1D0;
    shiftDiagonal(this->A, 1e-8);

    SparseMatrix<double> hodgeInvD1T = this->hodge1Inv * this->d1T;
    this->B = this->d1 * hodgeInvD1T;
    shiftDiagonal(this->B, 1e-12);
}

Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {
    if (omega.rows() != this->d0.rows()) return Vector<double>::Zero(mesh->nEdges());
    Vector<double> rhs = this->d0T * this->hodge1 * omega;
    SparseMatrix<double> system = this->A;
    Vector<double> alpha = solvePositiveDefinite(system, rhs);
    return this->d0 * alpha;
}

Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {
    if (omega.rows() != this->d1.cols()) return Vector<double>::Zero(mesh->nEdges());
    Vector<double> rhs = this->d1 * omega;
    SparseMatrix<double> system = this->B;
    Vector<double> dualPotential = solvePositiveDefinite(system, rhs);
    return this->hodge1Inv * this->d1T * dualPotential;
}

Vector<double> HodgeDecomposition::computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                                            const Vector<double>& deltaBeta) const {
    if (omega.rows() != dAlpha.rows() || omega.rows() != deltaBeta.rows()) {
        return Vector<double>::Zero(mesh->nEdges());
    }
    return omega - dAlpha - deltaBeta;
}
