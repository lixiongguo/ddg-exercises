#include "trivial-connections.h"

#include "geometrycentral/numerical/linear_solvers.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

namespace {

double wrapPi(double angle) {
    const double turn = 2.0 * PI;
    while (angle > PI) angle -= turn;
    while (angle <= -PI) angle += turn;
    return angle;
}

double edgeSign(Halfedge he) {
    return he == he.edge().halfedge() ? 1.0 : -1.0;
}

} // namespace

TrivialConnections::TrivialConnections(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
    mesh = inputMesh;
    geometry = inputGeo;
    geometry->requireVertexGaussianCurvatures();

    this->hodge1 = geometry->buildHodgeStar1Form();
    this->d0 = geometry->buildExteriorDerivative0Form();
    SparseMatrix<double> hodgeD0 = this->hodge1 * this->d0;
    this->A = this->d0.transpose() * hodgeD0;
    shiftDiagonal(this->A, 1e-12);

    if (mesh->genus() > 0) {
        TreeCotree tree(mesh, geometry);
        tree.buildGenerators();
        HodgeDecomposition decomposition(mesh, geometry);
        HarmonicBases harmonic(mesh, geometry);
        this->bases = harmonic.compute(tree.generators, decomposition);
    }
    this->P = this->buildPeriodMatrix();
}

SparseMatrix<double> TrivialConnections::buildPeriodMatrix() const {
    const size_t count = this->bases.size();
    if (count == 0) return SparseMatrix<double>(0, 0);

    TreeCotree tree(mesh, geometry);
    tree.buildGenerators();
    Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(count), static_cast<Eigen::Index>(count));
    const size_t loops = std::min(count, tree.generators.size());
    for (size_t i = 0; i < loops; ++i) {
        for (Halfedge he : tree.generators[i]) {
            const double sign = edgeSign(he);
            const size_t edge = he.edge().getIndex();
            for (size_t j = 0; j < count; ++j) {
                if (this->bases[j].rows() <= static_cast<Eigen::Index>(edge)) continue;
                dense(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) += sign * this->bases[j][edge];
            }
        }
    }
    return dense.sparseView();
}

bool TrivialConnections::satsifyGaussBonnet(const Vector<double>& singularity) const {
    return (abs(singularity.sum() - geometry->eulerCharacteristic()) < 1e-8);
}

Vector<double> TrivialConnections::computeCoExactComponent(const Vector<double>& singularity) const {
    Vector<double> rhs = Vector<double>::Zero(mesh->nVertices());
    for (Vertex vertex : mesh->vertices()) {
        const size_t index = vertex.getIndex();
        const double indexValue = index < static_cast<size_t>(singularity.rows()) ? singularity[index] : 0.0;
        rhs[index] = -geometry->vertexGaussianCurvatures[vertex] + 2.0 * PI * indexValue;
    }
    SparseMatrix<double> system = this->A;
    const Vector<double> potential = solvePositiveDefinite(system, rhs);
    return this->hodge1 * this->d0 * potential;
}

double TrivialConnections::transportNoRotation(Halfedge he, double alphaI) const {
    Vector3 u = geometry->halfedgeVector(he);

    Face fi = he.face();
    Face fj = he.twin().face();
    Vector3 e1 = geometry->halfedgeVector(fi.halfedge()).normalize();
    Vector3 e2 = cross(geometry->faceNormal(fi), e1);
    Vector3 f1 = geometry->halfedgeVector(fj.halfedge()).normalize();
    Vector3 f2 = cross(geometry->faceNormal(fj), f1);
    double thetaIJ = atan2(dot(u, e2), dot(u, e1));
    double thetaJI = atan2(dot(u, f2), dot(u, f1));

    return alphaI - thetaIJ + thetaJI;
}

Vector<double> TrivialConnections::computeHarmonicComponent(const Vector<double>& deltaBeta) const {
    Vector<double> gamma = Vector<double>::Zero(mesh->nEdges());
    if (this->bases.empty() || this->P.rows() == 0) return gamma;

    TreeCotree tree(mesh, geometry);
    tree.buildGenerators();
    const size_t count = this->bases.size();
    Vector<double> defect = Vector<double>::Zero(count);
    const size_t loops = std::min(count, tree.generators.size());
    for (size_t i = 0; i < loops; ++i) {
        double sum = 0.0;
        for (Halfedge he : tree.generators[i]) {
            const size_t edge = he.edge().getIndex();
            if (static_cast<Eigen::Index>(edge) >= deltaBeta.rows()) continue;
            sum += edgeSign(he) * deltaBeta[edge];
        }
        defect[i] = -wrapPi(sum);
    }

    Eigen::MatrixXd dense = Eigen::MatrixXd(this->P);
    const Eigen::VectorXd coefficients = dense.fullPivLu().solve(defect);
    for (size_t j = 0; j < count; ++j) {
        if (this->bases[j].rows() != gamma.rows()) continue;
        gamma += coefficients[static_cast<Eigen::Index>(j)] * this->bases[j];
    }
    return gamma;
}

Vector<double> TrivialConnections::computeConnections(const Vector<double>& singularity) const {
    if (!this->satsifyGaussBonnet(singularity)) {
        std::cerr << "Singularities do not add up to the Euler characteristic of the mesh" << std::endl;
        return Vector<double>::Zero(mesh->nEdges());
    }
    const Vector<double> coexact = this->computeCoExactComponent(singularity);
    const Vector<double> harmonic = this->computeHarmonicComponent(coexact);
    if (harmonic.rows() != coexact.rows()) return coexact;
    return coexact + harmonic;
}
