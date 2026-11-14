#include "harmonic-bases.h"

HarmonicBases::HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
    mesh = inputMesh;
    geometry = inputGeo;
}

Vector<double> HarmonicBases::buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const {
    Vector<double> omega = Vector<double>::Zero(mesh->nEdges());
    for (Halfedge he : generator) {
        const Edge edge = he.edge();
        const double sign = (he == edge.halfedge()) ? 1.0 : -1.0;
        omega[edge.getIndex()] += sign;
    }
    return omega;
}

std::vector<Vector<double>> HarmonicBases::compute(const std::vector<std::vector<Halfedge>>& generators,
                                                   const HodgeDecomposition& hodgeDecomposition) const {
    std::vector<Vector<double>> gammas;
    gammas.reserve(generators.size());
    for (const std::vector<Halfedge>& generator : generators) {
        const Vector<double> omega = buildClosedPrimalOneForm(generator);
        const Vector<double> exact = hodgeDecomposition.computeExactComponent(omega);
        const Vector<double> coexact = hodgeDecomposition.computeCoExactComponent(omega);
        gammas.push_back(hodgeDecomposition.computeHarmonicComponent(omega, exact, coexact));
    }
    return gammas;
}
