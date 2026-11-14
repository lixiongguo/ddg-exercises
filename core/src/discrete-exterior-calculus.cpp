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
    return identityMatrix<double>(mesh.nVertices());
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    SparseMatrix<double> star(mesh.nEdges(), mesh.nEdges());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.nEdges());
    for (Edge edge : mesh.edges()) {
        triplets.emplace_back(edge.getIndex(), edge.getIndex(), edgeCotanWeight(edge));
    }
    star.setFromTriplets(triplets.begin(), triplets.end());
    return star;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    SparseMatrix<double> star(mesh.nFaces(), mesh.nFaces());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.nFaces());
    for (Face face : mesh.faces()) {
        const double area = faceArea(face);
        triplets.emplace_back(face.getIndex(), face.getIndex(), area > 0.0 ? 1.0 / area : 0.0);
    }
    star.setFromTriplets(triplets.begin(), triplets.end());
    return star;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    SparseMatrix<double> derivative(mesh.nEdges(), mesh.nVertices());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.nEdges() * 2);
    for (Edge edge : mesh.edges()) {
        const Halfedge he = edge.halfedge();
        const size_t index = edge.getIndex();
        triplets.emplace_back(index, he.tipVertex().getIndex(), 1.0);
        triplets.emplace_back(index, he.tailVertex().getIndex(), -1.0);
    }
    derivative.setFromTriplets(triplets.begin(), triplets.end());
    return derivative;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    SparseMatrix<double> derivative(mesh.nFaces(), mesh.nEdges());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.nFaces() * 3);
    for (Face face : mesh.faces()) {
        for (Halfedge he : face.adjacentHalfedges()) {
            const double sign = (he == he.edge().halfedge()) ? 1.0 : -1.0;
            triplets.emplace_back(face.getIndex(), he.edge().getIndex(), sign);
        }
    }
    derivative.setFromTriplets(triplets.begin(), triplets.end());
    return derivative;
}

} // namespace surface
} // namespace geometrycentral