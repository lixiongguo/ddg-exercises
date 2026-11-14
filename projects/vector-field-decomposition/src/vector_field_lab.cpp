#include "vector_field_lab.h"

#include "harmonic-bases.h"
#include "hodge-decomposition.h"
#include "trivial-connections.h"
#include "tree-cotree.h"

#include "geometrycentral/surface/meshio.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

struct FieldProcess::State {
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    TreeCotree tree;
    HodgeDecomposition hodge;
    HarmonicBases harmonic;
    TrivialConnections connections;
    std::vector<Vector<double>> bases;
    Vector<double> omega;
    Vector<double> exact;
    Vector<double> coexact;
    Vector<double> harmonicForm;
    Vector<double> connection;
    Vector<double> singularity;
};

namespace {

Vector3 faceCentroid(VertexPositionGeometry* geometry, Face face) {
    Halfedge he = face.halfedge();
    Vector3 a = geometry->inputVertexPositions[he.vertex()];
    Vector3 b = geometry->inputVertexPositions[he.next().vertex()];
    Vector3 c = geometry->inputVertexPositions[he.next().next().vertex()];
    return (a + b + c) / 3.0;
}

void enableOneForm(polyscope::SurfaceVectorQuantity* quantity, bool enabled) {
    if (quantity == nullptr) return;
    quantity->setEnabled(enabled);
    quantity->setVectorColor({0.05, 0.05, 0.05});
}

} // namespace

FieldProcess& fieldProcess() {
    static FieldProcess process;
    return process;
}

double FieldProcess::singularitySum() const {
    if (state_ == nullptr) return 0.0;
    return state_->singularity.sum();
}

double* FieldProcess::singularityAt(int vertex) {
    if (state_ == nullptr || vertex < 0 || vertex >= state_->singularity.rows()) return nullptr;
    return &state_->singularity[vertex];
}

int FieldProcess::basisCount() const {
    if (state_ == nullptr) return 0;
    return static_cast<int>(state_->bases.size());
}

bool FieldProcess::prepare(const std::string& path, std::string& error) {
    if (ready_ && path_ == path && state_ != nullptr) return true;
    delete state_;
    state_ = nullptr;
    ready_ = false;
    drawn_ = nullptr;
    connectionForm_ = nullptr;
    try {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = readManifoldSurfaceMesh(path);
        if (mesh->nBoundaryLoops() > 0) {
            error = "生成树、Hodge 1-形式和奇异点联络需要闭网格";
            return false;
        }
        state_ = new State();
        state_->mesh = std::move(mesh);
        state_->geometry = std::move(geometry);
        state_->geometry->requireVertexGaussianCurvatures();
        state_->tree = TreeCotree(state_->mesh.get(), state_->geometry.get());
        state_->tree.buildGenerators();
        state_->hodge = HodgeDecomposition(state_->mesh.get(), state_->geometry.get());
        state_->harmonic = HarmonicBases(state_->mesh.get(), state_->geometry.get());
        state_->bases = state_->harmonic.compute(state_->tree.generators, state_->hodge);
        state_->connections = TrivialConnections(state_->mesh.get(), state_->geometry.get());
        state_->singularity = Vector<double>::Zero(state_->mesh->nVertices());
        state_->connection = Vector<double>::Zero(state_->mesh->nEdges());
        lengthScale_ = state_->geometry->meanEdgeLength();
        euler_ = state_->geometry->eulerCharacteristic();
        basisIndex_ = 0;
        path_ = path;
        newOneForm();
        const int genus = state_->mesh->genus();
        std::ostringstream text;
        text << "亏格 " << genus << "    基本环 " << state_->tree.generators.size()
             << "    调和基 " << state_->bases.size();
        note_ = text.str();
        ready_ = true;
        return true;
    } catch (const std::exception& exception) {
        error = exception.what();
        delete state_;
        state_ = nullptr;
        return false;
    }
}

void FieldProcess::decompose() {
    if (state_ == nullptr) return;
    state_->exact = state_->hodge.computeExactComponent(state_->omega);
    state_->coexact = state_->hodge.computeCoExactComponent(state_->omega);
    state_->harmonicForm = state_->hodge.computeHarmonicComponent(state_->omega, state_->exact, state_->coexact);
    if (state_->exact.rows() != state_->omega.rows()) state_->exact = Vector<double>::Zero(state_->omega.rows());
    if (state_->coexact.rows() != state_->omega.rows()) state_->coexact = Vector<double>::Zero(state_->omega.rows());
    if (state_->harmonicForm.rows() != state_->omega.rows()) {
        state_->harmonicForm = Vector<double>::Zero(state_->omega.rows());
    }
}

void FieldProcess::newOneForm() {
    if (state_ == nullptr) return;
    const int vertices = static_cast<int>(state_->mesh->nVertices());
    const int faces = static_cast<int>(state_->mesh->nFaces());
    Vector<double> vertexNoise = Vector<double>::Random(vertices);
    Vector<double> faceNoise = Vector<double>::Random(faces);
    state_->omega = state_->hodge.d0 * vertexNoise + state_->hodge.hodge1Inv * state_->hodge.d1T * faceNoise;
    for (const Vector<double>& basis : state_->bases) {
        if (basis.rows() == state_->omega.rows()) state_->omega += basis;
    }
    decompose();
}

void FieldProcess::nextBasis() {
    if (basisCount() == 0) return;
    basisIndex_ = (basisIndex_ + 1) % basisCount();
}

void FieldProcess::resetSingularities() {
    if (state_ == nullptr) return;
    state_->singularity = Vector<double>::Zero(state_->mesh->nVertices());
    state_->connection = Vector<double>::Zero(state_->mesh->nEdges());
    polyscope::state::subset.vertices.clear();
    polyscope::state::currVertexIndex = -1;
    polyscope::state::deleteVertexIndex = -1;
}

bool FieldProcess::computeConnection(std::string& error) {
    if (state_ == nullptr) {
        error = "还没有闭网格";
        return false;
    }
    if (!state_->connections.satsifyGaussBonnet(state_->singularity)) {
        error = "奇异点指标之和要等于欧拉示性数 " + std::to_string(euler_);
        return false;
    }
    state_->connection = state_->connections.computeConnections(state_->singularity);
    if (state_->connection.rows() != static_cast<Eigen::Index>(state_->mesh->nEdges()) ||
        state_->connection.norm() == 0.0) {
        error = "联络没有解出来";
        return false;
    }
    return true;
}

void FieldProcess::suspendPicking() { picking_ = false; }

void FieldProcess::showSelected() {
    if (!picking_ || drawn_ == nullptr || state_ == nullptr) return;
    const int removed = polyscope::state::deleteVertexIndex;
    if (removed >= 0 && removed < state_->singularity.rows()) state_->singularity[removed] = 0.0;

    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edges;
    for (size_t index : polyscope::state::subset.vertices) {
        if (static_cast<Eigen::Index>(index) >= state_->singularity.rows()) continue;
        if (state_->singularity[index] == 0.0) state_->singularity[index] = 1.0;
        positions.push_back(state_->geometry->inputVertexPositions[index]);
    }
    polyscope::SurfaceGraphQuantity* selected = drawn_->addSurfaceGraphQuantity("选中顶点", positions, edges);
    selected->setEnabled(true);
    selected->setRadius(lengthScale_ * 0.2);
    selected->setColor({1.0, 0.65, 0.0});

    const int current = polyscope::state::currVertexIndex;
    std::vector<Vector3> currentPosition;
    if (current >= 0 && current < state_->singularity.rows()) {
        currentPosition.push_back(state_->geometry->inputVertexPositions[current]);
    }
    polyscope::SurfaceGraphQuantity* active =
        drawn_->addSurfaceGraphQuantity("当前顶点", currentPosition, std::vector<std::array<size_t, 2>>());
    active->setEnabled(!currentPosition.empty());
    active->setRadius(lengthScale_ * 0.2);
    active->setColor({1.0, 0.0, 0.0});
}

void FieldProcess::drawDirection(polyscope::SurfaceMesh* mesh) {
    if (state_ == nullptr || state_->connection.rows() != static_cast<Eigen::Index>(state_->mesh->nEdges())) return;
    ManifoldSurfaceMesh* surface = state_->mesh.get();
    VertexPositionGeometry* geometry = state_->geometry.get();
    std::map<Face, Face> parent;
    for (Face face : surface->faces()) parent[face] = face;
    const Face root = surface->face(0);
    FaceData<double> angle(*surface, 0.0);
    std::vector<Face> pending;
    pending.push_back(root);
    while (!pending.empty()) {
        const Face face = pending.back();
        pending.pop_back();
        for (Halfedge he : face.adjacentHalfedges()) {
            const Face neighbor = he.twin().face();
            if (parent[neighbor] == neighbor && neighbor != root) {
                const double sign = he.edge().halfedge() == he ? 1.0 : -1.0;
                parent[neighbor] = face;
                const double turn = sign * state_->connection[he.edge().getIndex()];
                angle[neighbor] = state_->connections.transportNoRotation(he, angle[face] - turn);
                pending.push_back(neighbor);
            }
        }
    }
    std::vector<Vector3> field(surface->nFaces());
    for (Face face : surface->faces()) {
        const Vector3 u{std::cos(angle[face]), std::sin(angle[face]), 0.0};
        const Vector3 e1 = geometry->halfedgeVector(face.halfedge()).normalize();
        const Vector3 e2 = cross(geometry->faceNormal(face), e1);
        field[face.getIndex()] = e1 * u[0] + e2 * u[1];
    }
    polyscope::SurfaceVectorQuantity* arrows = mesh->addFaceVectorQuantity("方向场", field);
    arrows->setVectorColor({0.05, 0.05, 0.05});
    arrows->setVectorLengthScale(lengthScale_ * 0.5);
    arrows->setEnabled(true);
}

void FieldProcess::draw(int process, bool resetCamera, std::string& status) {
    picking_ = false;
    connectionForm_ = nullptr;
    drawn_ = nullptr;
    if (!ready_ || state_ == nullptr) {
        status = note_.empty() ? "这一步还没有结果" : note_;
        return;
    }
    ManifoldSurfaceMesh* surface = state_->mesh.get();
    VertexPositionGeometry* geometry = state_->geometry.get();
    polyscope::removeAllStructures();
    polyscope::SurfaceMesh* mesh = polyscope::registerSurfaceMesh(
        "向量场", geometry->inputVertexPositions, surface->getFaceVertexList(), polyscopePermutations(*surface));
    mesh->setSmoothShade(true);
    mesh->setSurfaceColor({0.34, 0.36, 0.40});
    mesh->setEdgeWidth(0.0);
    drawn_ = mesh;
    polyscope::state::facePickIndStart = surface->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + surface->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + surface->nEdges();

    const double radius = lengthScale_ * 0.04;
    if (process == StepTree) {
        std::vector<Vector3> primalNodes;
        std::vector<std::array<size_t, 2>> primalEdges;
        for (Vertex vertex : surface->vertices()) {
            const auto parent = state_->tree.vertexParent.find(vertex);
            if (parent == state_->tree.vertexParent.end()) continue;
            primalNodes.push_back(geometry->inputVertexPositions[vertex]);
            primalEdges.push_back({vertex.getIndex(), parent->second.getIndex()});
        }
        polyscope::SurfaceGraphQuantity* tree = mesh->addSurfaceGraphQuantity("生成树", primalNodes, primalEdges);
        tree->setEnabled(true);
        tree->setColor({0.15, 0.35, 0.95});
        tree->setRadius(radius);

        std::vector<Vector3> dualNodes;
        std::vector<std::array<size_t, 2>> dualEdges;
        for (Face face : surface->faces()) {
            const auto parent = state_->tree.faceParent.find(face);
            if (parent == state_->tree.faceParent.end()) continue;
            dualNodes.push_back(faceCentroid(geometry, face));
            dualEdges.push_back({face.getIndex(), parent->second.getIndex()});
        }
        polyscope::SurfaceGraphQuantity* cotree = mesh->addSurfaceGraphQuantity("余树", dualNodes, dualEdges);
        cotree->setEnabled(true);
        cotree->setColor({0.1, 0.75, 0.65});
        cotree->setRadius(radius);
    } else if (process == StepGenerators) {
        std::vector<Vector3> nodes;
        std::vector<std::array<size_t, 2>> edges;
        size_t offset = 0;
        for (const std::vector<Halfedge>& loop : state_->tree.generators) {
            for (size_t i = 0; i < loop.size(); ++i) {
                const Halfedge he = loop[i];
                nodes.push_back(faceCentroid(geometry, he.face()));
                nodes.push_back((geometry->inputVertexPositions[he.tailVertex()] +
                                 geometry->inputVertexPositions[he.tipVertex()]) *
                                0.5);
                nodes.push_back(faceCentroid(geometry, he.twin().face()));
                edges.push_back({offset + 3 * i, offset + 3 * i + 1});
                edges.push_back({offset + 3 * i + 1, offset + 3 * i + 2});
            }
            offset += 3 * loop.size();
        }
        polyscope::SurfaceGraphQuantity* generators = mesh->addSurfaceGraphQuantity("基本环", nodes, edges);
        generators->setEnabled(true);
        generators->setColor({0.9, 0.25, 0.2});
        generators->setRadius(radius);
    } else if (process == StepOmega || process == StepExact || process == StepCoexact || process == StepHarmonic) {
        polyscope::SurfaceVectorQuantity* omega =
            mesh->addOneFormIntrinsicVectorQuantity("ω", state_->omega, polyscopeEdgeOrientations(*surface));
        polyscope::SurfaceVectorQuantity* exact =
            mesh->addOneFormIntrinsicVectorQuantity("dα", state_->exact, polyscopeEdgeOrientations(*surface));
        polyscope::SurfaceVectorQuantity* coexact =
            mesh->addOneFormIntrinsicVectorQuantity("δβ", state_->coexact, polyscopeEdgeOrientations(*surface));
        polyscope::SurfaceVectorQuantity* harmonic =
            mesh->addOneFormIntrinsicVectorQuantity("γ", state_->harmonicForm, polyscopeEdgeOrientations(*surface));
        enableOneForm(omega, process == StepOmega);
        enableOneForm(exact, process == StepExact);
        enableOneForm(coexact, process == StepCoexact);
        enableOneForm(harmonic, process == StepHarmonic);
    } else if (process == StepBasis) {
        for (size_t i = 0; i < state_->bases.size(); ++i) {
            polyscope::SurfaceVectorQuantity* basis = mesh->addOneFormIntrinsicVectorQuantity(
                "γ" + std::to_string(i), state_->bases[i], polyscopeEdgeOrientations(*surface));
            enableOneForm(basis, static_cast<int>(i) == basisIndex_);
        }
    } else if (process == StepConnection) {
        picking_ = true;
        connectionForm_ = mesh->addOneFormIntrinsicVectorQuantity("φ", state_->connection, polyscopeEdgeOrientations(*surface));
        enableOneForm(connectionForm_, showConnectionForm_);
        if (state_->connection.norm() > 0.0) drawDirection(mesh);
        showSelected();
    }

    if (resetCamera) polyscope::view::resetCameraToHomeView();
    status = note_;
    if (process == StepBasis) status += "    当前 γ" + std::to_string(basisIndex_);
    if (process == StepConnection) {
        status += "    指标和 " + std::to_string(singularitySum()) + " / 欧拉 " + std::to_string(euler_);
    }
}
