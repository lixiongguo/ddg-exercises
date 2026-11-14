#include "tree-cotree.h"

#include "TreeCotreeBasis.h"

#include <iostream>
#include <queue>
#include <utility>
#include <vector>

namespace {

void parentsFromLinks(int count, int root, const std::vector<std::pair<int, int>>& links,
                      std::vector<int>& parent) {
    std::vector<std::vector<int>> adjacent(static_cast<size_t>(count));
    for (const std::pair<int, int>& link : links) {
        if (link.first < 0 || link.second < 0 || link.first >= count || link.second >= count) continue;
        adjacent[static_cast<size_t>(link.first)].push_back(link.second);
        adjacent[static_cast<size_t>(link.second)].push_back(link.first);
    }

    parent.assign(static_cast<size_t>(count), -1);
    if (count <= 0 || root < 0 || root >= count) return;

    std::queue<int> pending;
    parent[static_cast<size_t>(root)] = root;
    pending.push(root);
    while (!pending.empty()) {
        const int current = pending.front();
        pending.pop();
        for (int next : adjacent[static_cast<size_t>(current)]) {
            if (parent[static_cast<size_t>(next)] >= 0) continue;
            parent[static_cast<size_t>(next)] = current;
            pending.push(next);
        }
    }
}

} // namespace

TreeCotree::TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {
    mesh = inputMesh;
    geometry = inputGeo;
}

void TreeCotree::build() {
    if (built_) return;
    built_ = true;
    vertexParent.clear();
    faceParent.clear();
    generators.clear();
    if (mesh == nullptr) return;

    const int nV = static_cast<int>(mesh->nVertices());
    const int nF = static_cast<int>(mesh->nFaces());
    const int nE = static_cast<int>(mesh->nEdges());
    if (nV <= 0 || nF <= 0 || nE <= 0) return;

    std::vector<topology::TreeCotreeBasis::Edge> edges(static_cast<size_t>(nE));
    for (Edge edge : mesh->edges()) {
        const int edgeIndex = static_cast<int>(edge.getIndex());
        if (edgeIndex < 0 || edgeIndex >= nE) continue;
        Halfedge he = edge.halfedge();
        topology::TreeCotreeBasis::Edge record;
        record.v0 = static_cast<int>(he.tailVertex().getIndex());
        record.v1 = static_cast<int>(he.tipVertex().getIndex());
        if (record.v0 > record.v1) std::swap(record.v0, record.v1);
        if (edge.isBoundary()) {
            record.f0 = -1;
            record.f1 = -1;
        } else {
            record.f0 = static_cast<int>(he.face().getIndex());
            record.f1 = static_cast<int>(he.twin().face().getIndex());
        }
        edges[static_cast<size_t>(edgeIndex)] = record;
    }

    const int chi = nV - nE + nF;
    const int slack = 2 - chi;
    const int genus = (slack >= 0 && slack % 2 == 0) ? slack / 2 : 0;

    auto findEdge = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        for (int edgeIndex = 0; edgeIndex < nE; ++edgeIndex) {
            const topology::TreeCotreeBasis::Edge& record = edges[static_cast<size_t>(edgeIndex)];
            if (record.v0 == a && record.v1 == b) return edgeIndex;
        }
        return -1;
    };
    auto edgeSign = [&](int from, int to) {
        const int edgeIndex = findEdge(from, to);
        if (edgeIndex < 0) return 0;
        const topology::TreeCotreeBasis::Edge& record = edges[static_cast<size_t>(edgeIndex)];
        return (record.v0 == from && record.v1 == to) ? 1 : -1;
    };

    topology::TreeCotreeBasis::Parts parts;
    const std::vector<topology::TreeCotreeBasis::Cycle> cycles =
        topology::TreeCotreeBasis::buildClosedMeshBasis(
            nV, nF, edges, genus, findEdge, edgeSign, &parts);

    std::vector<std::pair<int, int>> primalLinks;
    primalLinks.reserve(parts.primalTree.size());
    for (int edgeIndex : parts.primalTree) {
        if (edgeIndex < 0 || edgeIndex >= nE) continue;
        const topology::TreeCotreeBasis::Edge& record = edges[static_cast<size_t>(edgeIndex)];
        primalLinks.emplace_back(record.v0, record.v1);
    }
    std::vector<int> parentVertex;
    parentsFromLinks(nV, 0, primalLinks, parentVertex);
    for (Vertex vertex : mesh->vertices()) {
        const int index = static_cast<int>(vertex.getIndex());
        const int parent = (index >= 0 && index < nV && parentVertex[static_cast<size_t>(index)] >= 0)
                               ? parentVertex[static_cast<size_t>(index)]
                               : index;
        vertexParent[vertex] = mesh->vertex(static_cast<size_t>(parent));
    }

    std::vector<std::pair<int, int>> dualLinks;
    dualLinks.reserve(parts.dualTree.size());
    for (int edgeIndex : parts.dualTree) {
        if (edgeIndex < 0 || edgeIndex >= nE) continue;
        const topology::TreeCotreeBasis::Edge& record = edges[static_cast<size_t>(edgeIndex)];
        if (record.f0 < 0 || record.f1 < 0) continue;
        dualLinks.emplace_back(record.f0, record.f1);
    }
    std::vector<int> parentFace;
    parentsFromLinks(nF, 0, dualLinks, parentFace);
    for (Face face : mesh->faces()) {
        const int index = static_cast<int>(face.getIndex());
        const int parent = (index >= 0 && index < nF && parentFace[static_cast<size_t>(index)] >= 0)
                               ? parentFace[static_cast<size_t>(index)]
                               : index;
        faceParent[face] = mesh->face(static_cast<size_t>(parent));
    }

    generators.reserve(cycles.size());
    for (const topology::TreeCotreeBasis::Cycle& cycle : cycles) {
        std::vector<Halfedge> loop;
        loop.reserve(cycle.size());
        for (const topology::TreeCotreeBasis::SignedEdge& signedEdge : cycle) {
            const int edgeIndex = signedEdge.first;
            const int sign = signedEdge.second;
            if (edgeIndex < 0 || edgeIndex >= nE || sign == 0) continue;
            const topology::TreeCotreeBasis::Edge& record = edges[static_cast<size_t>(edgeIndex)];
            Edge edge = mesh->edge(static_cast<size_t>(edgeIndex));
            Halfedge he = edge.halfedge();
            const int from = sign > 0 ? record.v0 : record.v1;
            const int to = sign > 0 ? record.v1 : record.v0;
            const int tail = static_cast<int>(he.tailVertex().getIndex());
            const int tip = static_cast<int>(he.tipVertex().getIndex());
            if (!(tail == from && tip == to)) he = he.twin();
            loop.push_back(he);
        }
        if (!loop.empty()) generators.push_back(std::move(loop));
    }
}

void TreeCotree::buildPrimalSpanningTree() { build(); }

bool TreeCotree::inPrimalSpanningTree(Halfedge he) {
    build();
    const auto tail = vertexParent.find(he.tailVertex());
    const auto tip = vertexParent.find(he.tipVertex());
    if (tail == vertexParent.end() || tip == vertexParent.end()) return false;
    return tail->second == he.tipVertex() || tip->second == he.tailVertex();
}

void TreeCotree::buildDualSpanningCoTree() { build(); }

bool TreeCotree::inDualSpanningCotree(Halfedge he) {
    build();
    if (he.edge().isBoundary()) return false;
    const Face face = he.face();
    const Face twinFace = he.twin().face();
    const auto faceParentIt = faceParent.find(face);
    const auto twinParentIt = faceParent.find(twinFace);
    if (faceParentIt == faceParent.end() || twinParentIt == faceParent.end()) return false;
    return faceParentIt->second == twinFace || twinParentIt->second == face;
}

Halfedge TreeCotree::sharedHalfedge(Face f, Face g) const {
    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) return he;
    }
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

void TreeCotree::buildGenerators() { build(); }
