#pragma once

namespace topology {

namespace detail {

inline bool extractClosedMeshEdges(
    const Mesh& mesh,
    std::vector<TreeCotreeBasis::Edge>& edges,
    std::vector<int>& meshEdgeIndexForBasisEdge) {

    edges.clear();
    meshEdgeIndexForBasisEdge.clear();
    edges.reserve(mesh.edges.size());
    meshEdgeIndexForBasisEdge.reserve(mesh.edges.size());

    std::vector<int> interiorFaceIndex(mesh.faces.size(), -1);
    int nInteriorFaces = 0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (!face->isBoundary()) {
            interiorFaceIndex[face->index] = nInteriorFaces++;
        }
    }

    for (EdgeCIter edge = mesh.edges.begin(); edge != mesh.edges.end(); ++edge) {
        if (edge->isBoundary()) {
            return false;
        }

        TreeCotreeBasis::Edge basisEdge;
        basisEdge.v0 = edge->he->vertex->index;
        basisEdge.v1 = edge->he->flip->vertex->index;
        if (basisEdge.v0 > basisEdge.v1) {
            std::swap(basisEdge.v0, basisEdge.v1);
        }

        basisEdge.f0 = interiorFaceIndex[edge->he->face->index];
        basisEdge.f1 = interiorFaceIndex[edge->he->flip->face->index];
        if (basisEdge.f0 < 0 || basisEdge.f1 < 0) {
            return false;
        }

        edges.push_back(basisEdge);
        meshEdgeIndexForBasisEdge.push_back(edge->index);
    }

    return !edges.empty() && nInteriorFaces > 0;
}

} // namespace detail

inline std::vector<TreeCotreeBasis::Cycle> TreeCotreeBasis::buildClosedMeshBasis(
    const Mesh& mesh) {

    std::vector<Edge> edges;
    std::vector<int> meshEdgeIndexForBasisEdge;
    if (!detail::extractClosedMeshEdges(mesh, edges, meshEdgeIndexForBasisEdge)) {
        return {};
    }

    int nInteriorFaces = 0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (!face->isBoundary()) {
            ++nInteriorFaces;
        }
    }

    const int genus = detail::computeClosedMeshGenus(
        static_cast<int>(mesh.vertices.size()),
        static_cast<int>(mesh.edges.size()),
        nInteriorFaces);
    if (genus <= 0) {
        return {};
    }

    auto findEdge = [&](int a, int b) {
        if (a > b) {
            std::swap(a, b);
        }

        for (size_t i = 0; i < edges.size(); ++i) {
            if (edges[i].v0 == a && edges[i].v1 == b) {
                return static_cast<int>(i);
            }
        }

        return -1;
    };

    auto edgeSign = [&](int from, int to) {
        const int edgeIndex = findEdge(from, to);
        if (edgeIndex < 0) {
            return 0;
        }

        return (edges[edgeIndex].v0 == from && edges[edgeIndex].v1 == to) ? 1 : -1;
    };

    return buildClosedMeshBasis(
        static_cast<int>(mesh.vertices.size()),
        nInteriorFaces,
        edges,
        genus,
        findEdge,
        edgeSign,
        nullptr);
}

inline bool TreeCotreeBasis::decompose(const Mesh& mesh, Decomposition& out) {
    out = Decomposition();
    out.vertexCount = static_cast<int>(mesh.vertices.size());

    std::vector<Edge> edges;
    std::vector<int> meshEdgeIndexForBasisEdge;
    if (!detail::extractClosedMeshEdges(mesh, edges, meshEdgeIndexForBasisEdge)) {
        return false;
    }

    int nInteriorFaces = 0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (!face->isBoundary()) ++nInteriorFaces;
    }
    out.faceCount = nInteriorFaces;
    out.edgeCount = static_cast<int>(edges.size());

    const int genus = detail::computeClosedMeshGenus(
        out.vertexCount,
        static_cast<int>(mesh.edges.size()),
        nInteriorFaces);
    if (genus < 0) return false;
    out.genus = genus;

    auto findEdge = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        for (size_t i = 0; i < edges.size(); ++i) {
            if (edges[i].v0 == a && edges[i].v1 == b) return static_cast<int>(i);
        }
        return -1;
    };
    auto edgeSign = [&](int from, int to) {
        const int edgeIndex = findEdge(from, to);
        if (edgeIndex < 0) return 0;
        return (edges[edgeIndex].v0 == from && edges[edgeIndex].v1 == to) ? 1 : -1;
    };

    Parts parts;
    const std::vector<Cycle> cycles = buildClosedMeshBasis(
        out.vertexCount, nInteriorFaces, edges, genus, findEdge, edgeSign, &parts);

    auto toMeshEdge = [&](int basisEdge) {
        if (basisEdge < 0 || basisEdge >= static_cast<int>(meshEdgeIndexForBasisEdge.size())) return -1;
        return meshEdgeIndexForBasisEdge[static_cast<size_t>(basisEdge)];
    };
    auto appendMapped = [&](const std::vector<int>& source, std::vector<int>& dest) {
        dest.reserve(source.size());
        for (int basisEdge : source) {
            const int meshEdge = toMeshEdge(basisEdge);
            if (meshEdge >= 0) dest.push_back(meshEdge);
        }
    };
    appendMapped(parts.primalTree, out.primalTree);
    appendMapped(parts.dualTree, out.dualTree);
    appendMapped(parts.homologyEdges, out.homologyEdges);
    out.visitedVertices = parts.visitedVertices;
    out.visitedFaces = parts.visitedFaces;
    out.cycles.reserve(cycles.size());
    for (const Cycle& cycle : cycles) {
        std::vector<int> meshCycle;
        meshCycle.reserve(cycle.size());
        for (const SignedEdge& signedEdge : cycle) {
            const int meshEdge = toMeshEdge(signedEdge.first);
            if (meshEdge >= 0) meshCycle.push_back(meshEdge);
        }
        if (!meshCycle.empty()) out.cycles.push_back(std::move(meshCycle));
    }
    return true;
}

inline bool TreeCotreeBasis::buildCutSeamMesh(
    const Mesh& mesh,
    CutSeamMesh& cutMesh) {

    std::vector<Edge> edges;
    std::vector<int> meshEdgeIndexForBasisEdge;
    if (!detail::extractClosedMeshEdges(mesh, edges, meshEdgeIndexForBasisEdge)) {
        return false;
    }

    int nInteriorFaces = 0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (!face->isBoundary()) {
            ++nInteriorFaces;
        }
    }

    const int genus = detail::computeClosedMeshGenus(
        static_cast<int>(mesh.vertices.size()),
        static_cast<int>(mesh.edges.size()),
        nInteriorFaces);
    if (genus < 0) {
        return false;
    }

    if (genus == 0) {
        return cutMesh.buildFromSeamEdges(mesh, std::vector<int>());
    }

    auto findEdge = [&](int a, int b) {
        if (a > b) {
            std::swap(a, b);
        }

        for (size_t i = 0; i < edges.size(); ++i) {
            if (edges[i].v0 == a && edges[i].v1 == b) {
                return static_cast<int>(i);
            }
        }

        return -1;
    };

    auto edgeSign = [&](int from, int to) {
        const int edgeIndex = findEdge(from, to);
        if (edgeIndex < 0) {
            return 0;
        }

        return (edges[edgeIndex].v0 == from && edges[edgeIndex].v1 == to) ? 1 : -1;
    };

    std::vector<Cycle> cycles = buildClosedMeshBasis(
        static_cast<int>(mesh.vertices.size()),
        nInteriorFaces,
        edges,
        genus,
        findEdge,
        edgeSign);

    std::vector<int> seamEdgeIndices;
    seamEdgeIndices.reserve(cycles.size() * 8);
    for (const Cycle& cycle : cycles) {
        for (const SignedEdge& signedEdge : cycle) {
            const int basisEdgeIndex = signedEdge.first;
            if (basisEdgeIndex < 0 ||
                basisEdgeIndex >= static_cast<int>(meshEdgeIndexForBasisEdge.size())) {
                continue;
            }

            seamEdgeIndices.push_back(meshEdgeIndexForBasisEdge[basisEdgeIndex]);
        }
    }

    std::sort(seamEdgeIndices.begin(), seamEdgeIndices.end());
    seamEdgeIndices.erase(
        std::unique(seamEdgeIndices.begin(), seamEdgeIndices.end()),
        seamEdgeIndices.end());

    return cutMesh.buildFromSeamEdges(mesh, seamEdgeIndices);
}

inline CutSeamMesh TreeCotreeBasis::buildCutSeamMesh(const Mesh& mesh) {
    CutSeamMesh cutMesh;
    buildCutSeamMesh(mesh, cutMesh);
    return cutMesh;
}

} // namespace topology
