// Native quad-meshing test: switch models and run naive merge, MIQ, or QuadCover.

#include "MIQQuad.h"
#include "Mesh.h"
#include "NaiveQuadMerge.h"
#include "QuadCover.h"
#include "models.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_parameterization_enums.h"
#include "polyscope/view.h"

#include "imgui.h"

#include <array>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

void showSelected() {}

enum class Algo { Naive, Miq, QuadCoverAlgo };

std::vector<ModelFile> gModels;
int gModel = 0;
Algo gAlgo = Algo::Naive;
std::string gStatus;

struct Triangles {
    std::vector<glm::vec3> positions;
    std::vector<std::array<size_t, 3>> faces;
    std::vector<double> flatPositions;
    std::vector<int> flatFaces;
};

Triangles collectTriangles(const Mesh& mesh) {
    Triangles out;
    out.positions.resize(mesh.vertices.size());
    out.flatPositions.resize(mesh.vertices.size() * 3);
    for (const Vertex& vertex : mesh.vertices) {
        out.positions[vertex.index] = glm::vec3(vertex.position.x(), vertex.position.y(), vertex.position.z());
        out.flatPositions[vertex.index * 3] = vertex.position.x();
        out.flatPositions[vertex.index * 3 + 1] = vertex.position.y();
        out.flatPositions[vertex.index * 3 + 2] = vertex.position.z();
    }
    for (const Face& face : mesh.faces) {
        if (face.isBoundary()) continue;
        const int a = face.he->vertex->index;
        const int b = face.he->next->vertex->index;
        const int c = face.he->next->next->vertex->index;
        out.faces.push_back({static_cast<size_t>(a), static_cast<size_t>(b), static_cast<size_t>(c)});
        out.flatFaces.push_back(a);
        out.flatFaces.push_back(b);
        out.flatFaces.push_back(c);
    }
    return out;
}

void drawCross(const std::vector<glm::vec3>& positions, const std::vector<std::array<size_t, 3>>& faces,
               const Eigen::MatrixXd& direction, const Eigen::MatrixXd& side) {
    const int nFaces = static_cast<int>(faces.size());
    if (direction.rows() < nFaces || side.rows() < nFaces) return;
    double edgeSum = 0;
    for (int face = 0; face < nFaces; ++face) {
        Eigen::Vector3d a(positions[faces[face][0]].x, positions[faces[face][0]].y, positions[faces[face][0]].z);
        Eigen::Vector3d b(positions[faces[face][1]].x, positions[faces[face][1]].y, positions[faces[face][1]].z);
        edgeSum += (b - a).norm();
    }
    const double half = nFaces > 0 ? 0.35 * edgeSum / nFaces : 0.02;
    std::vector<glm::vec3> majorNodes, minorNodes;
    std::vector<std::array<size_t, 2>> majorEdges, minorEdges;
    for (int face = 0; face < nFaces; ++face) {
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        for (int k = 0; k < 3; ++k) {
            const glm::vec3& p = positions[faces[face][k]];
            center += Eigen::Vector3d(p.x, p.y, p.z);
        }
        center /= 3.0;
        Eigen::Vector3d dir = direction.row(face);
        Eigen::Vector3d ortho = side.row(face);
        if (dir.norm() < 1e-8 || ortho.norm() < 1e-8) continue;
        dir.normalize();
        ortho.normalize();
        auto push = [&](std::vector<glm::vec3>& nodes, std::vector<std::array<size_t, 2>>& edges, const Eigen::Vector3d& axis) {
            const size_t base = nodes.size();
            Eigen::Vector3d a = center - axis * half;
            Eigen::Vector3d b = center + axis * half;
            nodes.push_back(glm::vec3(a.x(), a.y(), a.z()));
            nodes.push_back(glm::vec3(b.x(), b.y(), b.z()));
            edges.push_back({base, base + 1});
        };
        push(majorNodes, majorEdges, dir);
        push(minorNodes, minorEdges, ortho);
    }
    if (!majorEdges.empty()) {
        polyscope::CurveNetwork* major = polyscope::registerCurveNetwork("方向场", majorNodes, majorEdges);
        major->setColor({0.91f, 0.49f, 0.05f});
        major->setRadius(0.0014f, true);
    }
    if (!minorEdges.empty()) {
        polyscope::CurveNetwork* minor = polyscope::registerCurveNetwork("正交方向", minorNodes, minorEdges);
        minor->setColor({0.35f, 0.75f, 1.0f});
        minor->setRadius(0.0011f, true);
    }
}

void drawQuads(const std::vector<glm::vec3>& positions, const std::vector<int>& quads) {
    std::vector<std::array<size_t, 2>> edges;
    const int nQuads = static_cast<int>(quads.size() / 4);
    edges.reserve(nQuads * 4);
    for (int quad = 0; quad < nQuads; ++quad) {
        const int base = quad * 4;
        for (int k = 0; k < 4; ++k) {
            edges.push_back({
                static_cast<size_t>(quads[base + k]),
                static_cast<size_t>(quads[base + (k + 1) % 4]),
            });
        }
    }
    if (edges.empty()) return;
    polyscope::CurveNetwork* network = polyscope::registerCurveNetwork("四边形", positions, edges);
    network->setColor({0.91f, 0.49f, 0.05f});
    network->setRadius(0.0016f, true);
    network->setEnabled(true);
}

void drawParameterization(polyscope::SurfaceMesh* mesh, const Mesh& source) {
    std::vector<glm::vec2> coords(source.vertices.size(), glm::vec2(0.f));
    for (const Vertex& vertex : source.vertices) {
        coords[vertex.index] = glm::vec2(vertex.uv.x(), vertex.uv.y());
    }
    mesh->addVertexParameterizationQuantity("参数化", coords, polyscope::ParamCoordsType::WORLD)->setEnabled(true);
}

bool showCase(int modelIndex, Algo algo) {
    if (modelIndex < 0 || modelIndex >= static_cast<int>(gModels.size())) return false;
    const auto started = std::chrono::steady_clock::now();
    Mesh mesh;
    if (!mesh.read(gModels[modelIndex].path.string())) {
        gStatus = "无法读取 " + gModels[modelIndex].label;
        return false;
    }
    Triangles triangles = collectTriangles(mesh);
    polyscope::removeAllStructures();
    polyscope::SurfaceMesh* psMesh =
        polyscope::registerSurfaceMesh(gModels[modelIndex].label, triangles.positions, triangles.faces);
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({0.78, 0.78, 0.78});
    psMesh->setEdgeWidth(1.0);

    std::string summary;
    if (algo == Algo::Naive) {
        NaiveQuadMerge::Result result = NaiveQuadMerge::merge(
            triangles.flatPositions.data(),
            static_cast<int>(triangles.positions.size()),
            triangles.flatFaces.data(),
            static_cast<int>(triangles.faces.size()));
        drawQuads(triangles.positions, result.quads);
        summary = std::to_string(result.quads.size() / 4) + " 个四边形，剩余 " +
                  std::to_string(result.leftovers.size() / 3) + " 个三角形";
    } else if (algo == Algo::Miq) {
        MIQQuad miq(mesh);
        miq.setCrossFieldIterations(12);
        miq.setJumpRefinePasses(2);
        miq.parameterize();
        miq.ensureTargetDirs();
        drawCross(triangles.positions, triangles.faces, miq.faceDir1(), miq.faceDir2());
        drawParameterization(psMesh, mesh);
        summary = "MIQ 参数化";
    } else {
        QuadCover cover(mesh);
        cover.parameterize();
        QuadCover::PrincipalField field = cover.principalField();
        Eigen::MatrixXd side = field.directions;
        const int nFaces = static_cast<int>(triangles.faces.size());
        for (int face = 0; face < nFaces && face < field.directions.rows(); ++face) {
            Eigen::Vector3d dir = field.directions.row(face);
            Eigen::Vector3d normal = field.normals.row(face);
            side.row(face) = normal.cross(dir);
        }
        drawCross(triangles.positions, triangles.faces, field.directions, side);
        drawParameterization(psMesh, mesh);
        summary = "QuadCover 参数化";
    }

    polyscope::view::resetCameraToHomeView();
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    gModel = modelIndex;
    gAlgo = algo;
    gStatus = gModels[modelIndex].label + "    " + summary + "    " + std::to_string(seconds).substr(0, 4) + " 秒";
    std::cout << gStatus << std::endl;
    return true;
}

void callback() {
    ImGui::Text("四边形网格");
    if (gModels.empty()) {
        ImGui::Text("没有找到 obj 模型");
        return;
    }
    int nextModel = -1;
    if (ImGui::BeginCombo("模型", gModels[gModel].label.c_str())) {
        for (int i = 0; i < static_cast<int>(gModels.size()); ++i) {
            if (ImGui::Selectable(gModels[i].label.c_str(), i == gModel)) nextModel = i;
        }
        ImGui::EndCombo();
    }
    const char* algoNames[] = {"朴素合并", "MIQ", "QuadCover"};
    int nextAlgo = -1;
    if (ImGui::BeginCombo("算法", algoNames[static_cast<int>(gAlgo)])) {
        for (int i = 0; i < 3; ++i) {
            if (ImGui::Selectable(algoNames[i], i == static_cast<int>(gAlgo))) nextAlgo = i;
        }
        ImGui::EndCombo();
    }
    if (nextModel >= 0 && nextModel != gModel) showCase(nextModel, gAlgo);
    if (nextAlgo >= 0 && nextAlgo != static_cast<int>(gAlgo)) showCase(gModel, static_cast<Algo>(nextAlgo));
    if (!gStatus.empty()) ImGui::Text("%s", gStatus.c_str());
}

int main() {
    gModels = listModels();
    const int torus = modelIndexByName(gModels, "torus_F7680.obj");
    gModel = torus >= 0 ? torus : 0;
    if (gModels.empty()) {
        std::cerr << "没有找到 obj 模型" << std::endl;
        return 1;
    }
    polyscope::options::programName = "四边形网格";
    polyscope::init();
    polyscope::state::userCallback = callback;
    if (!showCase(gModel, gAlgo)) return 1;
    polyscope::show();
    return 0;
}
