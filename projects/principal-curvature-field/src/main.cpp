// QuadCover principal curvature field, before smoothing. Models can be switched live.

#include "Mesh.h"
#include "QuadCover.h"
#include "models.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"

#include "args/args.hxx"
#include "imgui.h"

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

void showSelected() {}

std::vector<ModelFile> gModels;
int gModel = 0;
std::string gStatus;

glm::vec3 at(const Eigen::MatrixXd& points, int row) {
    return glm::vec3(points(row, 0), points(row, 1), points(row, 2));
}

void drawField(polyscope::SurfaceMesh* mesh, const QuadCover::PrincipalField& field) {
    const int nFaces = static_cast<int>(field.directions.rows());
    double edgeSum = 0;
    for (int fi = 0; fi < nFaces; ++fi) {
        Eigen::Vector3d a = field.positions.row(field.triangles[fi * 3]);
        Eigen::Vector3d b = field.positions.row(field.triangles[fi * 3 + 1]);
        edgeSum += (b - a).norm();
    }
    const double half = nFaces > 0 ? 0.35 * edgeSum / nFaces : 0.02;

    std::vector<glm::vec3> majorNodes, minorNodes, arrows(nFaces);
    std::vector<std::array<size_t, 2>> majorEdges, minorEdges;
    for (int fi = 0; fi < nFaces; ++fi) {
        const int i0 = field.triangles[fi * 3];
        const int i1 = field.triangles[fi * 3 + 1];
        const int i2 = field.triangles[fi * 3 + 2];
        Eigen::Vector3d center = (field.positions.row(i0) + field.positions.row(i1) + field.positions.row(i2)) / 3.0;
        Eigen::Vector3d direction = field.directions.row(fi);
        Eigen::Vector3d normal = field.normals.row(fi);
        if (direction.norm() < 1e-8 || normal.norm() < 1e-8) continue;
        direction.normalize();
        normal.normalize();
        Eigen::Vector3d perpendicular = normal.cross(direction).normalized();
        auto push = [&](std::vector<glm::vec3>& nodes, std::vector<std::array<size_t, 2>>& edges, const Eigen::Vector3d& dir) {
            const size_t base = nodes.size();
            Eigen::Vector3d a = center - dir * half;
            Eigen::Vector3d b = center + dir * half;
            nodes.push_back(glm::vec3(a.x(), a.y(), a.z()));
            nodes.push_back(glm::vec3(b.x(), b.y(), b.z()));
            edges.push_back({base, base + 1});
        };
        push(majorNodes, majorEdges, direction);
        push(minorNodes, minorEdges, perpendicular);
        arrows[fi] = glm::vec3(direction.x(), direction.y(), direction.z());
    }

    if (!majorEdges.empty()) {
        polyscope::CurveNetwork* major = polyscope::registerCurveNetwork("主曲率方向", majorNodes, majorEdges);
        major->setColor({0.91f, 0.49f, 0.05f});
        major->setRadius(0.0014f, true);
        major->setEnabled(true);
    }
    if (!minorEdges.empty()) {
        polyscope::CurveNetwork* minor = polyscope::registerCurveNetwork("正交方向", minorNodes, minorEdges);
        minor->setColor({0.35f, 0.75f, 1.0f});
        minor->setRadius(0.0011f, true);
        minor->setEnabled(true);
    }
    mesh->addFaceVectorQuantity("主方向", arrows)->setEnabled(false);
}

bool showModel(int index) {
    if (index < 0 || index >= static_cast<int>(gModels.size())) return false;
    Mesh mesh;
    if (!mesh.read(gModels[index].path.string())) {
        gStatus = "无法读取 " + gModels[index].label;
        return false;
    }
    QuadCover cover(mesh);
    cover.pipelineReset();
    for (int stage = 0; stage <= 2; ++stage) {
        if (cover.pipelineRunStage(stage) != 0) {
            gStatus = std::string(cover.pipelineStageName(stage)) + " 失败";
            return false;
        }
    }
    QuadCover::PrincipalField field = cover.principalField();
    const int nFaces = static_cast<int>(field.directions.rows());
    std::vector<glm::vec3> positions(field.positions.rows());
    std::vector<std::array<size_t, 3>> triangles(nFaces);
    for (int i = 0; i < field.positions.rows(); ++i) positions[i] = at(field.positions, i);
    for (int face = 0; face < nFaces; ++face) {
        triangles[face] = {
            static_cast<size_t>(field.triangles[face * 3]),
            static_cast<size_t>(field.triangles[face * 3 + 1]),
            static_cast<size_t>(field.triangles[face * 3 + 2]),
        };
    }

    polyscope::removeAllStructures();
    polyscope::SurfaceMesh* psMesh = polyscope::registerSurfaceMesh(gModels[index].label, positions, triangles);
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({0.78, 0.78, 0.78});
    drawField(psMesh, field);
    polyscope::view::resetCameraToHomeView();
    gModel = index;
    gStatus = gModels[index].label + "    " + std::to_string(nFaces) + " 面";
    return true;
}

void callback() {
    ImGui::Text("QuadCover 主曲率场");
    ImGui::Text("光滑之前的中间结果");
    if (gModels.empty()) {
        ImGui::Text("没有找到 obj 模型");
        return;
    }
    const char* preview = gModels[gModel].label.c_str();
    int requested = -1;
    if (ImGui::BeginCombo("模型", preview)) {
        for (int i = 0; i < static_cast<int>(gModels.size()); ++i) {
            if (ImGui::Selectable(gModels[i].label.c_str(), i == gModel)) requested = i;
        }
        ImGui::EndCombo();
    }
    if (requested >= 0 && requested != gModel) showModel(requested);
    if (!gStatus.empty()) ImGui::Text("%s", gStatus.c_str());
    ImGui::Text("橙色是主方向，蓝色是切平面里的正交方向");
}

int main(int argc, char** argv) {
    args::ArgumentParser parser("QuadCover 主曲率场");
    args::Positional<std::string> inputFilename(parser, "mesh", "启动时打开的网格，默认识别 bunny.obj");
    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& error) {
        std::cerr << error.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    gModels = listModels();
    if (inputFilename) {
        ModelFile extra;
        extra.path = args::get(inputFilename);
        extra.label = extra.path.filename().string();
        const int existing = modelIndex(gModels, extra.path);
        if (existing >= 0) gModel = existing;
        else {
            gModels.insert(gModels.begin(), extra);
            gModel = 0;
        }
    } else {
        const int bunny = modelIndexByName(gModels, "bunny.obj");
        gModel = bunny >= 0 ? bunny : 0;
    }
    if (gModels.empty()) {
        std::cerr << "没有找到 obj 模型" << std::endl;
        return 1;
    }

    polyscope::options::programName = "QuadCover 主曲率场";
    polyscope::init();
    polyscope::state::userCallback = callback;
    if (!showModel(gModel)) return 1;
    polyscope::show();
    return 0;
}
