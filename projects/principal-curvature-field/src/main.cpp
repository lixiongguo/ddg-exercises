// QuadCover principal curvature field, before smoothing. Models can be switched live.

#include "Mesh.h"
#include "QuadCover.h"
#include "models.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"

#include "args/args.hxx"
#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

void showSelected() {}

std::vector<ModelFile> gModels;
int gModel = 0;
std::string gStatus;
std::string gCaption = "橙色流线顺着主曲率走，蓝色流线顺着正交方向走。红点是正奇异点，蓝点是负奇异点。";

glm::vec3 at(const Eigen::MatrixXd& points, int row) {
    return glm::vec3(points(row, 0), points(row, 1), points(row, 2));
}

glm::vec3 toGlm(const Eigen::Vector3d& point) {
    return glm::vec3(point.x(), point.y(), point.z());
}

constexpr double kPi = 3.14159265358979323846;

double foldQuarter(double delta) {
    const double period = kPi / 2.0;
    double wrapped = std::fmod(delta, period);
    if (wrapped > kPi / 4.0) wrapped -= period;
    else if (wrapped <= -kPi / 4.0) wrapped += period;
    return wrapped;
}

struct FieldPicture {
    std::vector<glm::vec3> majorNodes;
    std::vector<std::array<size_t, 2>> majorEdges;
    std::vector<glm::vec3> minorNodes;
    std::vector<std::array<size_t, 2>> minorEdges;
    std::vector<glm::vec3> arrows;
    std::vector<glm::vec3> singularPoints;
    std::vector<glm::vec3> singularColors;
    int majorCurves = 0;
    int minorCurves = 0;
    int positive = 0;
    int negative = 0;
};

struct Cell {
    int x = 0;
    int y = 0;
    int z = 0;
    bool operator==(const Cell& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct CellHash {
    size_t operator()(const Cell& cell) const {
        return static_cast<size_t>(cell.x) * 73856093u ^ static_cast<size_t>(cell.y) * 19349663u ^
               static_cast<size_t>(cell.z) * 83492791u;
    }
};

FieldPicture buildPicture(const QuadCover::PrincipalField& field) {
    FieldPicture picture;
    const int nFaces = static_cast<int>(field.directions.rows());
    const int nVerts = static_cast<int>(field.positions.rows());
    if (nFaces <= 0 || nVerts <= 0) return picture;

    std::vector<Eigen::Vector3d> position(nVerts);
    for (int vertex = 0; vertex < nVerts; ++vertex) position[vertex] = field.positions.row(vertex);

    std::vector<std::array<int, 3>> triangles(nFaces);
    std::vector<Eigen::Vector3d> direction(nFaces), normal(nFaces), minor(nFaces);
    std::vector<char> usable(nFaces, 0);
    double edgeSum = 0;
    int edgeCount = 0;
    std::unordered_map<uint64_t, std::array<int, 2>> neighbors;
    auto edgeKey = [](int a, int b) {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) | static_cast<uint32_t>(b);
    };
    std::vector<std::vector<int>> incident(nVerts);
    for (int face = 0; face < nFaces; ++face) {
        const int i0 = field.triangles[face * 3];
        const int i1 = field.triangles[face * 3 + 1];
        const int i2 = field.triangles[face * 3 + 2];
        triangles[face] = {i0, i1, i2};
        incident[i0].push_back(face);
        incident[i1].push_back(face);
        incident[i2].push_back(face);
        Eigen::Vector3d dir = field.directions.row(face);
        Eigen::Vector3d nrm = field.normals.row(face);
        if (dir.norm() > 1e-8 && nrm.norm() > 1e-8) {
            dir.normalize();
            nrm.normalize();
            dir -= nrm * dir.dot(nrm);
            if (dir.norm() > 1e-8) {
                dir.normalize();
                direction[face] = dir;
                normal[face] = nrm;
                minor[face] = nrm.cross(dir).normalized();
                usable[face] = 1;
                picture.arrows.push_back(toGlm(dir));
            } else {
                picture.arrows.push_back(glm::vec3(0.f));
            }
        } else {
            picture.arrows.push_back(glm::vec3(0.f));
        }
        const int corner[3] = {i0, i1, i2};
        for (int edge = 0; edge < 3; ++edge) {
            const int a = corner[edge];
            const int b = corner[(edge + 1) % 3];
            edgeSum += (position[a] - position[b]).norm();
            ++edgeCount;
            auto& slot = neighbors[edgeKey(a, b)];
            if (slot[0] == 0 && slot[1] == 0) slot = {-1, -1};
            if (slot[0] < 0) slot[0] = face;
            else if (slot[1] < 0 && slot[0] != face) slot[1] = face;
        }
    }
    const double meanEdge = edgeCount > 0 ? edgeSum / edgeCount : 0.02;
    const double spacing = std::max(meanEdge * 5.0, 1e-6);

    auto otherFace = [&](int face, int a, int b) {
        const auto found = neighbors.find(edgeKey(a, b));
        if (found == neighbors.end()) return -1;
        if (found->second[0] == face) return found->second[1];
        if (found->second[1] == face) return found->second[0];
        return -1;
    };

    auto cellOf = [&](const Eigen::Vector3d& point) {
        Cell cell;
        cell.x = static_cast<int>(std::floor(point.x() / spacing));
        cell.y = static_cast<int>(std::floor(point.y() / spacing));
        cell.z = static_cast<int>(std::floor(point.z() / spacing));
        return cell;
    };

    auto trace = [&](int startFace, const Eigen::Vector3d& start, Eigen::Vector3d dir, bool useMinor,
                     const std::unordered_map<Cell, int, CellHash>& taken) {
        std::vector<Eigen::Vector3d> points;
        std::vector<char> seen(nFaces, 0);
        if (!usable[startFace] || dir.norm() < 1e-8) return points;
        dir.normalize();
        int face = startFace;
        Eigen::Vector3d point = start;
        seen[startFace] = 1;
        for (int step = 0; step < 160; ++step) {
            const auto& tri = triangles[face];
            const Eigen::Vector3d& p0 = position[tri[0]];
            const Eigen::Vector3d& p1 = position[tri[1]];
            const Eigen::Vector3d& p2 = position[tri[2]];
            Eigen::Vector3d e1 = p1 - p0;
            Eigen::Vector3d e2 = p2 - p0;
            const double d00 = e1.dot(e1);
            const double d01 = e1.dot(e2);
            const double d11 = e2.dot(e2);
            const double denom = d00 * d11 - d01 * d01;
            if (std::abs(denom) < 1e-18) break;
            auto coords = [&](const Eigen::Vector3d& value, double& b0, double& b1, double& b2) {
                const double d20 = value.dot(e1);
                const double d21 = value.dot(e2);
                b1 = (d11 * d20 - d01 * d21) / denom;
                b2 = (d00 * d21 - d01 * d20) / denom;
                b0 = 1.0 - b1 - b2;
            };
            double b0 = 0, b1 = 0, b2 = 0, t0 = 0, t1 = 0, t2 = 0;
            coords(point - p0, b0, b1, b2);
            coords(dir, t0, t1, t2);
            const double db0 = -t1 - t2;
            const double bary[3] = {b0, b1, b2};
            const double velocity[3] = {db0, t1, t2};
            int hit = -1;
            double bestT = 1e100;
            for (int corner = 0; corner < 3; ++corner) {
                if (velocity[corner] < -1e-9) {
                    const double time = -bary[corner] / velocity[corner];
                    if (time > 1e-5 && time < bestT) {
                        bestT = time;
                        hit = corner;
                    }
                }
            }
            if (hit < 0 || bestT > meanEdge * 8.0) break;
            Eigen::Vector3d edgePoint = point + dir * bestT;
            edgePoint -= normal[face] * (edgePoint - p0).dot(normal[face]);
            const int va = tri[(hit + 1) % 3];
            const int vb = tri[(hit + 2) % 3];
            const int next = otherFace(face, va, vb);
            if ((edgePoint - point).norm() > meanEdge * 1e-4) points.push_back(edgePoint);
            if (next < 0 || !usable[next]) break;
            Eigen::Vector3d incoming = dir - normal[next] * dir.dot(normal[next]);
            if (incoming.norm() < 1e-8) break;
            incoming.normalize();
            Eigen::Vector3d fieldDir = useMinor ? minor[next] : direction[next];
            if (std::abs(fieldDir.dot(incoming)) < 0.5) break;
            if (fieldDir.dot(incoming) < 0) fieldDir = -fieldDir;
            const Eigen::Vector3d center = (position[triangles[next][0]] + position[triangles[next][1]] + position[triangles[next][2]]) / 3.0;
            const Eigen::Vector3d inside = edgePoint * 0.98 + center * 0.02;
            if (taken.find(cellOf(inside)) != taken.end()) break;
            if (seen[next]) break;
            seen[next] = 1;
            point = inside;
            dir = fieldDir;
            face = next;
        }
        return points;
    };

    auto addStream = [&](bool useMinor, std::vector<glm::vec3>& nodes, std::vector<std::array<size_t, 2>>& edges, int& curves) {
        std::unordered_map<Cell, int, CellHash> taken;
        for (int face = 0; face < nFaces; ++face) {
            if (!usable[face]) continue;
            const Eigen::Vector3d center =
                (position[triangles[face][0]] + position[triangles[face][1]] + position[triangles[face][2]]) / 3.0;
            if (taken.find(cellOf(center)) != taken.end()) continue;
            Eigen::Vector3d dir = useMinor ? minor[face] : direction[face];
            std::vector<Eigen::Vector3d> backward = trace(face, center, -dir, useMinor, taken);
            std::vector<Eigen::Vector3d> forward = trace(face, center, dir, useMinor, taken);
            std::vector<Eigen::Vector3d> curve;
            curve.reserve(backward.size() + forward.size() + 1);
            for (auto it = backward.rbegin(); it != backward.rend(); ++it) curve.push_back(*it);
            curve.push_back(center);
            curve.insert(curve.end(), forward.begin(), forward.end());
            if (curve.size() < 4) continue;
            double length = 0;
            for (size_t i = 1; i < curve.size(); ++i) length += (curve[i] - curve[i - 1]).norm();
            if (length < meanEdge * 3.0) continue;
            const size_t base = nodes.size();
            for (const Eigen::Vector3d& point : curve) {
                nodes.push_back(toGlm(point));
                taken[cellOf(point)] = 1;
            }
            for (size_t i = 1; i < curve.size(); ++i) edges.push_back({base + i - 1, base + i});
            ++curves;
        }
    };
    addStream(false, picture.majorNodes, picture.majorEdges, picture.majorCurves);
    addStream(true, picture.minorNodes, picture.minorEdges, picture.minorCurves);

    std::vector<double> indexAt(nVerts, 0);
    std::vector<Eigen::Vector3d> normalAt(nVerts, Eigen::Vector3d::Zero());
    std::vector<char> candidate(nVerts, 0);
    for (int vertex = 0; vertex < nVerts; ++vertex) {
        if (incident[vertex].size() < 3) continue;
        std::vector<int> ring;
        int face = incident[vertex][0];
        bool closed = false;
        for (size_t guard = 0; guard < incident[vertex].size() + 1; ++guard) {
            const auto& tri = triangles[face];
            int corner = -1;
            for (int k = 0; k < 3; ++k) {
                if (tri[k] == vertex) corner = k;
            }
            if (corner < 0) break;
            ring.push_back(face);
            const int next = otherFace(face, vertex, tri[(corner + 1) % 3]);
            if (next < 0) break;
            face = next;
            if (face == incident[vertex][0]) {
                closed = ring.size() >= 3;
                break;
            }
        }
        if (!closed) continue;
        double angleSum = 0;
        Eigen::Vector3d vertexNormal = Eigen::Vector3d::Zero();
        bool defined = true;
        for (int faceId : ring) {
            if (!usable[faceId]) {
                defined = false;
                break;
            }
            const auto& tri = triangles[faceId];
            int corner = 0;
            for (int k = 0; k < 3; ++k) {
                if (tri[k] == vertex) corner = k;
            }
            Eigen::Vector3d left = position[tri[(corner + 1) % 3]] - position[vertex];
            Eigen::Vector3d right = position[tri[(corner + 2) % 3]] - position[vertex];
            const double lengths = left.norm() * right.norm();
            if (lengths < 1e-14) {
                defined = false;
                break;
            }
            const double cosine = std::max(-1.0, std::min(1.0, left.dot(right) / lengths));
            angleSum += std::acos(cosine);
            vertexNormal += normal[faceId] * left.cross(right).norm();
        }
        if (!defined || vertexNormal.norm() < 1e-10) continue;
        vertexNormal.normalize();
        Eigen::Vector3d axis = std::abs(vertexNormal.x()) < 0.9 ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
        Eigen::Vector3d tangent = vertexNormal.cross(axis).normalized();
        Eigen::Vector3d bitangent = vertexNormal.cross(tangent);
        double turned = 0;
        double turnedAbs = 0;
        for (size_t i = 0; i < ring.size(); ++i) {
            auto angleOf = [&](int faceId) {
                Eigen::Vector3d projected = direction[faceId] - vertexNormal * direction[faceId].dot(vertexNormal);
                if (projected.norm() < 1e-8) projected = tangent;
                return std::atan2(projected.dot(bitangent), projected.dot(tangent));
            };
            const double step = foldQuarter(angleOf(ring[(i + 1) % ring.size()]) - angleOf(ring[i]));
            turned += step;
            turnedAbs += std::abs(step);
        }
        indexAt[vertex] = (turned - (2.0 * kPi - angleSum)) / (2.0 * kPi);
        normalAt[vertex] = vertexNormal;
        const double quarter = std::round(indexAt[vertex] * 4.0) / 4.0;
        const bool coherent = turnedAbs <= std::abs(turned) + 0.35;
        if (coherent && std::abs(quarter) >= 0.2 && std::abs(indexAt[vertex] - quarter) <= 0.08) candidate[vertex] = 1;
    }
    std::vector<std::vector<int>> vertexNeighbors(nVerts);
    for (int face = 0; face < nFaces; ++face) {
        const auto& tri = triangles[face];
        for (int corner = 0; corner < 3; ++corner) {
            vertexNeighbors[tri[corner]].push_back(tri[(corner + 1) % 3]);
        }
    }
    std::unordered_map<Cell, int, CellHash> strongest;
    const double singularSpacing = std::max(meanEdge * 12.0, 1e-6);
    auto singularCell = [&](const Eigen::Vector3d& point) {
        Cell cell;
        cell.x = static_cast<int>(std::floor(point.x() / singularSpacing));
        cell.y = static_cast<int>(std::floor(point.y() / singularSpacing));
        cell.z = static_cast<int>(std::floor(point.z() / singularSpacing));
        return cell;
    };
    for (int vertex = 0; vertex < nVerts; ++vertex) {
        if (!candidate[vertex]) continue;
        bool peak = true;
        for (int neighbor : vertexNeighbors[vertex]) {
            if (std::abs(indexAt[neighbor]) > std::abs(indexAt[vertex]) + 1e-6) peak = false;
        }
        if (!peak) continue;
        const Cell cell = singularCell(position[vertex]);
        const auto found = strongest.find(cell);
        if (found == strongest.end() || std::abs(indexAt[vertex]) > std::abs(indexAt[found->second])) {
            strongest[cell] = vertex;
        }
    }
    for (const auto& item : strongest) {
        const int vertex = item.second;
        const double quarter = std::round(indexAt[vertex] * 4.0) / 4.0;
        Eigen::Vector3d point = position[vertex] + normalAt[vertex] * meanEdge * 0.35;
        picture.singularPoints.push_back(toGlm(point));
        if (quarter > 0) {
            picture.singularColors.push_back(glm::vec3(0.95f, 0.22f, 0.12f));
            ++picture.positive;
        } else {
            picture.singularColors.push_back(glm::vec3(0.16f, 0.45f, 0.98f));
            ++picture.negative;
        }
    }
    return picture;
}

void drawField(polyscope::SurfaceMesh* mesh, const FieldPicture& picture) {
    if (picture.arrows.size() == mesh->nFaces()) {
        mesh->addFaceVectorQuantity("主方向", picture.arrows)->setEnabled(false);
    }
    if (!picture.majorEdges.empty()) {
        polyscope::CurveNetwork* major = polyscope::registerCurveNetwork("主曲率流线", picture.majorNodes, picture.majorEdges);
        major->setColor({0.91f, 0.49f, 0.05f});
        major->setRadius(0.0022f, true);
        major->setEnabled(true);
    }
    if (!picture.minorEdges.empty()) {
        polyscope::CurveNetwork* minor = polyscope::registerCurveNetwork("正交流线", picture.minorNodes, picture.minorEdges);
        minor->setColor({0.35f, 0.75f, 1.0f});
        minor->setRadius(0.0018f, true);
        minor->setEnabled(true);
    }
    if (!picture.singularPoints.empty()) {
        polyscope::PointCloud* points = polyscope::registerPointCloud("奇异点", picture.singularPoints);
        points->setPointRadius(0.0065f, true);
        points->addColorQuantity("指标", picture.singularColors)->setEnabled(true);
    }
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

    const FieldPicture picture = buildPicture(field);
    polyscope::removeAllStructures();
    polyscope::SurfaceMesh* psMesh = polyscope::registerSurfaceMesh(gModels[index].label, positions, triangles);
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({0.78, 0.78, 0.78});
    drawField(psMesh, picture);
    polyscope::view::resetCameraToHomeView();
    gModel = index;
    gStatus = gModels[index].label + "    " + std::to_string(nFaces) + " 面    流线 " +
              std::to_string(picture.majorCurves) + " / " + std::to_string(picture.minorCurves) + "    奇异点 +" +
              std::to_string(picture.positive) + " −" + std::to_string(picture.negative);
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
    if (!gStatus.empty()) ImGui::TextWrapped("%s", gStatus.c_str());
    ImGui::TextWrapped("%s", gCaption.c_str());
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
