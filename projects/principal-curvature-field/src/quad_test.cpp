// Native quad-meshing test: switch models and compare
//   naive triangle pairing (own), MIQ (own), QuadCover (own),
//   and Geogram's periodic global parameterization (PGP).

#include "MIQQuad.h"
#include "Mesh.h"
#include "NaiveQuadMerge.h"
#include "QuadCover.h"
#include "geogram_quadcover.h"
#include "models.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/surface_parameterization_enums.h"
#include "polyscope/surface_parameterization_quantity.h"
#include "polyscope/view.h"

#include "imgui.h"

#include <array>
#include <chrono>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

void showSelected() {}

enum class Algo { Naive, Miq, QuadCoverAlgo, Geogram };

const char* kAlgoNames[] = {"朴素合并（自研）", "MIQ（自研）", "QuadCover（自研）", "Geogram PGP"};
constexpr int kAlgoCount = 4;

std::vector<ModelFile> gModels;
int gModel = 0;
Algo gAlgo = Algo::Naive;
std::string gStatus;
// Geogram PGP scaling: quad edge length in units of the average input edge.
float gQuadSize = 1.0f;
// MIQ solver knobs, exposed in the 算法参数 group.
int gMiqIterations = 12;
int gMiqJumpPasses = 2;
// 显示开关；每次 showCase 重建结构后按这些值恢复。
bool gShowCross = false;
bool gShowParam = false;
bool gShowEdges = false;
polyscope::SurfaceMesh* gPsMesh = nullptr;
polyscope::SurfaceParameterizationQuantity* gParamQuantity = nullptr;

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

// One grid unit == one average input edge, which is the convention Geogram uses
// for PGP scaling = 1. Applying it to every parameterization keeps the quad size
// of the three solvers comparable.
void normalizeGridUV(const std::vector<glm::vec3>& positions,
                     const std::vector<std::array<size_t, 3>>& faces,
                     std::vector<double>& cornerUV) {
    if (cornerUV.size() < faces.size() * 6) return;
    double worldSum = 0.0;
    double uvSum = 0.0;
    int worldCount = 0;
    int uvCount = 0;
    for (size_t face = 0; face < faces.size(); ++face) {
        for (int edge = 0; edge < 3; ++edge) {
            const size_t a = edge;
            const size_t b = (edge + 1) % 3;
            worldSum += glm::length(positions[faces[face][b]] - positions[faces[face][a]]);
            ++worldCount;
            const double ua = cornerUV[face * 6 + a * 2];
            const double va = cornerUV[face * 6 + a * 2 + 1];
            const double ub = cornerUV[face * 6 + b * 2];
            const double vb = cornerUV[face * 6 + b * 2 + 1];
            if (!std::isfinite(ua) || !std::isfinite(va) || !std::isfinite(ub) || !std::isfinite(vb)) continue;
            uvSum += std::hypot(ub - ua, vb - va);
            ++uvCount;
        }
    }
    if (worldCount == 0 || uvCount == 0) return;
    const double world = worldSum / worldCount;
    const double uv = uvSum / uvCount;
    if (!(world > 1e-12) || !(uv > 1e-12)) return;
    const double factor = world / uv;
    for (double& value : cornerUV) {
        if (std::isfinite(value)) value *= factor;
    }
}

// Per-vertex UV (rows = vertices) -> three (u, v) pairs per triangle.
std::vector<double> cornerUVFromVertexUV(const Eigen::MatrixXd& uv,
                                         const std::vector<std::array<size_t, 3>>& faces) {
    std::vector<double> cornerUV(faces.size() * 6, std::numeric_limits<double>::quiet_NaN());
    if (uv.cols() < 2) return cornerUV;
    for (size_t face = 0; face < faces.size(); ++face) {
        for (size_t corner = 0; corner < 3; ++corner) {
            const size_t vertex = faces[face][corner];
            if (vertex >= static_cast<size_t>(uv.rows())) continue;
            cornerUV[face * 6 + corner * 2] = uv(static_cast<Eigen::Index>(vertex), 0);
            cornerUV[face * 6 + corner * 2 + 1] = uv(static_cast<Eigen::Index>(vertex), 1);
        }
    }
    return cornerUV;
}

// Fraction of interior edges across which the integer isolines do not meet.
// A quad grid is continuous on an edge when the two charts differ only by a
// 90-degree rotation plus an integer translation; any fraction of a cell left
// over shows up as a gap in the drawn grid, so this is the number to watch.
double seamJumpFraction(const std::vector<std::array<size_t, 3>>& faces,
                        const std::vector<double>& cornerUV) {
    struct EdgeUV {
        double u[2] = {0.0, 0.0};
        double v[2] = {0.0, 0.0};
    };
    std::map<std::pair<size_t, size_t>, EdgeUV> seen;
    int total = 0;
    int broken = 0;
    for (size_t face = 0; face < faces.size(); ++face) {
        for (int edge = 0; edge < 3; ++edge) {
            const size_t a = static_cast<size_t>(edge);
            const size_t b = static_cast<size_t>((edge + 1) % 3);
            size_t lo = faces[face][a];
            size_t hi = faces[face][b];
            double ulo = cornerUV[face * 6 + a * 2];
            double vlo = cornerUV[face * 6 + a * 2 + 1];
            double uhi = cornerUV[face * 6 + b * 2];
            double vhi = cornerUV[face * 6 + b * 2 + 1];
            if (lo > hi) {
                std::swap(lo, hi);
                std::swap(ulo, uhi);
                std::swap(vlo, vhi);
            }
            const std::pair<size_t, size_t> key(lo, hi);
            const auto found = seen.find(key);
            if (found == seen.end()) {
                seen.emplace(key, EdgeUV{{ulo, uhi}, {vlo, vhi}});
                continue;
            }
            ++total;
            const EdgeUV& other = found->second;
            auto residual = [](double value) {
                return std::isfinite(value) ? std::abs(value - std::round(value)) : 1.0;
            };
            double best = 1.0;
            for (int k = 0; k < 4; ++k) {
                const double c = (k == 1 || k == 3) ? 0.0 : (k == 2 ? -1.0 : 1.0);
                const double s = (k == 0 || k == 2) ? 0.0 : (k == 1 ? 1.0 : -1.0);
                const double pu0 = c * other.u[0] - s * other.v[0];
                const double pv0 = s * other.u[0] + c * other.v[0];
                const double pu1 = c * other.u[1] - s * other.v[1];
                const double pv1 = s * other.u[1] + c * other.v[1];
                const double r = std::max(std::max(residual(ulo - pu0), residual(uhi - pu1)),
                                          std::max(residual(vlo - pv0), residual(vhi - pv1)));
                best = std::min(best, r);
            }
            if (!(best < 1.0e-4)) ++broken;
        }
    }
    return total > 0 ? static_cast<double>(broken) / static_cast<double>(total) : 0.0;
}

std::string seamText(const std::vector<std::array<size_t, 3>>& faces, const std::vector<double>& cornerUV) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "，接缝断边 %.2f%%", 100.0 * seamJumpFraction(faces, cornerUV));
    return buffer;
}

// Per-face corner UV packed as [u0 v0 u1 v1 u2 v2] rows -> flat array.
std::vector<double> cornerUVFromMatrix(const Eigen::MatrixXd& cornerUV) {
    std::vector<double> flat(static_cast<size_t>(cornerUV.rows()) * 6,
                             std::numeric_limits<double>::quiet_NaN());
    if (cornerUV.cols() < 6) return flat;
    for (Eigen::Index face = 0; face < cornerUV.rows(); ++face) {
        for (int k = 0; k < 6; ++k) {
            flat[static_cast<size_t>(face) * 6 + static_cast<size_t>(k)] = cornerUV(face, k);
        }
    }
    return flat;
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
        major->setEnabled(gShowCross);
    }
    if (!minorEdges.empty()) {
        polyscope::CurveNetwork* minor = polyscope::registerCurveNetwork("正交方向", minorNodes, minorEdges);
        minor->setColor({0.35f, 0.75f, 1.0f});
        minor->setRadius(0.0011f, true);
        minor->setEnabled(gShowCross);
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
    network->setColor({1.0f, 0.55f, 0.08f});
    network->setRadius(0.0019f, true);
    network->setEnabled(true);
}

// Isolines of integer u and integer v inside each triangle, i.e. the quad grid
// induced by a seamless parameterization. Returns the number of drawn segments.
int drawIntegerQuads(const std::vector<glm::vec3>& positions, const std::vector<std::array<size_t, 3>>& faces,
                     const std::vector<double>& cornerUV) {
    const int nFaces = static_cast<int>(faces.size());
    if (static_cast<int>(cornerUV.size()) < nFaces * 6) return 0;
    double edgeSum = 0.0;
    for (int face = 0; face < nFaces; ++face) {
        const glm::vec3& a = positions[faces[face][0]];
        const glm::vec3& b = positions[faces[face][1]];
        edgeSum += glm::length(b - a);
    }
    // 抬升必须大于线管半径（下面 setRadius 用的是相对包围盒的值），
    // 否则线管一半埋在曲面里，曲率一变就被面片吞掉，看起来像一格一格的断裂。
    glm::vec3 lo(1e30f), hi(-1e30f);
    for (const glm::vec3& p : positions) {
        lo = glm::min(lo, p);
        hi = glm::max(hi, p);
    }
    const float diagonal = glm::length(hi - lo);
    const float lift = std::max(0.0032f * diagonal, 1e-7f);
    const int cornerOf[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    for (int face = 0; face < nFaces; ++face) {
        glm::vec3 position[3];
        glm::vec3 normal(0.f);
        for (int corner = 0; corner < 3; ++corner) position[corner] = positions[faces[face][corner]];
        normal = glm::cross(position[1] - position[0], position[2] - position[0]);
        const float normalLength = glm::length(normal);
        if (normalLength > 1e-8f) normal /= normalLength;
        const glm::vec3 offset = normal * lift;
        static const int axisFilter = [] {
            const char* env = std::getenv("QUAD_TEST_AXIS");
            return env ? std::atoi(env) : -1;
        }();
        for (int axis = 0; axis < 2; ++axis) {
            if (axisFilter >= 0 && axis != axisFilter) continue;
            const double value[3] = {
                cornerUV[static_cast<size_t>(face) * 6 + axis],
                cornerUV[static_cast<size_t>(face) * 6 + 2 + axis],
                cornerUV[static_cast<size_t>(face) * 6 + 4 + axis],
            };
            if (!std::isfinite(value[0]) || !std::isfinite(value[1]) || !std::isfinite(value[2])) continue;
            const double lo = std::min(value[0], std::min(value[1], value[2]));
            const double hi = std::max(value[0], std::max(value[1], value[2]));
            const int first = static_cast<int>(std::floor(lo)) + 1;
            const int last = static_cast<int>(std::ceil(hi)) - 1;
            if (last < first || last - first > 24) continue;
            for (int level = first; level <= last; ++level) {
                glm::vec3 hit[2];
                int count = 0;
                for (int edge = 0; edge < 3 && count < 2; ++edge) {
                    const int a = cornerOf[edge][0];
                    const int b = cornerOf[edge][1];
                    const double va = value[a];
                    const double vb = value[b];
                    if ((va - level) * (vb - level) > 0.0) continue;
                    if (std::abs(vb - va) < 1e-14) continue;
                    const double t = (level - va) / (vb - va);
                    if (t < -1e-6 || t > 1.0 + 1e-6) continue;
                    hit[count++] = position[a] * static_cast<float>(1.0 - t) + position[b] * static_cast<float>(t);
                }
                if (count == 2 && glm::length(hit[1] - hit[0]) > 1e-8f) {
                    const size_t base = nodes.size();
                    nodes.push_back(hit[0] + offset);
                    nodes.push_back(hit[1] + offset);
                    edges.push_back({base, base + 1});
                }
            }
        }
    }
    if (edges.empty()) return 0;
    polyscope::CurveNetwork* network = polyscope::registerCurveNetwork("四边形", nodes, edges);
    network->setColor({1.0f, 0.55f, 0.08f});
    network->setRadius(0.0014f, true);
    network->setEnabled(true);
    return static_cast<int>(edges.size());
}

void drawParameterization(polyscope::SurfaceMesh* mesh, const Mesh& source) {
    std::vector<glm::vec2> coords(source.vertices.size(), glm::vec2(0.f));
    for (const Vertex& vertex : source.vertices) {
        coords[vertex.index] = glm::vec2(vertex.uv.x(), vertex.uv.y());
    }
    gParamQuantity = mesh->addVertexParameterizationQuantity("参数化", coords, polyscope::ParamCoordsType::WORLD);
    gParamQuantity->setEnabled(gShowParam);
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
    gParamQuantity = nullptr;
    gPsMesh = polyscope::registerSurfaceMesh(gModels[modelIndex].label, triangles.positions, triangles.faces);
    gPsMesh->setSmoothShade(true);
    gPsMesh->setSurfaceColor({0.34, 0.36, 0.40});
    gPsMesh->setEdgeColor({0.55, 0.58, 0.62});
    gPsMesh->setEdgeWidth(gShowEdges ? 1.0 : 0.0);
    if (std::getenv("QUAD_TEST_NO_SURFACE")) gPsMesh->setEnabled(false);  // 调试：只看线
    polyscope::SurfaceMesh* psMesh = gPsMesh;

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
        miq.setCrossFieldIterations(gMiqIterations);
        miq.setJumpRefinePasses(gMiqJumpPasses);
        miq.parameterize();
        miq.ensureTargetDirs();
        drawCross(triangles.positions, triangles.faces, miq.faceDir1(), miq.faceDir2());
        drawParameterization(psMesh, mesh);
        std::vector<double> cornerUV = cornerUVFromMatrix(miq.cornerUV());
        normalizeGridUV(triangles.positions, triangles.faces, cornerUV);
        const int segments = drawIntegerQuads(triangles.positions, triangles.faces, cornerUV);
        summary = "MIQ 参数化，" + std::to_string(segments) + " 条四边形边" +
                  seamText(triangles.faces, cornerUV);
    } else if (algo == Algo::QuadCoverAlgo) {
        QuadCover cover(mesh);
        cover.parameterize();
        const QuadCover::PrincipalField field = cover.principalField();
        const int nFaces = static_cast<int>(triangles.faces.size());
        Eigen::MatrixXd side = field.directions;
        for (int face = 0; face < nFaces && face < field.directions.rows(); ++face) {
            Eigen::Vector3d dir = field.directions.row(face);
            Eigen::Vector3d normal = field.normals.row(face);
            side.row(face) = normal.cross(dir);
        }
        drawCross(triangles.positions, triangles.faces, field.directions, side);
        drawParameterization(psMesh, mesh);
        std::vector<double> cornerUV = cornerUVFromVertexUV(cover.chartUV(), triangles.faces);
        normalizeGridUV(triangles.positions, triangles.faces, cornerUV);
        const int segments = drawIntegerQuads(triangles.positions, triangles.faces, cornerUV);
        summary = "QuadCover 参数化，" + std::to_string(segments) + " 条四边形边" +
                  seamText(triangles.faces, cornerUV);
    } else {
        const GeogramQuadCoverResult cover = runGeogramQuadCover(
            triangles.flatPositions.data(),
            static_cast<int>(triangles.positions.size()),
            triangles.flatFaces.data(),
            static_cast<int>(triangles.faces.size()),
            gQuadSize);
        if (!cover.ok) {
            summary = "Geogram PGP 失败  " + cover.error;
        } else {
            const int nFaces = static_cast<int>(triangles.faces.size());
            Eigen::MatrixXd direction(nFaces, 3);
            Eigen::MatrixXd side(nFaces, 3);
            for (int face = 0; face < nFaces; ++face) {
                for (int axis = 0; axis < 3; ++axis) {
                    direction(face, axis) = cover.direction[static_cast<size_t>(face) * 3 + axis];
                    side(face, axis) = cover.side[static_cast<size_t>(face) * 3 + axis];
                }
            }
            drawCross(triangles.positions, triangles.faces, direction, side);
            // Geogram already returns integer-grid units, and the grid only stays
            // continuous modulo integer translations rotated by multiples of 90
            // degrees. Rescaling it by a real factor would break the grid at
            // every seam, so the quad size is chosen through PGP's own scaling
            // argument instead.
            std::vector<double> cornerUV = cover.cornerUV;
            if (const char* dump = std::getenv("QUAD_TEST_DUMP")) {
                std::ofstream out(dump);
                for (size_t face = 0; face < triangles.faces.size(); ++face) {
                    out << triangles.faces[face][0] << ' ' << triangles.faces[face][1] << ' '
                        << triangles.faces[face][2];
                    for (int k = 0; k < 3; ++k) {
                        out << ' ' << cornerUV[face * 6 + k * 2] << ' ' << cornerUV[face * 6 + k * 2 + 1];
                    }
                    out << '\n';
                }
                std::cout << "已导出 UV: " << dump << std::endl;
            }
            const int segments = drawIntegerQuads(triangles.positions, triangles.faces, cornerUV);
            summary = "Geogram PGP，" + std::to_string(segments) + " 条四边形边" +
                      seamText(triangles.faces, cornerUV);
        }
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
    if (gModels.empty()) {
        ImGui::Text("没有找到 obj 模型");
        return;
    }
    // 主控件：模型、算法
    int nextModel = -1;
    if (ImGui::BeginCombo("模型", gModels[gModel].label.c_str())) {
        for (int i = 0; i < static_cast<int>(gModels.size()); ++i) {
            if (ImGui::Selectable(gModels[i].label.c_str(), i == gModel)) nextModel = i;
        }
        ImGui::EndCombo();
    }
    int nextAlgo = -1;
    if (ImGui::BeginCombo("算法", kAlgoNames[static_cast<int>(gAlgo)])) {
        for (int i = 0; i < kAlgoCount; ++i) {
            if (ImGui::Selectable(kAlgoNames[i], i == static_cast<int>(gAlgo))) nextAlgo = i;
        }
        ImGui::EndCombo();
    }
    if (!gStatus.empty()) ImGui::TextWrapped("%s", gStatus.c_str());

    bool solveAgain = false;
    if (ImGui::CollapsingHeader("算法参数")) {
        if (gAlgo == Algo::Miq) {
            if (ImGui::SliderInt("交叉场迭代", &gMiqIterations, 0, 30) && ImGui::IsItemDeactivatedAfterEdit()) {
                solveAgain = true;
            }
            if (ImGui::SliderInt("跳跃细化", &gMiqJumpPasses, 0, 6) && ImGui::IsItemDeactivatedAfterEdit()) {
                solveAgain = true;
            }
        } else if (gAlgo == Algo::Geogram) {
            if (ImGui::SliderFloat("四边形边长", &gQuadSize, 0.25f, 4.0f, "%.2f × 平均边长") &&
                ImGui::IsItemDeactivatedAfterEdit()) {
                solveAgain = true;
            }
            ImGui::Text("数值越大，四边形越少");
        } else {
            ImGui::Text("该算法没有可调参数");
        }
    }
    if (ImGui::CollapsingHeader("显示")) {
        if (ImGui::Checkbox("方向场", &gShowCross)) {
            if (polyscope::CurveNetwork* major = polyscope::getCurveNetwork("方向场")) major->setEnabled(gShowCross);
            if (polyscope::CurveNetwork* minor = polyscope::getCurveNetwork("正交方向")) minor->setEnabled(gShowCross);
        }
        if (ImGui::Checkbox("参数化纹理", &gShowParam) && gParamQuantity) {
            gParamQuantity->setEnabled(gShowParam);
        }
        if (ImGui::Checkbox("三角网格", &gShowEdges) && gPsMesh) {
            gPsMesh->setEdgeWidth(gShowEdges ? 1.0 : 0.0);
        }
    }

    if (nextModel >= 0 && nextModel != gModel) {
        showCase(nextModel, gAlgo);
    } else if (nextAlgo >= 0 && nextAlgo != static_cast<int>(gAlgo)) {
        showCase(gModel, static_cast<Algo>(nextAlgo));
    } else if (solveAgain) {
        showCase(gModel, gAlgo);
    }
}

// Optional arguments: [model.obj] [algo index 0..3], handy for headless checks.
int main(int argc, char** argv) {
    gModels = listModels();
    const int torus = modelIndexByName(gModels, "torus_F7680.obj");
    gModel = torus >= 0 ? torus : 0;
    if (gModels.empty()) {
        std::cerr << "没有找到 obj 模型" << std::endl;
        return 1;
    }
    if (argc > 1) {
        const int model = modelIndexByName(gModels, argv[1]);
        if (model >= 0) gModel = model;
    }
    if (argc > 2) {
        const int algo = std::atoi(argv[2]);
        if (algo >= 0 && algo < kAlgoCount) gAlgo = static_cast<Algo>(algo);
    }
    if (argc > 3) {
        const float size = static_cast<float>(std::atof(argv[3]));
        if (size >= 0.25f && size <= 4.f) gQuadSize = size;
    }
    polyscope::options::programName = "四边形网格";
    // 深色底 + 关闭地面，让四边形线成为画面主体。
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::bgColor = {0.09f, 0.10f, 0.12f, 1.f};
    polyscope::init();
    polyscope::state::userCallback = callback;
    if (!showCase(gModel, gAlgo)) return 1;
    if (const char* shot = std::getenv("QUAD_TEST_SHOT")) {
        const std::string path = shot;
        polyscope::state::userCallback = [path]() {
            static bool done = false;
            if (done) return;
            done = true;
            polyscope::screenshot(path, false);
            polyscope::shutdown();
        };
    }
    polyscope::show();
    return 0;
}
