// Principal curvature field from geometry-central's dihedral-angle curvature tensor.

#include "QuadCover.h"
#include "models.h"

#include "ComplexPolyField.h"
#include "MIQQuad.h"
#include "Mesh.h"
#include "NRosyVectorFields.h"

#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"

#include "args/args.hxx"
#include "imgui.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

void showSelected() {}

struct FieldAlgo {
    const char* name;
    const char* caption;
};

const FieldAlgo kAlgos[] = {
    {"曲率对齐", "曲率对齐的光滑十字场。橙色和蓝色是两族互相垂直的流线，圆点是奇异点。"},
    {"主曲率", "每个三角面的主曲率方向，没有再做平滑。"},
    {"平凡联络", "平凡联络生成的十字场。奇异点放在曲率绝对值最大的顶点上，个数由拓扑决定。"},
    {"复多项式", "用复多项式表示的光滑十字场，奇异点由平滑能量决定。"},
};

std::vector<ModelFile> gModels;
int gModel = 0;
int gAlgo = 0;
int gQuadMethod = 0;
float gDensity = 20.f;
bool gShowEdges = true;
bool gReady = false;
polyscope::SurfaceMesh* gPsMesh = nullptr;
Mesh gSurface;
std::string gStatus;
std::string gCaption = kAlgos[0].caption;

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

QuadCover::PrincipalField gField;
FieldPicture gSingularities;

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

FieldPicture buildPicture(const QuadCover::PrincipalField& field, float spacingEdges) {
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
    const double spacing = std::max(meanEdge * static_cast<double>(spacingEdges), 1e-6);

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
        for (int step = 0; step < 420; ++step) {
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
            if (std::abs(fieldDir.dot(incoming)) < 0.08) break;
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
            if (curve.size() < 3) continue;
            double length = 0;
            for (size_t i = 1; i < curve.size(); ++i) length += (curve[i] - curve[i - 1]).norm();
            if (length < meanEdge) continue;
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
        major->setColor({0.72f, 0.34f, 0.05f});
        major->setMaterial("flat");
        major->setRadius(0.00018f, true);
        major->setEnabled(true);
    }
    if (!picture.minorEdges.empty()) {
        polyscope::CurveNetwork* minor = polyscope::registerCurveNetwork("正交流线", picture.minorNodes, picture.minorEdges);
        minor->setColor({0.12f, 0.42f, 0.68f});
        minor->setMaterial("flat");
        minor->setRadius(0.00015f, true);
        minor->setEnabled(true);
    }
    if (!picture.singularPoints.empty()) {
        polyscope::PointCloud* points = polyscope::registerPointCloud("奇异点", picture.singularPoints);
        points->setMaterial("flat");
        points->setPointRadius(0.0026f, true);
        points->addColorQuantity("指标", picture.singularColors)->setEnabled(true);
    }
}

double meanEdgeLength(const Mesh& mesh) {
    double sum = 0;
    int count = 0;
    for (EdgeCIter edge = mesh.edges.begin(); edge != mesh.edges.end(); ++edge) {
        sum += edge->length();
        ++count;
    }
    return count > 0 ? sum / count : 0.01;
}

void paintSingularities(const std::vector<Eigen::Vector3d>& points, const std::vector<Eigen::Vector3d>& normals,
                        const std::vector<char>& boundary, const std::vector<int>& index, double lift,
                        FieldPicture& singularities) {
    singularities = {};
    const int nVerts = static_cast<int>(points.size());
    for (int vertex = 0; vertex < nVerts; ++vertex) {
        if (boundary[vertex] || index[vertex] == 0) continue;
        const Eigen::Vector3d point = points[vertex] + normals[vertex] * lift;
        singularities.singularPoints.push_back(toGlm(point));
        if (index[vertex] > 0) {
            singularities.singularColors.push_back(glm::vec3(0.95f, 0.22f, 0.12f));
            ++singularities.positive;
        } else {
            singularities.singularColors.push_back(glm::vec3(0.16f, 0.45f, 0.98f));
            ++singularities.negative;
        }
    }
}

void paintMeshSingularities(const Mesh& mesh, const std::vector<int>& index, FieldPicture& singularities) {
    const int nVerts = static_cast<int>(mesh.vertices.size());
    std::vector<Eigen::Vector3d> points(nVerts), normals(nVerts);
    std::vector<char> boundary(nVerts, 0);
    for (VertexCIter vertex = mesh.vertices.begin(); vertex != mesh.vertices.end(); ++vertex) {
        points[vertex->index] = vertex->position;
        normals[vertex->index] = vertex->normal();
        boundary[vertex->index] = vertex->isBoundary() || vertex->isIsolated();
    }
    paintSingularities(points, normals, boundary, index, meanEdgeLength(mesh) * 0.12, singularities);
}

bool exportInterior(const Mesh& mesh, QuadCover::PrincipalField& field, std::string& error) {
    const int nVerts = static_cast<int>(mesh.vertices.size());
    int nFaces = 0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (!face->isBoundary()) ++nFaces;
    }
    if (nVerts <= 0 || nFaces <= 0) {
        error = "网格没有面";
        return false;
    }
    field.positions.resize(nVerts, 3);
    field.triangles.resize(nFaces * 3);
    field.directions.resize(nFaces, 3);
    field.normals.resize(nFaces, 3);
    field.directions.setZero();
    for (VertexCIter vertex = mesh.vertices.begin(); vertex != mesh.vertices.end(); ++vertex) {
        field.positions.row(vertex->index) = vertex->position;
    }
    int faceIndex = 0;
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        field.triangles[faceIndex * 3] = face->he->vertex->index;
        field.triangles[faceIndex * 3 + 1] = face->he->next->vertex->index;
        field.triangles[faceIndex * 3 + 2] = face->he->next->next->vertex->index;
        Eigen::Vector3d normal = face->normal();
        const double length = normal.norm();
        if (length > 1e-12) normal /= length;
        field.normals.row(faceIndex) = normal;
        ++faceIndex;
    }
    return true;
}

bool singularitiesFromDirections(const QuadCover::PrincipalField& field, FieldPicture& singularities, std::string& error) {
    namespace gc = geometrycentral;
    namespace surface = geometrycentral::surface;
    const int nFaces = static_cast<int>(field.directions.rows());
    Eigen::MatrixXi faces(nFaces, 3);
    for (int face = 0; face < nFaces; ++face) {
        faces.row(face) << field.triangles[face * 3], field.triangles[face * 3 + 1], field.triangles[face * 3 + 2];
    }
    std::unique_ptr<surface::ManifoldSurfaceMesh> mesh;
    std::unique_ptr<surface::VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = surface::makeManifoldSurfaceMeshAndGeometry(field.positions, faces);
    if (static_cast<int>(mesh->nVertices()) != field.positions.rows() || static_cast<int>(mesh->nFaces()) != nFaces) {
        error = "网格顶点和面无法对应";
        return false;
    }
    geometry->requireFaceTangentBasis();
    geometry->requireVertexNormals();
    surface::FaceData<gc::Vector2> powered(*mesh);
    for (surface::Face face : mesh->faces()) {
        const int index = face.getIndex();
        const gc::Vector3 direction{field.directions(index, 0), field.directions(index, 1), field.directions(index, 2)};
        const gc::Vector3 basisX = geometry->faceTangentBasis[face][0];
        const gc::Vector3 basisY = geometry->faceTangentBasis[face][1];
        const double theta = std::atan2(geometrycentral::dot(direction, basisY), geometrycentral::dot(direction, basisX));
        powered[face] = gc::Vector2::fromAngle(4.0 * theta);
    }
    const surface::VertexData<int> index = surface::computeVertexIndex(*geometry, powered, 4);
    std::vector<Eigen::Vector3d> points(mesh->nVertices()), normals(mesh->nVertices());
    std::vector<char> boundary(mesh->nVertices(), 0);
    std::vector<int> values(mesh->nVertices(), 0);
    for (surface::Vertex vertex : mesh->vertices()) {
        const size_t id = vertex.getIndex();
        const gc::Vector3 position = geometry->inputVertexPositions[vertex];
        const gc::Vector3 normal = geometry->vertexNormals[vertex];
        points[id] = Eigen::Vector3d(position.x, position.y, position.z);
        normals[id] = Eigen::Vector3d(normal.x, normal.y, normal.z);
        boundary[id] = vertex.isBoundary();
        values[id] = index[vertex];
    }
    double edgeSum = 0;
    int edgeCount = 0;
    for (surface::Edge edge : mesh->edges()) {
        edgeSum += geometry->edgeLength(edge);
        ++edgeCount;
    }
    paintSingularities(points, normals, boundary, values, (edgeCount > 0 ? edgeSum / edgeCount : 0.01) * 0.12, singularities);
    return true;
}

bool computeAligned(QuadCover::PrincipalField& field, FieldPicture& singularities, std::string& error) {
    namespace gc = geometrycentral;
    namespace surface = geometrycentral::surface;
    const int nFaces = static_cast<int>(field.directions.rows());
    Eigen::MatrixXi faces(nFaces, 3);
    for (int face = 0; face < nFaces; ++face) {
        faces.row(face) << field.triangles[face * 3], field.triangles[face * 3 + 1], field.triangles[face * 3 + 2];
    }
    std::unique_ptr<surface::ManifoldSurfaceMesh> mesh;
    std::unique_ptr<surface::VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = surface::makeManifoldSurfaceMeshAndGeometry(field.positions, faces);
    geometry->requireFaceTangentBasis();
    geometry->requireFaceNormals();
    const surface::FaceData<gc::Vector2> aligned = surface::computeCurvatureAlignedFaceDirectionField(*geometry, 4);
    for (surface::Face face : mesh->faces()) {
        const int index = face.getIndex();
        const double theta = std::atan2(aligned[face].y, aligned[face].x) * 0.25;
        const gc::Vector3 basisX = geometry->faceTangentBasis[face][0];
        const gc::Vector3 basisY = geometry->faceTangentBasis[face][1];
        gc::Vector3 direction = basisX * std::cos(theta) + basisY * std::sin(theta);
        const double length = direction.norm();
        if (length > 1e-12) direction = direction * (1.0 / length);
        field.directions.row(index) << direction.x, direction.y, direction.z;
        const gc::Vector3 normal = geometry->faceNormals[face];
        field.normals.row(index) << normal.x, normal.y, normal.z;
    }
    return singularitiesFromDirections(field, singularities, error);
}

bool computePrincipal(QuadCover::PrincipalField& field, FieldPicture& singularities, std::string& error) {
    NRosyVectorFields rosy(gSurface, 4);
    if (!rosy.computeFromCurvature()) {
        error = "主曲率方向计算失败";
        return false;
    }
    if (rosy.getField().rows() != field.directions.rows()) {
        error = "主曲率方向和面数不一致";
        return false;
    }
    field.directions = rosy.getField();
    return singularitiesFromDirections(field, singularities, error);
}

bool computeTrivial(QuadCover::PrincipalField& field, FieldPicture& singularities, std::string& error) {
    NRosyVectorFields rosy(gSurface, 4);
    rosy.autoSingularities();
    std::vector<int> prescribed(static_cast<size_t>(rosy.numVertices()));
    for (int vertex = 0; vertex < rosy.numVertices(); ++vertex) prescribed[vertex] = rosy.getSingularityIndices()(vertex);
    if (!rosy.computeTrivialConnection(prescribed)) {
        error = "平凡联络求解失败";
        return false;
    }
    if (rosy.getField().rows() != field.directions.rows()) {
        error = "平凡联络和面数不一致";
        return false;
    }
    field.directions = rosy.getField();
    paintMeshSingularities(gSurface, prescribed, singularities);
    return true;
}

bool computePoly(QuadCover::PrincipalField& field, FieldPicture& singularities, std::string& error) {
    ComplexPolyField poly(gSurface, 4);
    poly.compute();
    if (poly.numFaces() != field.directions.rows()) {
        error = "复多项式场和面数不一致";
        return false;
    }
    for (int face = 0; face < poly.numFaces(); ++face) field.directions.row(face) = poly.direction(face, 0);
    return singularitiesFromDirections(field, singularities, error);
}

bool computeSelected(int algo, std::string& error) {
    QuadCover::PrincipalField field;
    FieldPicture singularities;
    try {
        if (!exportInterior(gSurface, field, error)) return false;
        bool ok = false;
        if (algo == 0) ok = computeAligned(field, singularities, error);
        else if (algo == 1) ok = computePrincipal(field, singularities, error);
        else if (algo == 2) ok = computeTrivial(field, singularities, error);
        else ok = computePoly(field, singularities, error);
        if (!ok) return false;
    } catch (const std::exception& exception) {
        error = exception.what();
        return false;
    }
    gField = std::move(field);
    gSingularities = std::move(singularities);
    gAlgo = algo;
    gCaption = kAlgos[gAlgo].caption;
    return true;
}

float spacingEdgesForDensity(float density) {
    return std::max(0.22f, 10.f / std::max(density, 0.5f));
}

void present(bool resetCamera) {
    if (!gReady) return;
    const int nFaces = static_cast<int>(gField.directions.rows());
    std::vector<glm::vec3> positions(gField.positions.rows());
    std::vector<std::array<size_t, 3>> triangles(nFaces);
    for (int i = 0; i < gField.positions.rows(); ++i) positions[i] = at(gField.positions, i);
    for (int face = 0; face < nFaces; ++face) {
        triangles[face] = {
            static_cast<size_t>(gField.triangles[face * 3]),
            static_cast<size_t>(gField.triangles[face * 3 + 1]),
            static_cast<size_t>(gField.triangles[face * 3 + 2]),
        };
    }
    FieldPicture picture = buildPicture(gField, spacingEdgesForDensity(gDensity));
    picture.singularPoints = gSingularities.singularPoints;
    picture.singularColors = gSingularities.singularColors;
    picture.positive = gSingularities.positive;
    picture.negative = gSingularities.negative;
    polyscope::removeAllStructures();
    gPsMesh = polyscope::registerSurfaceMesh(gModels[gModel].label, positions, triangles);
    gPsMesh->setSmoothShade(true);
    gPsMesh->setSurfaceColor({0.78, 0.78, 0.78});
    gPsMesh->setEdgeColor({0.15, 0.15, 0.15});
    gPsMesh->setEdgeWidth(gShowEdges ? 1.0 : 0.0);
    drawField(gPsMesh, picture);
    if (resetCamera) polyscope::view::resetCameraToHomeView();
    gStatus = gModels[gModel].label + "    " + kAlgos[gAlgo].name + "    " + std::to_string(nFaces) + " 面    流线 " +
              std::to_string(picture.majorCurves) + " / " + std::to_string(picture.minorCurves) + "    奇异点 +" +
              std::to_string(picture.positive) + " −" + std::to_string(picture.negative);
}

bool showModel(int index) {
    if (index < 0 || index >= static_cast<int>(gModels.size())) return false;
    if (!gSurface.read(gModels[index].path.string())) {
        gStatus = "无法读取 " + gModels[index].label;
        return false;
    }
    std::string error;
    if (!computeSelected(gAlgo, error)) {
        gStatus = "无法计算 " + gModels[index].label + (error.empty() ? "" : "  " + error);
        return false;
    }
    gModel = index;
    gReady = true;
    present(true);
    return true;
}

void addIntegerSegment(const Eigen::Vector2d corners[3], const Eigen::Vector3d positions[3], const Eigen::Vector3d& lift,
                      std::vector<glm::vec3>& nodes, std::vector<std::array<size_t, 2>>& edges) {
    const int cornerOf[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int axis = 0; axis < 2; ++axis) {
        const double value[3] = {corners[0][axis], corners[1][axis], corners[2][axis]};
        if (!std::isfinite(value[0]) || !std::isfinite(value[1]) || !std::isfinite(value[2])) return;
        const double lo = std::min(value[0], std::min(value[1], value[2]));
        const double hi = std::max(value[0], std::max(value[1], value[2]));
        const int first = static_cast<int>(std::floor(lo)) + 1;
        const int last = static_cast<int>(std::ceil(hi)) - 1;
        if (last - first > 24) continue;
        for (int level = first; level <= last; ++level) {
            Eigen::Vector3d hit[2];
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
                hit[count++] = positions[a] * (1.0 - t) + positions[b] * t;
            }
            if (count == 2 && (hit[1] - hit[0]).norm() > 1e-8) {
                const size_t base = nodes.size();
                nodes.push_back(toGlm(hit[0] + lift));
                nodes.push_back(toGlm(hit[1] + lift));
                edges.push_back({base, base + 1});
            }
        }
    }
}

void drawQuadMesh(const std::vector<std::array<Eigen::Vector2d, 3>>& charts) {
    const int nFaces = static_cast<int>(gField.directions.rows());
    if (static_cast<int>(charts.size()) != nFaces) return;
    double edgeSum = 0;
    int edgeCount = 0;
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    for (int face = 0; face < nFaces; ++face) {
        Eigen::Vector3d position[3];
        for (int corner = 0; corner < 3; ++corner) {
            position[corner] = gField.positions.row(gField.triangles[face * 3 + corner]);
        }
        edgeSum += (position[1] - position[0]).norm();
        ++edgeCount;
        Eigen::Vector3d normal = gField.normals.row(face);
        if (normal.norm() > 1e-8) normal.normalize();
        const Eigen::Vector2d chart[3] = {charts[face][0], charts[face][1], charts[face][2]};
        addIntegerSegment(chart, position, normal * (edgeCount > 0 ? edgeSum / edgeCount : 0.01) * 0.04, nodes, edges);
    }
    polyscope::removeStructure("四边形网格", false);
    if (edges.empty()) {
        gStatus = std::string(gQuadMethod == 0 ? "MIQ" : "QuadCover") + " 没有诱导出四边形边";
        return;
    }
    polyscope::CurveNetwork* network = polyscope::registerCurveNetwork("四边形网格", nodes, edges);
    network->setColor({0.08f, 0.38f, 0.24f});
    network->setMaterial("flat");
    network->setRadius(0.00022f, true);
    network->setEnabled(true);
    gStatus = std::string(gQuadMethod == 0 ? "MIQ" : "QuadCover") + "    四边形边 " + std::to_string(edges.size());
}

void quadrangulate() {
    if (!gReady || gField.directions.rows() <= 0) return;
    const auto started = std::chrono::steady_clock::now();
    std::vector<std::array<Eigen::Vector2d, 3>> charts;
    try {
        const int nFaces = static_cast<int>(gField.directions.rows());
        if (gQuadMethod == 0) {
            MIQQuad miq(gSurface);
            miq.setCrossField(gField.directions);
            miq.parameterize();
            const Eigen::MatrixXd& corner = miq.cornerUV();
            if (corner.rows() != nFaces || corner.cols() < 6) {
                gStatus = "MIQ 没有得到参数坐标";
                return;
            }
            charts.resize(nFaces);
            for (int face = 0; face < nFaces; ++face) {
                for (int cornerIndex = 0; cornerIndex < 3; ++cornerIndex) {
                    charts[face][cornerIndex] = Eigen::Vector2d(corner(face, cornerIndex * 2), corner(face, cornerIndex * 2 + 1));
                }
            }
        } else {
            QuadCover cover(gSurface);
            cover.setCrossField(gField.directions);
            cover.setSmoothIterations(0);
            if (!cover.parameterizeFull()) {
                gStatus = "QuadCover 参数化失败";
                return;
            }
            const Eigen::MatrixXd& chart = cover.chartUV();
            if (chart.rows() < gField.positions.rows() || chart.cols() < 2) {
                gStatus = "QuadCover 没有得到参数坐标";
                return;
            }
            charts.resize(nFaces);
            for (int face = 0; face < nFaces; ++face) {
                for (int cornerIndex = 0; cornerIndex < 3; ++cornerIndex) {
                    const int vertex = gField.triangles[face * 3 + cornerIndex];
                    charts[face][cornerIndex] = chart.row(vertex);
                }
            }
        }
    } catch (const std::exception& exception) {
        gStatus = std::string("四边形化失败  ") + exception.what();
        return;
    }
    drawQuadMesh(charts);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    gStatus += "    " + std::to_string(seconds).substr(0, 4) + " 秒";
}

void callback() {
    ImGui::Text("向量场");
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
    if (ImGui::BeginCombo("算法", kAlgos[gAlgo].name)) {
        for (int algo = 0; algo < static_cast<int>(sizeof(kAlgos) / sizeof(kAlgos[0])); ++algo) {
            if (ImGui::Selectable(kAlgos[algo].name, algo == gAlgo) && algo != gAlgo) {
                std::string error;
                if (!computeSelected(algo, error)) gStatus = "算法失败  " + error;
                else present(false);
            }
        }
        ImGui::EndCombo();
    }
    if (ImGui::Checkbox("三角网格", &gShowEdges) && gPsMesh) gPsMesh->setEdgeWidth(gShowEdges ? 1.0 : 0.0);
    if (ImGui::SliderFloat("流线密度", &gDensity, 1.f, 40.f, "%.1f")) present(false);
    const char* quadNames[] = {"MIQ", "QuadCover"};
    if (ImGui::BeginCombo("四边形化", quadNames[gQuadMethod])) {
        for (int method = 0; method < 2; ++method) {
            if (ImGui::Selectable(quadNames[method], method == gQuadMethod)) gQuadMethod = method;
        }
        ImGui::EndCombo();
    }
    if (ImGui::Button("生成四边形网格") && gReady) quadrangulate();
    if (!gStatus.empty()) ImGui::TextWrapped("%s", gStatus.c_str());
    ImGui::TextWrapped("%s", gCaption.c_str());
}

int main(int argc, char** argv) {
    args::ArgumentParser parser("主曲率场");
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

    polyscope::options::programName = "主曲率场";
    polyscope::init();
    polyscope::state::userCallback = callback;
    if (!showModel(gModel)) return 1;
    polyscope::show();
    return 0;
}
