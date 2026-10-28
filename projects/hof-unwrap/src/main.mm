#define GL_SILENCE_DEPRECATION
#define GLFW_INCLUDE_NONE
#include <OpenGL/gl3.h>
#include <GLFW/glfw3.h>
#include <AppKit/AppKit.h>
#include <CoreText/CoreText.h>

#include "HolomorphicOneForm.h"
#include "Mesh.h"
#include "QcError.h"
#include "models.h"

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct DrawVertex {
    float px, py, pz, nx, ny, nz, u, v, r, g, b;
};

struct TextImage {
    GLuint texture = 0;
    int width = 1;
    int height = 1;
};

struct Button {
    float x = 0, y = 0, w = 0, h = 0;
    int id = 0;
    bool on = false;
    std::string label;
};

struct Mat4 {
    float m[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
};

struct Vec3 {
    float x = 0, y = 0, z = 0;
};

struct Panel {
    float x, y, w, h;
};

enum class Texture { Checker, Distortion };

Mesh gMesh;
std::vector<ModelFile> gModels;
int gModel = 0;
Texture gTexture = Texture::Checker;
bool gHasResult = false;
bool gHasTreeCotree = false;
bool gShowTree = false;
bool gShowCotree = false;
bool gShowHomology = false;
bool gShowCycles = false;
bool gDrag = false;
int gDragPanel = -1;
double gLastX = 0, gLastY = 0;
float gYaw = 0.7f, gPitch = 0.35f, gRadius = 3.2f;
float gUvZoom = 1.f, gUvPanX = 0.f, gUvPanY = 0.f;
Vec3 gCenter;
std::string gStatus = "选择模型后点「展开」。全纯 1-形式直接走 C++。";
double gSeconds = 0;
int gFlips = 0;
int gSeams = 0;
int gGenus = 0;
int gWidth = 1440, gHeight = 860;
float gScale = 1.f;
float gStatusX = 12.f;
float gStatusY = 86.f;
std::vector<Button> gButtons;
std::unordered_map<std::string, TextImage> gText;
std::vector<DrawVertex> gTriangles, gUvTriangles, gFlipTris, gFlips3d, gSeamLines;
std::vector<DrawVertex> gTreeLines, gCotreeLines, gHomologyLines, gCycleLines;

GLuint gProgram = 0, gFlatProgram = 0, gVao = 0, gBuffer = 0;

Mat4 multiply(const Mat4& a, const Mat4& b) {
    Mat4 out;
    for (int column = 0; column < 4; ++column) {
        for (int row = 0; row < 4; ++row) {
            out.m[column * 4 + row] = 0;
            for (int k = 0; k < 4; ++k) out.m[column * 4 + row] += a.m[k * 4 + row] * b.m[column * 4 + k];
        }
    }
    return out;
}

Mat4 perspective(float fov, float aspect, float nearPlane, float farPlane) {
    const float f = 1.f / std::tan(fov * 0.5f);
    Mat4 out{};
    out.m[0] = f / std::max(aspect, 1e-4f);
    out.m[5] = f;
    out.m[10] = (farPlane + nearPlane) / (nearPlane - farPlane);
    out.m[11] = -1.f;
    out.m[14] = (2.f * farPlane * nearPlane) / (nearPlane - farPlane);
    return out;
}

Mat4 lookAt(Vec3 eye, Vec3 center) {
    Vec3 forward{center.x - eye.x, center.y - eye.y, center.z - eye.z};
    const float fl = std::sqrt(forward.x * forward.x + forward.y * forward.y + forward.z * forward.z);
    forward.x /= fl; forward.y /= fl; forward.z /= fl;
    const Vec3 worldUp{0, 1, 0};
    Vec3 side{forward.y * worldUp.z - forward.z * worldUp.y, forward.z * worldUp.x - forward.x * worldUp.z,
              forward.x * worldUp.y - forward.y * worldUp.x};
    const float sl = std::sqrt(side.x * side.x + side.y * side.y + side.z * side.z);
    side.x /= sl; side.y /= sl; side.z /= sl;
    const Vec3 up{side.y * forward.z - side.z * forward.y, side.z * forward.x - side.x * forward.z,
                  side.x * forward.y - side.y * forward.x};
    Mat4 out;
    out.m[0] = side.x; out.m[4] = side.y; out.m[8] = side.z;
    out.m[12] = -(side.x * eye.x + side.y * eye.y + side.z * eye.z);
    out.m[1] = up.x; out.m[5] = up.y; out.m[9] = up.z;
    out.m[13] = -(up.x * eye.x + up.y * eye.y + up.z * eye.z);
    out.m[2] = -forward.x; out.m[6] = -forward.y; out.m[10] = -forward.z;
    out.m[14] = forward.x * eye.x + forward.y * eye.y + forward.z * eye.z;
    out.m[15] = 1;
    return out;
}

Mat4 ortho(float left, float right, float bottom, float top) {
    Mat4 out{};
    out.m[0] = 2.f / (right - left);
    out.m[5] = 2.f / (top - bottom);
    out.m[10] = -1.f;
    out.m[12] = -(right + left) / (right - left);
    out.m[13] = -(top + bottom) / (top - bottom);
    out.m[15] = 1.f;
    return out;
}

Vec3 eyePosition() {
    const float cp = std::cos(gPitch), sp = std::sin(gPitch);
    return {gCenter.x + gRadius * cp * std::sin(gYaw), gCenter.y + gRadius * sp,
            gCenter.z + gRadius * cp * std::cos(gYaw)};
}

void fitCamera() {
    if (gMesh.vertices.empty()) return;
    Eigen::Vector3d lo = gMesh.vertices[0].position, hi = lo, sum = Eigen::Vector3d::Zero();
    for (VertexCIter vertex = gMesh.vertices.begin(); vertex != gMesh.vertices.end(); ++vertex) {
        lo = lo.cwiseMin(vertex->position);
        hi = hi.cwiseMax(vertex->position);
        sum += vertex->position;
    }
    const Eigen::Vector3d mid = sum / std::max(1, static_cast<int>(gMesh.vertices.size()));
    gCenter = {static_cast<float>(mid.x()), static_cast<float>(mid.y()), static_cast<float>(mid.z())};
    gRadius = static_cast<float>((hi - lo).norm()) * 1.15f;
    if (gRadius < 0.2f) gRadius = 2.f;
}

int meshGenus() {
    int faces = 0;
    for (FaceCIter face = gMesh.faces.begin(); face != gMesh.faces.end(); ++face) {
        if (!face->isBoundary()) ++faces;
    }
    const int chi = static_cast<int>(gMesh.vertices.size()) - static_cast<int>(gMesh.edges.size()) + faces;
    if (!gMesh.boundaries.empty()) return 0;
    return std::max(0, (2 - chi) / 2);
}

TextImage rasterText(const std::string& text, float pixelHeight, bool dark) {
    const std::string key = text + "|" + std::to_string(static_cast<int>(pixelHeight)) + (dark ? "|d" : "|w");
    const auto found = gText.find(key);
    if (found != gText.end()) return found->second;
    TextImage image;
    @autoreleasepool {
        NSString* string = [NSString stringWithUTF8String:text.c_str()];
        if (string == nil) string = @"";
        NSFont* font = [NSFont systemFontOfSize:pixelHeight];
        NSColor* color = dark ? [NSColor colorWithCalibratedWhite:0.12 alpha:1.0]
                              : [NSColor colorWithCalibratedWhite:0.95 alpha:1.0];
        NSDictionary* attributes = @{NSFontAttributeName: font, NSForegroundColorAttributeName: color};
        const NSSize measured = [string sizeWithAttributes:attributes];
        image.width = std::max(1, static_cast<int>(std::ceil(measured.width)) + 4);
        image.height = std::max(1, static_cast<int>(std::ceil(measured.height)) + 4);
        std::vector<unsigned char> pixels(static_cast<size_t>(image.width * image.height * 4), 0);
        CGColorSpaceRef space = CGColorSpaceCreateDeviceRGB();
        CGContextRef context = CGBitmapContextCreate(pixels.data(), image.width, image.height, 8, image.width * 4, space,
                                                     kCGImageAlphaPremultipliedLast);
        CGContextTranslateCTM(context, 0, image.height);
        CGContextScaleCTM(context, 1, -1);
        NSGraphicsContext* graphics = [NSGraphicsContext graphicsContextWithCGContext:context flipped:YES];
        [NSGraphicsContext saveGraphicsState];
        [NSGraphicsContext setCurrentContext:graphics];
        [string drawAtPoint:NSMakePoint(2, 2) withAttributes:attributes];
        [NSGraphicsContext restoreGraphicsState];
        CGContextRelease(context);
        CGColorSpaceRelease(space);
        glGenTextures(1, &image.texture);
        glBindTexture(GL_TEXTURE_2D, image.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    }
    gText[key] = image;
    return image;
}

void pushVertex(std::vector<DrawVertex>& target, const Eigen::Vector3d& p, const Eigen::Vector3d& n,
                const Eigen::Vector2d& uv, const Eigen::Vector3d& color) {
    target.push_back(DrawVertex{static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()),
                                 static_cast<float>(n.x()), static_cast<float>(n.y()), static_cast<float>(n.z()),
                                 static_cast<float>(uv.x()), static_cast<float>(uv.y()),
                                 static_cast<float>(color.x()), static_cast<float>(color.y()), static_cast<float>(color.z())});
}

bool negativeUv(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c) {
    return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()) < 0;
}

void rebuildFromMesh(const Mesh& mesh, bool useUv) {
    gTriangles.clear();
    gUvTriangles.clear();
    gFlipTris.clear();
    gFlips3d.clear();
    gFlips = 0;
    std::vector<Eigen::Vector3d> distortion(mesh.vertices.size(), Eigen::Vector3d(0.72, 0.78, 0.84));
    if (useUv && gTexture == Texture::Distortion) {
        std::vector<double> weight(mesh.vertices.size(), 0);
        std::vector<Eigen::Vector3d> color(mesh.vertices.size(), Eigen::Vector3d::Zero());
        for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
            if (face->isBoundary()) continue;
            HalfEdgeCIter he = face->he;
            VertexCIter verts[3] = {he->vertex, he->next->vertex, he->next->next->vertex};
            std::vector<Eigen::Vector3d> p = {verts[0]->position, verts[1]->position, verts[2]->position};
            std::vector<Eigen::Vector3d> q;
            for (int corner = 0; corner < 3; ++corner) q.emplace_back(verts[corner]->uv.x(), verts[corner]->uv.y(), 0);
            const Eigen::Vector3d faceColor = QuasiConformalError::color(QuasiConformalError::compute(p, q));
            for (int corner = 0; corner < 3; ++corner) {
                if (verts[corner]->index < 0 || verts[corner]->index >= static_cast<int>(color.size())) continue;
                color[verts[corner]->index] += faceColor;
                weight[verts[corner]->index] += 1;
            }
        }
        for (size_t index = 0; index < color.size(); ++index) {
            if (weight[index] > 0) distortion[index] = color[index] / weight[index];
        }
    }
    for (FaceCIter face = mesh.faces.begin(); face != mesh.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        HalfEdgeCIter he = face->he;
        VertexCIter verts[3] = {he->vertex, he->next->vertex, he->next->next->vertex};
        Eigen::Vector3d positions[3] = {verts[0]->position, verts[1]->position, verts[2]->position};
        Eigen::Vector3d normal = (positions[1] - positions[0]).cross(positions[2] - positions[0]);
        if (normal.norm() > 1e-12) normal.normalize();
        else normal = Eigen::Vector3d(0, 1, 0);
        const bool flipped = useUv && negativeUv(verts[0]->uv, verts[1]->uv, verts[2]->uv);
        if (flipped) ++gFlips;
        for (int corner = 0; corner < 3; ++corner) {
            Eigen::Vector3d color(0.72, 0.78, 0.84);
            if (useUv && gTexture == Texture::Distortion && verts[corner]->index >= 0 &&
                verts[corner]->index < static_cast<int>(distortion.size())) {
                color = distortion[verts[corner]->index];
            }
            const Eigen::Vector2d uv = useUv ? verts[corner]->uv : Eigen::Vector2d::Zero();
            pushVertex(gTriangles, positions[corner], normal, uv, color);
            if (useUv) {
                pushVertex(gUvTriangles, Eigen::Vector3d(uv.x(), uv.y(), 0), Eigen::Vector3d(0, 0, 1), uv, color);
                if (flipped) {
                    pushVertex(gFlipTris, Eigen::Vector3d(uv.x(), uv.y(), 0.001), Eigen::Vector3d(0, 0, 1), uv,
                               Eigen::Vector3d(0.95, 0.78, 0.15));
                    pushVertex(gFlips3d, positions[corner] + normal * 0.002, normal, uv, Eigen::Vector3d(0.95, 0.78, 0.15));
                }
            }
        }
    }
}

void appendVertexPairs(std::vector<DrawVertex>& lines, const std::vector<std::pair<int, int>>& pairs,
                      const Eigen::Vector3d& color, float liftAmount) {
    std::vector<Eigen::Vector3d> position(gMesh.vertices.size(), Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> normal(gMesh.vertices.size(), Eigen::Vector3d(0, 1, 0));
    std::vector<char> known(gMesh.vertices.size(), 0);
    for (VertexCIter vertex = gMesh.vertices.begin(); vertex != gMesh.vertices.end(); ++vertex) {
        if (vertex->index < 0 || vertex->index >= static_cast<int>(position.size())) continue;
        position[static_cast<size_t>(vertex->index)] = vertex->position;
        Eigen::Vector3d n = vertex->normal();
        if (n.norm() > 1e-8) n.normalize();
        else n = Eigen::Vector3d(0, 1, 0);
        normal[static_cast<size_t>(vertex->index)] = n;
        known[static_cast<size_t>(vertex->index)] = 1;
    }
    for (const std::pair<int, int>& pair : pairs) {
        if (pair.first < 0 || pair.second < 0 ||
            pair.first >= static_cast<int>(known.size()) || pair.second >= static_cast<int>(known.size()) ||
            !known[static_cast<size_t>(pair.first)] || !known[static_cast<size_t>(pair.second)]) {
            continue;
        }
        const Eigen::Vector3d sum = normal[static_cast<size_t>(pair.first)] + normal[static_cast<size_t>(pair.second)];
        const Eigen::Vector3d n = sum.norm() > 1e-8 ? (sum / sum.norm()).eval() : Eigen::Vector3d(0, 1, 0);
        const Eigen::Vector3d lift = n * static_cast<double>(liftAmount);
        pushVertex(lines, position[static_cast<size_t>(pair.first)] + lift, n, Eigen::Vector2d::Zero(), color);
        pushVertex(lines, position[static_cast<size_t>(pair.second)] + lift, n, Eigen::Vector2d::Zero(), color);
    }
}

void clearTreeCotreeView() {
    gHasTreeCotree = false;
    gShowTree = false;
    gShowCotree = false;
    gShowHomology = false;
    gShowCycles = false;
    gTreeLines.clear();
    gCotreeLines.clear();
    gHomologyLines.clear();
    gCycleLines.clear();
}

void captureTreeCotree(const HolomorphicOneForm::TreeCotree& treeCotree) {
    gTreeLines.clear();
    gCotreeLines.clear();
    gHomologyLines.clear();
    gCycleLines.clear();
    gHasTreeCotree = !treeCotree.primalTree.empty() || !treeCotree.dualTree.empty() || !treeCotree.homologyEdges.empty();
    if (!gHasTreeCotree) return;
    appendVertexPairs(gTreeLines, treeCotree.primalTree, Eigen::Vector3d(0.86, 0.72, 0.13), 0.003f);
    appendVertexPairs(gCotreeLines, treeCotree.dualTree, Eigen::Vector3d(0.12, 0.78, 0.86), 0.004f);
    appendVertexPairs(gHomologyLines, treeCotree.homologyEdges, Eigen::Vector3d(0.35, 0.86, 0.38), 0.006f);
    std::vector<std::pair<int, int>> cycleEdges;
    for (const std::vector<std::pair<int, int>>& cycle : treeCotree.cycles) {
        cycleEdges.insert(cycleEdges.end(), cycle.begin(), cycle.end());
    }
    appendVertexPairs(gCycleLines, cycleEdges, Eigen::Vector3d(0.93, 0.32, 0.55), 0.008f);
}

void collectSeams(const CutSeamMesh& cut) {
    gSeamLines.clear();
    gSeams = 0;
    for (EdgeCIter edge = cut.edges.begin(); edge != cut.edges.end(); ++edge) {
        if (!cut.isCutBoundaryEdge(edge)) continue;
        const Eigen::Vector3d a = edge->he->vertex->position;
        const Eigen::Vector3d b = edge->he->flip->vertex->position;
        Eigen::Vector3d n = edge->he->vertex->normal() + edge->he->flip->vertex->normal();
        if (n.norm() > 1e-8) n.normalize();
        else n = Eigen::Vector3d(0, 1, 0);
        const Eigen::Vector3d lift = n * 0.004;
        const Eigen::Vector3d color(0.95, 0.25, 0.45);
        pushVertex(gSeamLines, a + lift, n, Eigen::Vector2d::Zero(), color);
        pushVertex(gSeamLines, b + lift, n, Eigen::Vector2d::Zero(), color);
        ++gSeams;
    }
}

void rebuildGeometry() {
    rebuildFromMesh(gMesh, false);
    gSeamLines.clear();
    gFlips = 0;
}

void loadModel(int index) {
    if (gModels.empty()) return;
    gModel = (index % static_cast<int>(gModels.size()) + static_cast<int>(gModels.size())) % static_cast<int>(gModels.size());
    gMesh.read(gModels[gModel].path.string());
    gHasResult = false;
    clearTreeCotreeView();
    gSeconds = 0;
    gFlips = 0;
    gSeams = 0;
    gGenus = meshGenus();
    gUvZoom = 1.f;
    gUvPanX = gUvPanY = 0.f;
    fitCamera();
    std::ostringstream text;
    text << gModels[gModel].label << " · 亏格 " << gGenus << " · " << gMesh.faces.size()
         << " 面。点「展开」才计算全纯 1-形式。";
    gStatus = text.str();
    rebuildGeometry();
}

void unfold() {
    if (gMesh.vertices.empty()) return;
    const auto started = std::chrono::steady_clock::now();
    std::string error;
    CutSeamMesh* cut = nullptr;
    try {
        HolomorphicOneForm solver(gMesh);
        solver.parameterize();
        captureTreeCotree(solver.treeCotree());
        cut = &solver.cutMesh();
        int placed = 0;
        double minU = 1e9, maxU = -1e9, minV = 1e9, maxV = -1e9;
        for (VertexCIter vertex = cut->vertices.begin(); vertex != cut->vertices.end(); ++vertex) {
            if (!std::isfinite(vertex->uv.x()) || !std::isfinite(vertex->uv.y())) continue;
            ++placed;
            minU = std::min(minU, vertex->uv.x());
            maxU = std::max(maxU, vertex->uv.x());
            minV = std::min(minV, vertex->uv.y());
            maxV = std::max(maxV, vertex->uv.y());
        }
        gHasResult = placed > 2 && (maxU - minU) > 1e-6 && (maxV - minV) > 1e-6;
        if (!gHasResult) {
            if (static_cast<int>(cut->edges.size()) > 12000) {
                error = "割缝后的边数超过 12000，全纯 1-形式跳过了求解。换一个更小的模型。";
            } else {
                error = "全纯 1-形式没有得到坐标。tree-cotree 割缝或周期约束没有给出可用的 1-形式。";
            }
        } else {
            rebuildFromMesh(*cut, true);
            collectSeams(*cut);
        }
    } catch (const std::exception& exception) {
        error = exception.what();
        gHasResult = false;
    }
    gSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    if (!error.empty()) {
        if (gHasTreeCotree) error += " 第二行可以看这次的 tree-cotree。";
        gStatus = error;
        rebuildGeometry();
        return;
    }
    std::ostringstream text;
    text.setf(std::ios::fixed);
    text.precision(3);
    text << "全纯 1-形式 · 亏格 " << gGenus << " · " << gSeconds << " 秒 · 翻转 " << gFlips
         << " · 割缝 " << gSeams << "。第二行是这次的中间结果。";
    gStatus = text.str();
}

float toolbarHeight() { return 118.f; }

Panel leftPanel() {
    const float gap = 8.f;
    const float width = (static_cast<float>(gWidth) - 3.f * gap) * 0.5f;
    return {gap, toolbarHeight(), width, static_cast<float>(gHeight) - toolbarHeight() - gap};
}

Panel rightPanel() {
    const Panel left = leftPanel();
    return {left.x + left.w + 8.f, left.y, left.w, left.h};
}

Mat4 modelViewProjection(const Panel& panel) {
    const float aspect = std::max(panel.w, 1.f) / std::max(panel.h, 1.f);
    return multiply(perspective(0.9f, aspect, 0.01f, 100.f), lookAt(eyePosition(), gCenter));
}

void addButton(float& cursor, float y, float height, int id, const std::string& label, bool on) {
    const TextImage image = rasterText(label, 15.f * gScale, true);
    const float width = std::max(28.f, image.width / gScale + 16.f);
    gButtons.push_back(Button{cursor, y, width, height, id, on, label});
    cursor += width + 6.f;
}

void layoutButtons() {
    gButtons.clear();
    float cursor = 12.f;
    const float y = 12.f;
    const float height = 28.f;
    addButton(cursor, y, height, 1, "<", false);
    addButton(cursor, y, height, 0, gModels.empty() ? "没有模型" : gModels[gModel].label, false);
    addButton(cursor, y, height, 2, ">", false);
    cursor += 8.f;
    addButton(cursor, y, height, 40, "展开", false);
    addButton(cursor, y, height, 41, "清除", false);
    addButton(cursor, y, height, 50, "棋盘格", gTexture == Texture::Checker);
    addButton(cursor, y, height, 51, "共形畸变", gTexture == Texture::Distortion);
    cursor = 12.f;
    const float row2 = 46.f;
    addButton(cursor, row2, height, 60, "生成树", gShowTree);
    addButton(cursor, row2, height, 61, "余树", gShowCotree);
    addButton(cursor, row2, height, 62, "同调边", gShowHomology);
    addButton(cursor, row2, height, 63, "基本环", gShowCycles);
    gStatusX = 12.f;
    gStatusY = 84.f;
}

void recolorDistortion() {
    if (!gHasResult || gTriangles.size() != gUvTriangles.size() || gTriangles.size() % 3 != 0) return;
    for (size_t face = 0; face < gTriangles.size(); face += 3) {
        std::vector<Eigen::Vector3d> positions, uvs;
        for (int corner = 0; corner < 3; ++corner) {
            const DrawVertex& vertex = gTriangles[face + corner];
            const DrawVertex& uv = gUvTriangles[face + corner];
            positions.emplace_back(vertex.px, vertex.py, vertex.pz);
            uvs.emplace_back(uv.px, uv.py, 0);
        }
        const Eigen::Vector3d color = QuasiConformalError::color(QuasiConformalError::compute(positions, uvs));
        for (int corner = 0; corner < 3; ++corner) {
            DrawVertex& vertex = gTriangles[face + corner];
            DrawVertex& uv = gUvTriangles[face + corner];
            vertex.r = uv.r = static_cast<float>(color.x());
            vertex.g = uv.g = static_cast<float>(color.y());
            vertex.b = uv.b = static_cast<float>(color.z());
        }
    }
}

void clickButton(int id) {
    if (id == 1) loadModel(gModel - 1);
    else if (id == 2) loadModel(gModel + 1);
    else if (id == 40) unfold();
    else if (id == 41) {
        gHasResult = false;
        gSeconds = 0;
        gFlips = 0;
        gSeams = 0;
        clearTreeCotreeView();
        gStatus = "已清除展开结果。点「展开」再算全纯 1-形式。";
        rebuildGeometry();
    } else if (id == 50) {
        gTexture = Texture::Checker;
        if (!gHasResult) rebuildGeometry();
    } else if (id == 51) {
        gTexture = Texture::Distortion;
        if (gHasResult) recolorDistortion();
        else rebuildGeometry();
    } else if (id >= 60 && id <= 63) {
        if (!gHasTreeCotree) {
            gStatus = "先点「展开」。生成树、余树、同调边和基本环来自这一次计算。";
            return;
        }
        if (id == 60) gShowTree = !gShowTree;
        else if (id == 61) gShowCotree = !gShowCotree;
        else if (id == 62) gShowHomology = !gShowHomology;
        else gShowCycles = !gShowCycles;
    }
}

GLuint compileShader(GLenum type, const char* source) {
    const GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    return shader;
}

void setupGL() {
    const char* vertexSource = R"(#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec2 aUv;
layout(location=3) in vec3 aCol;
uniform mat4 uMvp;
out vec3 vNrm; out vec2 vUv; out vec3 vCol;
void main(){ gl_Position = uMvp * vec4(aPos,1.0); vNrm=aNrm; vUv=aUv; vCol=aCol; })";
    const char* fragmentSource = R"(#version 330 core
in vec3 vNrm; in vec2 vUv; in vec3 vCol; out vec4 frag;
uniform int uShade; uniform vec3 uLight;
void main(){
 vec3 n = normalize(vNrm);
 float light = 0.42 + 0.58*max(dot(n, normalize(uLight)), 0.0);
 vec3 color = vCol;
 if(uShade==0){
   float cell = step(0.5, fract(vUv.x*8.0)) == step(0.5, fract(vUv.y*8.0)) ? 0.92 : 0.28;
   color = vec3(cell);
 }
 if(uShade==3) color = vec3(0.95,0.78,0.15);
 if(uShade==4){ frag = vec4(vCol, 1.0); return; }
 frag = vec4(color * light, 1.0);
})";
    gProgram = glCreateProgram();
    glAttachShader(gProgram, compileShader(GL_VERTEX_SHADER, vertexSource));
    glAttachShader(gProgram, compileShader(GL_FRAGMENT_SHADER, fragmentSource));
    glLinkProgram(gProgram);
    const char* flatVertex = R"(#version 330 core
layout(location=0) in vec2 aPos; layout(location=1) in vec2 aUv; uniform vec2 uScreen; out vec2 vUv;
void main(){ vUv=aUv; gl_Position = vec4(aPos.x/uScreen.x*2.0-1.0, 1.0-aPos.y/uScreen.y*2.0, 0.0, 1.0); })";
    const char* flatFragment = R"(#version 330 core
uniform vec4 uColor; uniform sampler2D uTex; uniform int uTextured; in vec2 vUv; out vec4 frag;
void main(){ vec4 sample = texture(uTex, vUv); frag = uTextured==1 ? sample : uColor; })";
    gFlatProgram = glCreateProgram();
    glAttachShader(gFlatProgram, compileShader(GL_VERTEX_SHADER, flatVertex));
    glAttachShader(gFlatProgram, compileShader(GL_FRAGMENT_SHADER, flatFragment));
    glLinkProgram(gFlatProgram);
    glGenVertexArrays(1, &gVao);
    glGenBuffers(1, &gBuffer);
}

void drawMesh(const std::vector<DrawVertex>& vertices, const Mat4& mvp, int shade, GLenum mode) {
    if (vertices.empty()) return;
    glUseProgram(gProgram);
    glUniformMatrix4fv(glGetUniformLocation(gProgram, "uMvp"), 1, GL_FALSE, mvp.m);
    glUniform1i(glGetUniformLocation(gProgram, "uShade"), shade);
    const Vec3 eye = eyePosition();
    glUniform3f(glGetUniformLocation(gProgram, "uLight"), eye.x - gCenter.x, eye.y - gCenter.y, eye.z - gCenter.z);
    glBindVertexArray(gVao);
    glBindBuffer(GL_ARRAY_BUFFER, gBuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(DrawVertex), vertices.data(), GL_DYNAMIC_DRAW);
    const GLsizei stride = sizeof(DrawVertex);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(12));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(24));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(32));
    glDrawArrays(mode, 0, static_cast<GLsizei>(vertices.size()));
}

void fillRect(float x, float y, float w, float h, float r, float g, float b) {
    const float data[] = {x, y, 0, 0, x + w, y, 1, 0, x + w, y + h, 1, 1, x, y, 0, 0, x + w, y + h, 1, 1, x, y + h, 0, 1};
    glUseProgram(gFlatProgram);
    glUniform2f(glGetUniformLocation(gFlatProgram, "uScreen"), static_cast<float>(gWidth), static_cast<float>(gHeight));
    glUniform4f(glGetUniformLocation(gFlatProgram, "uColor"), r, g, b, 1);
    glUniform1i(glGetUniformLocation(gFlatProgram, "uTextured"), 0);
    glBindVertexArray(gVao);
    glBindBuffer(GL_ARRAY_BUFFER, gBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 16, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 16, reinterpret_cast<void*>(8));
    glDrawArrays(GL_TRIANGLES, 0, 6);
}

void drawLabel(const std::string& text, float x, float y) {
    const TextImage image = rasterText(text, 15.f * gScale, true);
    const float w = image.width / gScale, h = image.height / gScale;
    const float data[] = {x, y, 0, 0, x + w, y, 1, 0, x + w, y + h, 1, 1, x, y, 0, 0, x + w, y + h, 1, 1, x, y + h, 0, 1};
    glUseProgram(gFlatProgram);
    glUniform2f(glGetUniformLocation(gFlatProgram, "uScreen"), static_cast<float>(gWidth), static_cast<float>(gHeight));
    glUniform1i(glGetUniformLocation(gFlatProgram, "uTextured"), 1);
    glUniform1i(glGetUniformLocation(gFlatProgram, "uTex"), 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, image.texture);
    glBindVertexArray(gVao);
    glBindBuffer(GL_ARRAY_BUFFER, gBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 16, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 16, reinterpret_cast<void*>(8));
    glDrawArrays(GL_TRIANGLES, 0, 6);
}

Mat4 uvMatrix(const Panel& panel) {
    float minU = 0, minV = 0, maxU = 1, maxV = 1;
    if (!gUvTriangles.empty()) {
        minU = maxU = gUvTriangles[0].px;
        minV = maxV = gUvTriangles[0].py;
        for (const DrawVertex& vertex : gUvTriangles) {
            minU = std::min(minU, vertex.px); maxU = std::max(maxU, vertex.px);
            minV = std::min(minV, vertex.py); maxV = std::max(maxV, vertex.py);
        }
    }
    const float cx = (minU + maxU) * 0.5f + gUvPanX;
    const float cy = (minV + maxV) * 0.5f + gUvPanY;
    const float span = std::max(maxU - minU, maxV - minV) * 0.62f / std::max(gUvZoom, 0.05f);
    const float aspect = panel.w / std::max(panel.h, 1.f);
    return ortho(cx - span * aspect, cx + span * aspect, cy - span, cy + span);
}

void drawFrame() {
    layoutButtons();
    glViewport(0, 0, static_cast<int>(gWidth * gScale), static_cast<int>(gHeight * gScale));
    glClearColor(0.93f, 0.94f, 0.95f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    fillRect(0, 0, static_cast<float>(gWidth), toolbarHeight(), 0.97f, 0.97f, 0.98f);
    for (const Button& button : gButtons) {
        fillRect(button.x, button.y, button.w, button.h, button.on ? 0.78f : 0.88f, button.on ? 0.86f : 0.89f, button.on ? 0.95f : 0.91f);
        drawLabel(button.label, button.x + 8.f, button.y + 5.f);
    }
    drawLabel(gStatus, gStatusX, gStatusY);

    glEnable(GL_DEPTH_TEST);
    const Panel left = leftPanel();
    const Panel right = rightPanel();
    auto viewportOf = [](const Panel& panel) {
        glViewport(static_cast<int>(panel.x * gScale), static_cast<int>((gHeight - panel.y - panel.h) * gScale),
                   static_cast<int>(panel.w * gScale), static_cast<int>(panel.h * gScale));
        glScissor(static_cast<int>(panel.x * gScale), static_cast<int>((gHeight - panel.y - panel.h) * gScale),
                  static_cast<int>(panel.w * gScale), static_cast<int>(panel.h * gScale));
    };
    glEnable(GL_SCISSOR_TEST);
    viewportOf(left);
    glClearColor(0.86f, 0.89f, 0.92f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    const int shade = gHasResult && gTexture == Texture::Checker ? 0 : 2;
    const Mat4 mvp = modelViewProjection(left);
    drawMesh(gTriangles, mvp, shade, GL_TRIANGLES);
    drawMesh(gFlips3d, mvp, 3, GL_TRIANGLES);
    if (gShowCotree) drawMesh(gCotreeLines, mvp, 4, GL_LINES);
    if (gShowTree) drawMesh(gTreeLines, mvp, 4, GL_LINES);
    if (gShowHomology) drawMesh(gHomologyLines, mvp, 4, GL_LINES);
    if (gShowCycles) drawMesh(gCycleLines, mvp, 4, GL_LINES);
    drawMesh(gSeamLines, mvp, 4, GL_LINES);
    viewportOf(right);
    glClearColor(0.96f, 0.96f, 0.94f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (gHasResult) {
        const Mat4 uv = uvMatrix(right);
        drawMesh(gUvTriangles, uv, shade, GL_TRIANGLES);
        glDisable(GL_DEPTH_TEST);
        drawMesh(gFlipTris, uv, 3, GL_TRIANGLES);
        glEnable(GL_DEPTH_TEST);
    }
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_DEPTH_TEST);
    glViewport(0, 0, static_cast<int>(gWidth * gScale), static_cast<int>(gHeight * gScale));
    drawLabel("金生成树 · 青余树 · 绿同调边 · 品红基本环 · 粉红割缝", left.x + 10.f, left.y + 8.f);
    drawLabel(gHasResult ? "UV" : "UV · 点展开", right.x + 10.f, right.y + 8.f);
}

void mouseButton(GLFWwindow*, int button, int action, int) {
    if (button != GLFW_MOUSE_BUTTON_LEFT) return;
    if (action == GLFW_RELEASE) { gDrag = false; return; }
    glfwGetCursorPos(glfwGetCurrentContext(), &gLastX, &gLastY);
    for (const Button& item : gButtons) {
        if (gLastX >= item.x && gLastX <= item.x + item.w && gLastY >= item.y && gLastY <= item.y + item.h) {
            clickButton(item.id);
            return;
        }
    }
    const Panel left = leftPanel();
    const Panel right = rightPanel();
    if (gLastX >= left.x && gLastX <= left.x + left.w && gLastY >= left.y && gLastY <= left.y + left.h) {
        gDrag = true; gDragPanel = 0;
    } else if (gLastX >= right.x && gLastX <= right.x + right.w && gLastY >= right.y && gLastY <= right.y + right.h) {
        gDrag = true; gDragPanel = 1;
    }
}

void cursor(GLFWwindow*, double x, double y) {
    if (!gDrag) return;
    const double dx = x - gLastX, dy = y - gLastY;
    gLastX = x; gLastY = y;
    if (gDragPanel == 0) {
        gYaw += static_cast<float>(dx) * 0.008f;
        gPitch = std::clamp(gPitch + static_cast<float>(dy) * 0.008f, -1.2f, 1.2f);
    } else {
        gUvPanX -= static_cast<float>(dx) * 0.002f / std::max(gUvZoom, 0.05f);
        gUvPanY += static_cast<float>(dy) * 0.002f / std::max(gUvZoom, 0.05f);
    }
}

void scroll(GLFWwindow*, double, double yoff) {
    double x = 0, y = 0;
    glfwGetCursorPos(glfwGetCurrentContext(), &x, &y);
    const Panel right = rightPanel();
    if (x >= right.x && x <= right.x + right.w && y >= right.y && y <= right.y + right.h) {
        gUvZoom = std::clamp(gUvZoom * (yoff > 0 ? 1.1f : 0.9f), 0.2f, 40.f);
    } else {
        gRadius = std::clamp(gRadius * (yoff > 0 ? 0.9f : 1.1f), 0.2f, 40.f);
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc > 2 && std::string(argv[1]) == "--check") {
        if (!gMesh.read(argv[2])) return 2;
        gGenus = meshGenus();
        unfold();
        std::cout << gStatus << std::endl;
        std::cout << "tree " << gTreeLines.size() / 2
                  << " cotree " << gCotreeLines.size() / 2
                  << " homology " << gHomologyLines.size() / 2
                  << " cycle-edges " << gCycleLines.size() / 2 << std::endl;
        return gHasResult ? 0 : 1;
    }
    gModels = listModels();
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
    GLFWwindow* window = glfwCreateWindow(gWidth, gHeight, "全纯1-形式", nullptr, nullptr);
    if (!window) return 1;
    glfwMakeContextCurrent(window);
    glfwSetMouseButtonCallback(window, mouseButton);
    glfwSetCursorPosCallback(window, cursor);
    glfwSetScrollCallback(window, scroll);
    glfwSwapInterval(1);
    setupGL();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if (!gModels.empty()) loadModel(0);
    else gStatus = "没有找到模型。";
    while (!glfwWindowShouldClose(window)) {
        int framebufferWidth = 0, framebufferHeight = 0;
        glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
        glfwGetWindowSize(window, &gWidth, &gHeight);
        gScale = gWidth > 0 ? static_cast<float>(framebufferWidth) / static_cast<float>(gWidth) : 1.f;
        drawFrame();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
