/**
 * Standalone demo: Abel–Jacobi map and divisor check on a closed mesh.
 *
 * Build (from repo cpp/):
 *   g++ -std=c++17 -O2 -I conformal-parameterization -I deps/eigen-3.4.0 \
 *       conformal-parameterization/Abel_Jacoi/abel_jacobi_demo.cpp \
 *       conformal-parameterization/Abel_Jacoi/AbelJacobi.cpp \
 *       conformal-parameterization/Abel_Jacoi/AbelJacobiParameterization.cpp \
 *       conformal-parameterization/Mesh.cpp MeshIO.cpp Parameterization.cpp \
 *       conformal-parameterization/Vertex.cpp Edge.cpp Face.cpp HalfEdge.cpp \
 *       conformal-parameterization/Solver.cpp QcError.cpp \
 *       conformal-parameterization/uv_unwrap_simple/Lscm.cpp \
 *       -o abel_jacobi_demo
 */
#include "AbelJacobiParameterization.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>

static bool isClosed(const Mesh& m)
{
    if (!m.boundaries.empty()) return false;
    for (EdgeCIter e = m.edges.begin(); e != m.edges.end(); ++e) {
        if (e->isBoundary()) return false;
    }
    return true;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh.obj> [base_vertex]\n";
        return 1;
    }

    Mesh mesh;
    if (!mesh.read(argv[1])) {
        std::cerr << "Failed to read mesh: " << argv[1] << "\n";
        return 1;
    }

    if (!isClosed(mesh)) {
        std::cerr << "Mesh must be closed (no boundary) for Abel–Jacobi.\n";
        return 1;
    }

    AbelJacobi aj(mesh);
    if (!aj.build()) {
        std::cerr << "Abel–Jacobi build failed.\n";
        return 1;
    }

    std::cout << "亏格 g = " << aj.genus() << "\n";
    std::cout << "正规基 a-周期误差 = " << aj.canonicalError() << "\n";
    std::cout << "周期矩阵 Π = A^{-1} B\n" << aj.periodMatrix() << "\n";

    const int base = (argc >= 3) ? std::atoi(argv[2]) : 0;
    const int nV = static_cast<int>(mesh.vertices.size());
    for (int vi = 0; vi < std::min(nV, 5); ++vi) {
        const Eigen::MatrixXd mu = aj.abelJacobiMap(base, vi);
        std::cout << "mu(" << vi << ") mod Γ = " << mu.transpose() << "\n";
    }

    if (!aj.quantizeAndParameterize(base)) {
        std::cerr << "积分 φ_1 失败。\n";
        return 1;
    }
    std::cout << "已把 φ_1 的积分写入顶点 UV。\n";

    AbelJacobiParameterization param(mesh);
    std::vector<Eigen::Vector3d> quadLines;
    int cells = 0;
    if (!param.buildQuadLines(quadLines, nullptr, nullptr, &cells)) {
        std::cerr << "四边形网格提取失败。\n";
        return 1;
    }
    std::cout << "四边形网格段数 = " << (quadLines.size() / 2) << " 格子 = " << cells << "\n";

    const std::string out = "abel_jacobi_uv.obj";
    if (mesh.write(out)) {
        std::cout << "Wrote " << out << "\n";
    }
    return 0;
}
