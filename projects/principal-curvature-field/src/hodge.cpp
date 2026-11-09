#include "hodge.h"

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <queue>
#include <utility>

namespace {

constexpr double kPi = 3.14159265358979323846;

Eigen::Vector3d pointAt(const Eigen::MatrixXd& points, int row) {
    return Eigen::Vector3d(points(row, 0), points(row, 1), points(row, 2));
}

// 把 v 绕单位轴 axis 转 angle。
Eigen::Vector3d rotateAround(const Eigen::Vector3d& axis, const Eigen::Vector3d& v, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return v * c + axis.cross(v) * s + axis * (axis.dot(v)) * (1.0 - c);
}

// b 相对 a 绕 normal 的有符号夹角。
double angleAround(const Eigen::Vector3d& normal, const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return std::atan2(normal.dot(a.cross(b)), a.dot(b));
}

// 十字场每面给出的只是一个代表方向，面与面之间可能差 90° 的整数倍。
// 沿面片做并行输运（绕公共边转过去），每一步选最接近的分支，得到尽量连续的单个向量场。
std::vector<Eigen::Vector3d> trivialize(const std::vector<Eigen::Vector3d>& raw,
                                        const std::vector<Eigen::Vector3d>& normals,
                                        const Eigen::MatrixXd& positions,
                                        const Eigen::MatrixXi& faces,
                                        int fold) {
    const int nFaces = static_cast<int>(faces.rows());
    std::vector<Eigen::Vector3d> field(raw);
    std::vector<char> visited(static_cast<size_t>(nFaces), 0);
    if (fold <= 1) {
        std::fill(visited.begin(), visited.end(), 1);
        return field;
    }

    // 无向边 -> 相邻面
    std::map<std::pair<int, int>, std::vector<int>> edges;
    for (int face = 0; face < nFaces; ++face) {
        for (int a = 0; a < 3; ++a) {
            const int i = faces(face, a);
            const int j = faces(face, (a + 1) % 3);
            edges[{std::min(i, j), std::max(i, j)}].push_back(face);
        }
    }

    const double step = 2.0 * kPi / fold;
    std::queue<int> pending;
    for (int seed = 0; seed < nFaces; ++seed) {
        if (visited[seed]) continue;
        visited[seed] = 1;
        pending.push(seed);
        while (!pending.empty()) {
            const int face = pending.front();
            pending.pop();
            for (int a = 0; a < 3; ++a) {
                const int i = faces(face, a);
                const int j = faces(face, (a + 1) % 3);
                auto it = edges.find({std::min(i, j), std::max(i, j)});
                if (it == edges.end()) continue;
                for (int other : it->second) {
                    if (other == face || visited[other]) continue;
                    const Eigen::Vector3d edge = (pointAt(positions, j) - pointAt(positions, i)).normalized();
                    const Eigen::Vector3d source = field[face];
                    // 输运：沿边分量不变，面内垂直分量换到相邻面的对应方向
                    const Eigen::Vector3d carried =
                        edge * source.dot(edge) + normals[other].cross(edge) * source.dot(normals[face].cross(edge));
                    const double delta = angleAround(normals[other], carried, raw[other]);
                    const double turn = std::llround(delta / step) * step - delta;
                    field[other] = rotateAround(normals[other], raw[other], turn).normalized();
                    visited[other] = 1;
                    pending.push(other);
                }
            }
        }
    }
    return field;
}

}  // namespace

HodgeSplit decomposeHodge(const Eigen::MatrixXd& positions,
                          const Eigen::MatrixXi& faces,
                          const Eigen::MatrixXd& directions,
                          int fold) {
    HodgeSplit split;
    const int nVerts = static_cast<int>(positions.rows());
    const int nFaces = static_cast<int>(faces.rows());
    if (nVerts < 4 || nFaces < 2 || directions.rows() != nFaces) {
        split.note = "网格或向量场尺寸不匹配";
        return split;
    }

    std::vector<Eigen::Vector3d> normals(static_cast<size_t>(nFaces));
    std::vector<Eigen::Vector3d> raw(static_cast<size_t>(nFaces));
    std::vector<double> area(static_cast<size_t>(nFaces));
    for (int face = 0; face < nFaces; ++face) {
        const Eigen::Vector3d p0 = pointAt(positions, faces(face, 0));
        const Eigen::Vector3d p1 = pointAt(positions, faces(face, 1));
        const Eigen::Vector3d p2 = pointAt(positions, faces(face, 2));
        const Eigen::Vector3d cross = (p1 - p0).cross(p2 - p0);
        const double length = cross.norm();
        if (length < 1e-14) {
            split.note = "存在退化三角形";
            return split;
        }
        area[face] = 0.5 * length;
        normals[face] = cross / length;
        Eigen::Vector3d direction = pointAt(directions, face);
        direction -= normals[face] * direction.dot(normals[face]);  // 投影回切平面
        if (direction.norm() < 1e-12) {
            split.note = "向量场存在零向量";
            return split;
        }
        raw[face] = direction.normalized();
    }

    const std::vector<Eigen::Vector3d> field = trivialize(raw, normals, positions, faces, fold);

    // hat 基梯度：∇φ_i = n × (p_k - p_j) / (2A)
    std::vector<std::array<Eigen::Vector3d, 3>> hatGrad(static_cast<size_t>(nFaces));
    for (int face = 0; face < nFaces; ++face) {
        const Eigen::Vector3d p0 = pointAt(positions, faces(face, 0));
        const Eigen::Vector3d p1 = pointAt(positions, faces(face, 1));
        const Eigen::Vector3d p2 = pointAt(positions, faces(face, 2));
        const double twiceArea = 2.0 * area[face];
        hatGrad[face][0] = normals[face].cross(p2 - p1) / twiceArea;
        hatGrad[face][1] = normals[face].cross(p0 - p2) / twiceArea;
        hatGrad[face][2] = normals[face].cross(p1 - p0) / twiceArea;
    }

    // 边界顶点：流函数 v 在边界上取 0（Dirichlet），
    // 这样 n×∇v 的法向通量为零，梯度部分与旋度部分才不会互相重叠、能量占比才有意义。
    std::map<std::pair<int, int>, int> edgeFaces;
    for (int face = 0; face < nFaces; ++face) {
        for (int a = 0; a < 3; ++a) {
            const int i = faces(face, a);
            const int j = faces(face, (a + 1) % 3);
            ++edgeFaces[{std::min(i, j), std::max(i, j)}];
        }
    }
    std::vector<char> onBoundary(static_cast<size_t>(nVerts), 0);
    for (const auto& entry : edgeFaces) {
        if (entry.second == 1) {
            onBoundary[entry.first.first] = 1;
            onBoundary[entry.first.second] = 1;
        }
    }

    std::vector<int> dof(static_cast<size_t>(2 * nVerts), -1);
    int size = 0;
    for (int vertex = 0; vertex < nVerts; ++vertex) dof[vertex] = size++;
    for (int vertex = 0; vertex < nVerts; ++vertex) {
        if (!onBoundary[vertex]) dof[static_cast<size_t>(nVerts + vertex)] = size++;
    }

    // 解 min ‖ξ - ∇u - n×∇v‖²，法方程是 [[L, M], [Mᵀ, L]]。
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(nFaces) * 36);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(size);
    for (int face = 0; face < nFaces; ++face) {
        const Eigen::Vector3d normal = normals[face];
        for (int a = 0; a < 3; ++a) {
            const int i = faces(face, a);
            const Eigen::Vector3d gi = hatGrad[face][a];
            const Eigen::Vector3d jgi = normal.cross(gi);
            if (dof[i] >= 0) rhs[dof[i]] += area[face] * gi.dot(field[face]);
            if (dof[nVerts + i] >= 0) rhs[dof[nVerts + i]] += area[face] * jgi.dot(field[face]);
            for (int b = 0; b < 3; ++b) {
                const int j = faces(face, b);
                const Eigen::Vector3d gj = hatGrad[face][b];
                const double shared = area[face] * gi.dot(gj);
                const double mixed = area[face] * gi.dot(normal.cross(gj));
                const int rowU = dof[i];
                const int colU = dof[j];
                const int rowV = dof[nVerts + i];
                const int colV = dof[nVerts + j];
                if (rowU >= 0 && colU >= 0) triplets.emplace_back(rowU, colU, shared);
                if (rowU >= 0 && colV >= 0) triplets.emplace_back(rowU, colV, mixed);
                if (rowV >= 0 && colU >= 0) triplets.emplace_back(rowV, colU, -mixed);  // Mᵀ = -M
                if (rowV >= 0 && colV >= 0) triplets.emplace_back(rowV, colV, shared);
            }
        }
    }

    Eigen::SparseMatrix<double> system(size, size);
    system.setFromTriplets(triplets.begin(), triplets.end());
    system.makeCompressed();

    // 两个子空间在常量场上重叠，方程本身奇异；加一点 Tikhonov 项取零均值解。
    double scale = 0.0;
    for (int k = 0; k < size; ++k) scale = std::max(scale, std::abs(system.coeff(k, k)));
    const double lambda = (scale > 0.0 ? scale : 1.0) * 1e-9;
    for (int k = 0; k < size; ++k) system.coeffRef(k, k) += lambda;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(system);
    if (solver.info() != Eigen::Success) {
        split.note = "线性系统分解失败";
        return split;
    }
    const Eigen::VectorXd solution = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        split.note = "线性系统求解失败";
        return split;
    }

    split.potentialU.resize(static_cast<size_t>(nVerts));
    split.streamV.resize(static_cast<size_t>(nVerts));
    for (int vertex = 0; vertex < nVerts; ++vertex) {
        split.potentialU[vertex] = solution[dof[vertex]];
        const int vIndex = dof[nVerts + vertex];
        split.streamV[vertex] = vIndex >= 0 ? solution[vIndex] : 0.0;  // 边界上 v ≡ 0
    }

    split.gradient.resize(static_cast<size_t>(nFaces));
    split.curl.resize(static_cast<size_t>(nFaces));
    split.harmonic.resize(static_cast<size_t>(nFaces));
    for (int face = 0; face < nFaces; ++face) {
        Eigen::Vector3d gradU = Eigen::Vector3d::Zero();
        Eigen::Vector3d gradV = Eigen::Vector3d::Zero();
        for (int a = 0; a < 3; ++a) {
            const int vertex = faces(face, a);
            gradU += split.potentialU[vertex] * hatGrad[face][a];
            gradV += split.streamV[vertex] * hatGrad[face][a];
        }
        const Eigen::Vector3d curlPart = normals[face].cross(gradV);
        const Eigen::Vector3d harmonic = field[face] - gradU - curlPart;
        split.gradient[face] = glm::vec3(static_cast<float>(gradU.x()), static_cast<float>(gradU.y()),
                                         static_cast<float>(gradU.z()));
        split.curl[face] = glm::vec3(static_cast<float>(curlPart.x()), static_cast<float>(curlPart.y()),
                                     static_cast<float>(curlPart.z()));
        split.harmonic[face] = glm::vec3(static_cast<float>(harmonic.x()), static_cast<float>(harmonic.y()),
                                         static_cast<float>(harmonic.z()));
        split.total += area[face] * field[face].squaredNorm();
        split.gradientEnergy += area[face] * gradU.squaredNorm();
        split.curlEnergy += area[face] * curlPart.squaredNorm();
        split.harmonicEnergy += area[face] * harmonic.squaredNorm();
    }

    // 逐顶点散度 / 旋度，用来检查残差确实两头都不占
    auto divergence = [&](const std::vector<glm::vec3>& data) {
        Eigen::VectorXd out = Eigen::VectorXd::Zero(nVerts);
        for (int face = 0; face < nFaces; ++face) {
            const Eigen::Vector3d value(data[face].x, data[face].y, data[face].z);
            for (int a = 0; a < 3; ++a) out[faces(face, a)] += area[face] * hatGrad[face][a].dot(value);
        }
        return out;
    };
    auto rotation = [&](const std::vector<glm::vec3>& data) {
        Eigen::VectorXd out = Eigen::VectorXd::Zero(nVerts);
        for (int face = 0; face < nFaces; ++face) {
            const Eigen::Vector3d value(data[face].x, data[face].y, data[face].z);
            for (int a = 0; a < 3; ++a) {
                out[faces(face, a)] += area[face] * normals[face].dot(hatGrad[face][a].cross(value));
            }
        }
        return out;
    };

    std::vector<glm::vec3> original(static_cast<size_t>(nFaces));
    for (int face = 0; face < nFaces; ++face) {
        original[face] = glm::vec3(static_cast<float>(field[face].x()), static_cast<float>(field[face].y()),
                                   static_cast<float>(field[face].z()));
    }
    const Eigen::VectorXd divField = divergence(original);
    const Eigen::VectorXd divHarmonic = divergence(split.harmonic);
    const Eigen::VectorXd curlField = rotation(original);
    const Eigen::VectorXd curlHarmonic = rotation(split.harmonic);
    split.divResidual = divField.norm() > 1e-12 ? divHarmonic.norm() / divField.norm() : 0.0;
    split.curlResidual = curlField.norm() > 1e-12 ? curlHarmonic.norm() / curlField.norm() : 0.0;
    split.ok = true;
    split.note =
        "十字场没有连续单值分支：分解前先沿面片取最接近的 90° 分支，奇异点处消不掉的跳变只能落进旋度/调和分量。"
        "梯度占比越高，这个场越接近可积。";
    return split;
}
