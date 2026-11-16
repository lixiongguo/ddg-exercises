#include "AbelJacobiParameterization.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <set>

AbelJacobiParameterization::AbelJacobiParameterization(Mesh& mesh0)
    : CutSeamParameterization(mesh0)
    , aj_(mesh)
{
}

void AbelJacobiParameterization::setSingularities(
    const std::vector<AbelJacobi::SingularPoint>& singularities)
{
    singularities_ = singularities;
}

AbelJacobi::DivisorCheck AbelJacobiParameterization::checkSingularities() const
{
    return aj_.checkDivisor(baseVertex_, singularities_);
}

void AbelJacobiParameterization::parameterize()
{
    if (!aj_.build()) {
        std::cerr << "[AbelJacobiParameterization] build failed.\n";
        return;
    }

    if (!singularities_.empty()) {
        const auto check = checkSingularities();
        if (!check.poincareHopfOk) {
            std::cerr << "[AbelJacobiParameterization] Warning: Poincare-Hopf constraint not satisfied.\n";
        }
        if (!check.abelJacobiOk) {
            std::cerr << "[AbelJacobiParameterization] Warning: Abel-Jacobi divisor not principal "
                      << "(residual=" << check.latticeResidual << ").\n";
        }
    }

    if (!aj_.quantizeAndParameterize(baseVertex_)) {
        std::cerr << "[AbelJacobiParameterization] quantization failed.\n";
        return;
    }

    normalize();
}

bool AbelJacobiParameterization::buildQuadLines(std::vector<Eigen::Vector3d>& points,
                                                std::vector<Eigen::Vector3d>* trianglePositions,
                                                std::vector<Eigen::Vector3d>* triangleColors,
                                                int* cellCount)
{
    points.clear();
    if (trianglePositions) trianglePositions->clear();
    if (triangleColors) triangleColors->clear();
    if (cellCount) *cellCount = 0;
    if (!aj_.isReady() && !aj_.build()) return false;
    const std::vector<std::pair<int, int>> seams = aj_.homologySeams();
    if (seams.empty()) return false;
    setSeamVertexPairs(seams);
    if (!prepareCutMesh()) return false;

    int nCut = 0;
    for (VertexCIter vertex = cutMesh_.vertices.begin(); vertex != cutMesh_.vertices.end(); ++vertex) {
        nCut = std::max(nCut, vertex->index + 1);
    }
    if (nCut <= 0) return false;

    // ω is the real canonical form. Its conjugate for the quad grid is the
    // pointwise 90° rotation of grad u, averaged onto edges. The harmonic
    // projection stored in φ is the right object for periods, but on this mesh
    // it is several times weaker than ω and far from orthogonal, so its integer
    // levels never cross into quads.
    struct EdgeValue {
        double omega = 0.0;
        bool hasOmega = false;
        double star = 0.0;
        int starCount = 0;
        double length = 0.0;
    };
    std::map<std::pair<int, int>, EdgeValue> edges;

    for (FaceCIter face = cutMesh_.faces.begin(); face != cutMesh_.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        HalfEdgeCIter he = face->he;
        int index[3];
        int original[3];
        Eigen::Vector3d pos[3];
        bool cutEdge[3];
        bool usable = true;
        for (int corner = 0; corner < 3; ++corner) {
            index[corner] = he->vertex->index;
            original[corner] = cutMesh_.originalVertexIndex(index[corner]);
            pos[corner] = he->vertex->position;
            cutEdge[corner] = cutMesh_.isCutBoundaryEdge(he->edge);
            if (index[corner] < 0 || index[corner] >= nCut || original[corner] < 0) usable = false;
            he = he->next;
        }
        if (!usable) continue;

        const double omega[3] = {
            aj_.canonicalEdge(original[0], original[1]).x(),
            aj_.canonicalEdge(original[1], original[2]).x(),
            aj_.canonicalEdge(original[2], original[0]).x(),
        };
        const Eigen::Vector3d e0 = pos[1] - pos[0];
        const Eigen::Vector3d e1 = pos[2] - pos[0];
        const double e0Len = e0.norm();
        Eigen::Vector3d normal = e0.cross(e1);
        const double area2 = normal.norm();
        if (e0Len < 1e-12 || area2 < 1e-12) continue;
        normal /= area2;
        const Eigen::Vector3d xAxis = e0 / e0Len;
        Eigen::Vector3d yAxis = normal.cross(e0);
        yAxis.normalize();
        const double e1x = e1.dot(xAxis);
        const double e1y = e1.dot(yAxis);
        const double det = e0Len * e1y;
        if (std::abs(det) < 1e-14) continue;
        const double along02 = -omega[2];
        const double gradX = (omega[0] * e1y) / det;
        const double gradY = (e0Len * along02 - e1x * omega[0]) / det;
        const Eigen::Vector3d gradU = xAxis * gradX + yAxis * gradY;
        const Eigen::Vector3d gradV = normal.cross(gradU);
        const double star[3] = {
            gradV.dot(pos[1] - pos[0]),
            gradV.dot(pos[2] - pos[1]),
            gradV.dot(pos[0] - pos[2]),
        };

        for (int corner = 0; corner < 3; ++corner) {
            if (cutEdge[corner]) continue;
            const int from = index[corner];
            const int to = index[(corner + 1) % 3];
            if (from == to) continue;
            const int low = std::min(from, to);
            const int high = std::max(from, to);
            const double sign = (from == low) ? 1.0 : -1.0;
            EdgeValue& edge = edges[{low, high}];
            if (!edge.hasOmega) {
                edge.omega = sign * omega[corner];
                edge.length = (pos[corner] - pos[(corner + 1) % 3]).norm();
                edge.hasOmega = true;
            }
            edge.star += sign * star[corner];
            ++edge.starCount;
        }
    }
    if (edges.empty()) return false;

    struct Link {
        int low;
        int high;
        double du;
        double dv;
    };
    std::vector<Link> links;
    links.reserve(edges.size());
    std::vector<std::vector<int>> neighbors(static_cast<size_t>(nCut));
    double sumAbsU = 0.0;
    double sumAbsV = 0.0;
    for (const auto& item : edges) {
        const EdgeValue& edge = item.second;
        if (!edge.hasOmega || edge.starCount <= 0) continue;
        const double dv = edge.star / static_cast<double>(edge.starCount);
        Link link;
        link.low = item.first.first;
        link.high = item.first.second;
        link.du = edge.omega;
        link.dv = dv;
        neighbors[static_cast<size_t>(link.low)].push_back(link.high);
        neighbors[static_cast<size_t>(link.high)].push_back(link.low);
        sumAbsU += std::abs(link.du);
        sumAbsV += std::abs(link.dv);
        links.push_back(link);
    }
    if (links.empty() || sumAbsU < 1e-14 || sumAbsV < 1e-14) return false;

    std::vector<int> pins;
    std::vector<char> seen(static_cast<size_t>(nCut), 0);
    for (int start = 0; start < nCut; ++start) {
        if (seen[static_cast<size_t>(start)] || neighbors[static_cast<size_t>(start)].empty()) continue;
        pins.push_back(start);
        std::queue<int> queue;
        queue.push(start);
        seen[static_cast<size_t>(start)] = 1;
        while (!queue.empty()) {
            const int current = queue.front();
            queue.pop();
            for (int next : neighbors[static_cast<size_t>(current)]) {
                if (seen[static_cast<size_t>(next)]) continue;
                seen[static_cast<size_t>(next)] = 1;
                queue.push(next);
            }
        }
    }
    if (pins.empty()) return false;

    std::vector<int> compact(static_cast<size_t>(nCut), -1);
    int nUsed = 0;
    for (int index = 0; index < nCut; ++index) {
        if (!neighbors[static_cast<size_t>(index)].empty()) compact[static_cast<size_t>(index)] = nUsed++;
    }
    if (nUsed == 0) return false;

    const int nRow = static_cast<int>(links.size() + pins.size());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(nRow) * 2);
    Eigen::VectorXd rhsU(nRow);
    Eigen::VectorXd rhsV(nRow);
    int row = 0;
    for (const Link& link : links) {
        triplets.emplace_back(row, compact[static_cast<size_t>(link.low)], -1.0);
        triplets.emplace_back(row, compact[static_cast<size_t>(link.high)], 1.0);
        rhsU(row) = link.du;
        rhsV(row) = link.dv;
        ++row;
    }
    for (int pin : pins) {
        triplets.emplace_back(row, compact[static_cast<size_t>(pin)], 1.0);
        rhsU(row) = 0.0;
        rhsV(row) = 0.0;
        ++row;
    }

    Eigen::SparseMatrix<double> system(nRow, nUsed);
    system.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    const Eigen::SparseMatrix<double> normal = system.transpose() * system;
    solver.compute(normal);
    if (solver.info() != Eigen::Success) return false;
    const Eigen::VectorXd solvedU = solver.solve(system.transpose() * rhsU);
    const Eigen::VectorXd solvedV = solver.solve(system.transpose() * rhsV);
    if (solver.info() != Eigen::Success || solvedU.size() != nUsed || solvedV.size() != nUsed) return false;

    // One integer step is about four mesh edges, in each family separately,
    // so the two directions form quads instead of one sparse family of loops.
    constexpr double kEdgesPerCell = 4.0;
    const double meanAbsU = sumAbsU / static_cast<double>(links.size());
    const double meanAbsV = sumAbsV / static_cast<double>(links.size());
    const double scaleU = 1.0 / (kEdgesPerCell * meanAbsU);
    const double scaleV = 1.0 / (kEdgesPerCell * meanAbsV);

    std::vector<Eigen::Vector2d> uv(static_cast<size_t>(nCut), Eigen::Vector2d::Zero());
    for (int index = 0; index < nCut; ++index) {
        const int column = compact[static_cast<size_t>(index)];
        if (column < 0 || !seen[static_cast<size_t>(index)]) continue;
        if (!std::isfinite(solvedU(column)) || !std::isfinite(solvedV(column))) {
            seen[static_cast<size_t>(index)] = 0;
            continue;
        }
        uv[static_cast<size_t>(index)] = Eigen::Vector2d(scaleU * solvedU(column), scaleV * solvedV(column));
    }

    Eigen::Vector3d lo = Eigen::Vector3d::Constant(1e30);
    Eigen::Vector3d hi = Eigen::Vector3d::Constant(-1e30);
    for (VertexCIter vertex = cutMesh_.vertices.begin(); vertex != cutMesh_.vertices.end(); ++vertex) {
        lo = lo.cwiseMin(vertex->position);
        hi = hi.cwiseMax(vertex->position);
    }
    const double diagonal = (hi - lo).norm();
    const double lineLift = std::max(0.004 * diagonal, 1e-6);
    const double faceLift = std::max(0.0015 * diagonal, 1e-6);

    const int cornerOf[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    auto trace = [&](const Eigen::Vector3d pos[3], const Eigen::Vector2d chart[3], const Eigen::Vector3d& normal, int axis) {
        const double value[3] = {chart[0][axis], chart[1][axis], chart[2][axis]};
        const double lowValue = std::min(value[0], std::min(value[1], value[2]));
        const double highValue = std::max(value[0], std::max(value[1], value[2]));
        const int first = static_cast<int>(std::floor(lowValue)) + 1;
        const int last = static_cast<int>(std::ceil(highValue)) - 1;
        if (last < first || last - first > 4) return;
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
                hit[count++] = pos[a] * (1.0 - t) + pos[b] * t;
            }
            if (count == 2 && (hit[1] - hit[0]).norm() > 1e-8) {
                points.push_back(hit[0] + normal * lineLift);
                points.push_back(hit[1] + normal * lineLift);
            }
        }
    };

    std::set<std::pair<int, int>> cells;
    const Eigen::Vector3d cream(0.97, 0.84, 0.58);
    const Eigen::Vector3d orange(0.86, 0.46, 0.18);
    for (FaceCIter face = cutMesh_.faces.begin(); face != cutMesh_.faces.end(); ++face) {
        if (face->isBoundary()) continue;
        HalfEdgeCIter he = face->he;
        VertexCIter verts[3] = {he->vertex, he->next->vertex, he->next->next->vertex};
        Eigen::Vector3d pos[3];
        Eigen::Vector2d chart[3];
        bool ready = true;
        for (int corner = 0; corner < 3; ++corner) {
            const int index = verts[corner]->index;
            if (index < 0 || index >= nCut || !seen[static_cast<size_t>(index)]) {
                ready = false;
                break;
            }
            pos[corner] = verts[corner]->position;
            chart[corner] = uv[static_cast<size_t>(index)];
        }
        if (!ready) continue;
        Eigen::Vector3d normal = (pos[1] - pos[0]).cross(pos[2] - pos[0]);
        const double length = normal.norm();
        if (length < 1e-12) continue;
        normal /= length;
        const Eigen::Vector2d centroid = (chart[0] + chart[1] + chart[2]) / 3.0;
        if (!std::isfinite(centroid.x()) || !std::isfinite(centroid.y())) continue;
        const int cellU = static_cast<int>(std::floor(centroid.x()));
        const int cellV = static_cast<int>(std::floor(centroid.y()));
        cells.emplace(cellU, cellV);
        if (trianglePositions && triangleColors) {
            const Eigen::Vector3d color = ((cellU + cellV) & 1) ? orange : cream;
            for (int corner = 0; corner < 3; ++corner) {
                trianglePositions->push_back(pos[corner] + normal * faceLift);
                triangleColors->push_back(color);
            }
        }
        trace(pos, chart, normal, 0);
        trace(pos, chart, normal, 1);
    }
    if (cellCount) *cellCount = static_cast<int>(cells.size());
    return points.size() >= 8 && cells.size() >= 4;
}
