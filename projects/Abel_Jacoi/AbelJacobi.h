#ifndef ABEL_JACOBI_H
#define ABEL_JACOBI_H

#include "Mesh.h"
#include "TreeCotreeBasis.h"
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <map>
#include <vector>
#include <utility>

/**
 * Discrete Abel–Jacobi map for a closed triangle mesh of genus g >= 1.
 *
 * Gu's construction:
 *   1. Tree–cotree homology basis {a_i, b_i}
 *   2. Real harmonic basis {η_k} dual to that basis, ∫_{γ_i} η_j = δ_ij
 *   3. Hodge star inside the harmonic space, *η_k = Σ λ η
 *   4. Canonical holomorphic basis φ = ψ A^{-1}, ∫_{a_i} φ_j = δ_ij
 *   5. Period matrix Π = B, lattice Γ = Z^g + Π Z^g
 *   6. μ(p) = ∫_{p0}^p φ  mod Γ
 */
class AbelJacobi {
public:
    /** Quad-mesh singularity: degree = 4 - valence (3-valence → +1, 5-valence → -1). */
    struct SingularPoint {
        int vertexIndex = -1;
        int valence = 4;
        int degree = 0;
    };

    struct DivisorCheck {
        bool poincareHopfOk = false;
        bool abelJacobiOk = false;
        /** μ(D) in C^g: row j = (Re, Im) of component j. */
        Eigen::MatrixXd muD;
        /** Nearest integer coefficients on {λ_{a_k}, λ_{b_k}} basis. */
        Eigen::VectorXi latticeCoeffs;
        double latticeResidual = 0.0;
    };

    explicit AbelJacobi(Mesh& mesh);

    /** Extract topology, build canonical basis and period lattice. */
    bool build();

    int genus() const { return genus_; }
    int numEdges() const { return nE_; }
    bool isReady() const { return ready_; }

    /** Poincaré–Hopf: Σ (4-k) n_k = 4χ for quad valences. */
    static bool checkPoincareHopf(int genus, int nBoundary,
                                  const std::vector<SingularPoint>& singularities);

    /** μ(D) = Σ d_i μ(p_i) and membership in Γ. */
    DivisorCheck checkDivisor(int baseVertex,
                              const std::vector<SingularPoint>& singularities) const;

    /** g×2 matrix: row j = (Re μ(p), Im μ(p)), reduced into the fundamental domain of Γ. */
    Eigen::MatrixXd abelJacobiMap(int baseVertex, int targetVertex) const;

    /** Period lattice generators. Columns (2k, 2k+1) are λ_{a_k}; columns (2g+2k, 2g+2k+1) are λ_{b_k}. */
    const Eigen::MatrixXd& latticeGenerators() const { return latticeGens_; }

    /** Riemann matrix Π after a-periods are normalized to δ_ij. Entry (i, j) is ∫_{b_i} φ_j. */
    const Eigen::MatrixXcd& periodMatrix() const { return period_; }

    /** Max |∫_{a_i} φ_j − δ_ij| after canonicalization. */
    double canonicalError() const { return canonicalError_; }

    /** Integrate Σ_j (cRe_j + i cIm_j) φ_j along BFS tree from baseVertex. */
    bool integrateHolomorphic(int baseVertex,
                              const Eigen::VectorXd& coeffRe,
                              const Eigen::VectorXd& coeffIm,
                              std::vector<Eigen::Vector2d>& uv) const;

    /** Integrate the first canonical form φ_1 and write (Re, Im) into mesh UVs. */
    bool quantizeAndParameterize(int baseVertex = 0);

    /** Edges of the a- and b-cycles, as vertex pairs on the original mesh. */
    std::vector<std::pair<int, int>> homologySeams() const;

    /** (Re, Im) of φ_which along the oriented edge from → to. */
    Eigen::Vector2d canonicalEdge(int from, int to, int which = 0) const;

private:
    struct EdgeInfo {
        int v0 = 0;
        int v1 = 0;
        int f0 = -1;
        int f1 = -1;
        double length = 0.0;
    };

    struct HoloForm {
        Eigen::VectorXd omega;
        Eigen::VectorXd starOmega;
    };

    void extractMesh();
    int findEdge(int a, int b) const;
    int edgeSign(int from, int to) const;
    double cotan(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c) const;
    double cotanWeightAtVertex(int vi, int edgeIdx) const;
    void collectHomologyBasis();
    bool solveHarmonic(const std::vector<double>& periodRhs, Eigen::VectorXd& omega) const;
    double periodOf(const Eigen::VectorXd& omega, const topology::TreeCotreeBasis::Cycle& cycle) const;
    double triangleArea(int face) const;
    void faceEdgeValues(int face, const Eigen::VectorXd& omega, Eigen::Vector3d& values) const;
    double faceWedge(const Eigen::Vector3d& left, const Eigen::Vector3d& right) const;
    double faceInner(int face, const Eigen::Vector3d& left, const Eigen::Vector3d& right) const;
    Eigen::Vector2d periodOnCycle(const HoloForm& form,
                                  const topology::TreeCotreeBasis::Cycle& cycle) const;
    void buildLattice();
    Eigen::Vector2d integrateEdge(int from, int to, const HoloForm& form) const;
    Eigen::MatrixXd mapAlongPath(int baseVertex, int targetVertex) const;
    Eigen::MatrixXd reduceModLattice(const Eigen::MatrixXd& mu, Eigen::VectorXi& coeffs, double& residual) const;
    bool nearestLatticePoint(const Eigen::MatrixXd& mu, Eigen::VectorXi& coeffs, double& residual) const;

    Mesh& mesh_;
    Eigen::MatrixXd V_;
    Eigen::VectorXi F_;
    int nV_ = 0;
    int nF_ = 0;
    int nE_ = 0;
    int nB_ = 0;
    int genus_ = 0;
    bool ready_ = false;

    std::vector<EdgeInfo> edges_;
    std::map<std::pair<int, int>, int> edgeMap_;
    std::vector<std::vector<std::pair<int, int>>> vtxInc_;

    std::vector<topology::TreeCotreeBasis::Cycle> aCycles_;
    std::vector<topology::TreeCotreeBasis::Cycle> bCycles_;
    std::vector<HoloForm> canonicalBasis_;
    /** 2g generators, each stored as g×2 (Re, Im per component). */
    Eigen::MatrixXd latticeGens_;
    Eigen::MatrixXcd period_;
    double canonicalError_ = 0.0;
};

#endif
