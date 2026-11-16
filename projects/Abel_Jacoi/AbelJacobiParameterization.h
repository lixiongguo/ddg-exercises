#ifndef ABEL_JACOBI_PARAMETERIZATION_H
#define ABEL_JACOBI_PARAMETERIZATION_H

#include "CutSeamParameterization.h"
#include "AbelJacobi.h"
#include <Eigen/Dense>
#include <vector>

/**
 * Global conformal parameterization via Abel–Jacobi / period lattice (closed, g >= 1).
 */
class AbelJacobiParameterization : public CutSeamParameterization {
public:
    AbelJacobiParameterization(Mesh& mesh0);

    void parameterize() override;

    void setBaseVertex(int v) { baseVertex_ = v; }
    void setSingularities(const std::vector<AbelJacobi::SingularPoint>& singularities);

    /** Run Poincaré–Hopf + Abel divisor check before parameterization. */
    AbelJacobi::DivisorCheck checkSingularities() const;

    const AbelJacobi& engine() const { return aj_; }
    AbelJacobi& engine() { return aj_; }

    /**
     * Quad mesh of φ_1: real part, plus the 90° rotation of its gradient as the
     * conjugate. `points` holds line endpoints in consecutive pairs.
     * Optional triangle buffers are the same mesh, colored by quad-cell parity.
     */
    bool buildQuadLines(std::vector<Eigen::Vector3d>& points,
                        std::vector<Eigen::Vector3d>* trianglePositions = nullptr,
                        std::vector<Eigen::Vector3d>* triangleColors = nullptr,
                        int* cellCount = nullptr);

private:
    AbelJacobi aj_;
    int baseVertex_ = 0;
    std::vector<AbelJacobi::SingularPoint> singularities_;
};

#endif
