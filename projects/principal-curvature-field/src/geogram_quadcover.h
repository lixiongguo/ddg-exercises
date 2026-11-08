#pragma once

#include <string>
#include <vector>

/**
 * Quadrangulation parameters from Geogram's QuadCover / PGP solver
 * (GEO::GlobalParam2d::frame_field + GEO::GlobalParam2d::PGP).
 *
 * cornerUV holds three (u, v) pairs per triangle, in face order.
 * direction and side are one tangent vector per triangle.
 */
struct GeogramQuadCoverResult {
    bool ok = false;
    std::string error;
    std::vector<double> cornerUV;
    std::vector<double> direction;
    std::vector<double> side;
};

/**
 * @param scaling  Quad size in units of the average edge length. 1 keeps
 *                 about one quad per input edge.
 */
GeogramQuadCoverResult runGeogramQuadCover(
    const double* positions, int nVerts,
    const int* faces, int nFaces,
    double scaling = 1.0);
