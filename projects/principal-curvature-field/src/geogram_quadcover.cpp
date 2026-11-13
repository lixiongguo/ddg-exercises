#include "geogram_quadcover.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/parameterization/mesh_PGP_2d.h>

namespace {

void ensureGeogram() {
    static const bool ready = [] {
        GEO::initialize(GEO::GEOGRAM_INSTALL_LOCALE | GEO::GEOGRAM_INSTALL_ERRNO);
        return true;
    }();
    (void)ready;
}

GEO::vec3 triangleNormal(const GEO::Mesh& mesh, GEO::index_t face) {
    const GEO::vec3 a = mesh.vertices.point<3>(mesh.facets.vertex(face, 0));
    const GEO::vec3 b = mesh.vertices.point<3>(mesh.facets.vertex(face, 1));
    const GEO::vec3 c = mesh.vertices.point<3>(mesh.facets.vertex(face, 2));
    return GEO::cross(b - a, c - a);
}

} // namespace

GeogramQuadCoverResult runGeogramQuadCover(
    const double* positions, int nVerts,
    const int* faces, int nFaces,
    double scaling) {
    GeogramQuadCoverResult result;
    if (positions == nullptr || faces == nullptr || nVerts < 3 || nFaces < 1 || !(scaling > 0.0)) {
        result.error = "网格为空";
        return result;
    }
    for (int face = 0; face < nFaces; ++face) {
        for (int corner = 0; corner < 3; ++corner) {
            const int vertex = faces[face * 3 + corner];
            if (vertex < 0 || vertex >= nVerts) {
                result.error = "三角形索引越界";
                return result;
            }
        }
    }

    try {
        ensureGeogram();
        GEO::Mesh mesh;
        mesh.vertices.set_dimension(3);
        mesh.vertices.create_vertices(static_cast<GEO::index_t>(nVerts));
        for (int vertex = 0; vertex < nVerts; ++vertex) {
            double* point = mesh.vertices.point_ptr(static_cast<GEO::index_t>(vertex));
            point[0] = positions[vertex * 3];
            point[1] = positions[vertex * 3 + 1];
            point[2] = positions[vertex * 3 + 2];
        }
        mesh.facets.create_triangles(static_cast<GEO::index_t>(nFaces));
        for (int face = 0; face < nFaces; ++face) {
            for (int corner = 0; corner < 3; ++corner) {
                mesh.facets.set_vertex(
                    static_cast<GEO::index_t>(face),
                    static_cast<GEO::index_t>(corner),
                    static_cast<GEO::index_t>(faces[face * 3 + corner]));
            }
        }
        mesh.facets.connect();

        GEO::Attribute<GEO::vec3> frame(mesh.facets.attributes(), "B");
        GEO::GlobalParam2d::frame_field(&mesh, frame, 45.0);
        GEO::Attribute<GEO::vec2> param(mesh.facet_corners.attributes(), "U");
        // The bundled OpenNL has no usable CHOLMOD, so use Jacobi PCG.
        GEO::GlobalParam2d::PGP(&mesh, frame, param, scaling, false, false, 2.0);

        result.cornerUV.assign(static_cast<size_t>(nFaces) * 6, 0.0);
        result.direction.assign(static_cast<size_t>(nFaces) * 3, 0.0);
        result.side.assign(static_cast<size_t>(nFaces) * 3, 0.0);
        int finiteCorners = 0;
        for (int face = 0; face < nFaces; ++face) {
            const GEO::index_t index = static_cast<GEO::index_t>(face);
            GEO::vec3 dir = frame[index];
            GEO::vec3 normal = triangleNormal(mesh, index);
            const double normalLength = GEO::length(normal);
            if (normalLength > 1e-14) normal = normal / normalLength;
            dir = dir - normal * GEO::dot(dir, normal);
            const double dirLength = GEO::length(dir);
            GEO::vec3 side(0.0, 0.0, 0.0);
            if (dirLength > 1e-14) {
                dir = dir / dirLength;
                side = GEO::cross(normal, dir);
            }
            result.direction[static_cast<size_t>(face) * 3] = dir.x;
            result.direction[static_cast<size_t>(face) * 3 + 1] = dir.y;
            result.direction[static_cast<size_t>(face) * 3 + 2] = dir.z;
            result.side[static_cast<size_t>(face) * 3] = side.x;
            result.side[static_cast<size_t>(face) * 3 + 1] = side.y;
            result.side[static_cast<size_t>(face) * 3 + 2] = side.z;
            for (int corner = 0; corner < 3; ++corner) {
                const GEO::index_t cornerIndex = mesh.facets.corner(index, static_cast<GEO::index_t>(corner));
                const double u = param[cornerIndex].x;
                const double v = param[cornerIndex].y;
                result.cornerUV[static_cast<size_t>(face) * 6 + static_cast<size_t>(corner) * 2] = u;
                result.cornerUV[static_cast<size_t>(face) * 6 + static_cast<size_t>(corner) * 2 + 1] = v;
                if (std::isfinite(u) && std::isfinite(v)) ++finiteCorners;
            }
        }
        if (finiteCorners == 0) {
            result.error = "没有得到参数坐标";
            return result;
        }
        result.ok = true;
    } catch (const std::exception& exception) {
        result.error = exception.what();
    }
    return result;
}
