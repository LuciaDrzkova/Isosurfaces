#pragma once
// Header-only classic Marching Cubes for a scalar field sampled on a regular grid.
//
// The implementation uses the standard Lorensen-Cline / Bourke-style 256-case
// table, linear interpolation along cell edges, and a per-grid-edge cache so
// adjacent cells share the same mesh vertex.
//
// Coordinates returned by extract() are in the caller-provided world-space bounds.
// The scalar field is sampled at nx * ny * nz grid points.
//
// Classic Marching Cubes has ambiguous configurations. This is intentionally the
// classic 256-case algorithm, not MC33/Lewiner.
//
// Dependencies: C++17 standard library only.

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iso::mc {

using Index = std::uint32_t;
using Point = std::array<double, 3>;
using Triangle = std::array<Index, 3>;

struct Bounds {
    Point min;
    Point max;
};

struct Mesh {
    std::vector<Point> vertices;
    std::vector<Triangle> triangles;
};

struct Options {
    std::size_t nx = 96;
    std::size_t ny = 96;
    std::size_t nz = 96;
    double isovalue = 0.0;
};

namespace detail {

// Corner order MUST match the Lorensen-Cline / Paul Bourke lookup table
// stored below. The table convention is:
//
//   7------6
//  /|     /|
// 4------5 |
// | 3----|-2
// |/     |/
// 0------1
//
// with the coordinate assignment:
// 0=(0,0,0)  1=(1,0,0)  2=(1,0,1)  3=(0,0,1)
// 4=(0,1,0)  5=(1,1,0)  6=(1,1,1)  7=(0,1,1)
//
// Edge numbering:
// 0:0-1  1:1-2  2:2-3  3:3-0
// 4:4-5  5:5-6  6:6-7  7:7-4
// 8:0-4  9:1-5 10:2-6 11:3-7
//
// This is effectively the conventional table with the y/z axes named in
// accordance with the cell layout used by the reference implementation.
// The important requirement is that the corner offsets and edge cache below
// use exactly this same numbering.

inline constexpr int corner_offset[8][3] = {
    {0,0,0}, {1,0,0}, {1,0,1}, {0,0,1},
    {0,1,0}, {1,1,0}, {1,1,1}, {0,1,1}
};

inline constexpr int edge_corners[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0},
    {4,5}, {5,6}, {6,7}, {7,4},
    {0,4}, {1,5}, {2,6}, {3,7}
};

// Each case stores the triangle-table edge indices as 16 nibbles.
// 0xF is the end marker. This is the standard 256-case classic MC table
// (the same convention used by the linked nihaljn implementation).
inline constexpr std::uint64_t tri_table[256] = {
    0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFF380ULL, 0xFFFFFFFFFFFFF910ULL, 0xFFFFFFFFFF189381ULL,
    0xFFFFFFFFFFFFFA21ULL, 0xFFFFFFFFFFA21380ULL, 0xFFFFFFFFFF920A29ULL, 0xFFFFFFF89A8A2382ULL,
    0xFFFFFFFFFFFFF2B3ULL, 0xFFFFFFFFFF0B82B0ULL, 0xFFFFFFFFFFB32091ULL, 0xFFFFFFFB89B912B1ULL,
    0xFFFFFFFFFF3AB1A3ULL, 0xFFFFFFFAB8A801A0ULL, 0xFFFFFFF9AB9B3093ULL, 0xFFFFFFFFFFB8AA89ULL,
    0xFFFFFFFFFFFFF874ULL, 0xFFFFFFFFFF437034ULL, 0xFFFFFFFFFF748910ULL, 0xFFFFFFF137174914ULL,
    0xFFFFFFFFFF748A21ULL, 0xFFFFFFFA21403743ULL, 0xFFFFFFF748209A29ULL, 0xFFFF4973727929A2ULL,
    0xFFFFFFFFFF2B3748ULL, 0xFFFFFFF40242B74BULL, 0xFFFFFFFB32748109ULL, 0xFFFF1292B9B49B74ULL,
    0xFFFFFFF487AB31A3ULL, 0xFFFF4B7401B41AB1ULL, 0xFFFF30BAB9B09874ULL, 0xFFFFFFFAB99B4B74ULL,
    0xFFFFFFFFFFFFF459ULL, 0xFFFFFFFFFF380459ULL, 0xFFFFFFFFFF051450ULL, 0xFFFFFFF513538458ULL,
    0xFFFFFFFFFF459A21ULL, 0xFFFFFFF594A21803ULL, 0xFFFFFFF204245A25ULL, 0xFFFF8434535235A2ULL,
    0xFFFFFFFFFFB32459ULL, 0xFFFFFFF594B802B0ULL, 0xFFFFFFFB32510450ULL, 0xFFFF584B82852512ULL,
    0xFFFFFFF45931AB3AULL, 0xFFFFAB81A8180594ULL, 0xFFFF30BAB5B05045ULL, 0xFFFFFFFB8AA85845ULL,
    0xFFFFFFFFFF975879ULL, 0xFFFFFFF375359039ULL, 0xFFFFFFF751710870ULL, 0xFFFFFFFFFF753351ULL,
    0xFFFFFFF21A759879ULL, 0xFFFF37503505921AULL, 0xFFFF25A758528208ULL, 0xFFFFFFF7533525A2ULL,
    0xFFFFFFF2B3987597ULL, 0xFFFFB72029279759ULL, 0xFFFF751871810B32ULL, 0xFFFFFFF51771B12BULL,
    0xFFFFB3A31A758859ULL, 0xF0ABA010B7905075ULL, 0xF07570805A30B0ABULL, 0xFFFFFFFFFF5B75ABULL,
    0xFFFFFFFFFFFFF56AULL, 0xFFFFFFFFFF6A5380ULL, 0xFFFFFFFFFF6A5109ULL, 0xFFFFFFF6A5891381ULL,
    0xFFFFFFFFFF162561ULL, 0xFFFFFFF803621561ULL, 0xFFFFFFF620609569ULL, 0xFFFF823625285895ULL,
    0xFFFFFFFFFF56AB32ULL, 0xFFFFFFF56A02B80BULL, 0xFFFFFFF6A5B32910ULL, 0xFFFFB892B92916A5ULL,
    0xFFFFFFF315356B36ULL, 0xFFFF6B51505B0B80ULL, 0xFFFF9505606306B3ULL, 0xFFFFFFF89BB96956ULL,
    0xFFFFFFFFFF8746A5ULL, 0xFFFFFFFA56374034ULL, 0xFFFFFFF7486A5091ULL, 0xFFFF49737179156AULL,
    0xFFFFFFF874156216ULL, 0xFFFF743403625521ULL, 0xFFFF620560509748ULL, 0xF962695923497937ULL,
    0xFFFFFFF56A4872B3ULL, 0xFFFFB720242746A5ULL, 0xFFFF6A5B32874910ULL, 0xF6A54B7B492B9129ULL,
    0xFFFF6B51535B3748ULL, 0xFB404B7B016B5B15ULL, 0xF74836B630560950ULL, 0xFFFF9B7974B96956ULL,
    0xFFFFFFFFFFA4694AULL, 0xFFFFFFF380A946A4ULL, 0xFFFFFFF04606A10AULL, 0xFFFFA16468618138ULL,
    0xFFFFFFF462421941ULL, 0xFFFF462942921803ULL, 0xFFFFFFFFFF624420ULL, 0xFFFFFFF624428238ULL,
    0xFFFFFFF32B46A94AULL, 0xFFFF6A4A94B82280ULL, 0xFFFFA164606102B3ULL, 0xF1B8B12184A16146ULL,
    0xFFFF36B319639469ULL, 0xF14641916B0181B8ULL, 0xFFFFFFF4600636B3ULL, 0xFFFFFFFFFF86B846ULL,
    0xFFFFFFFA98A876A7ULL, 0xFFFFA76A907A0370ULL, 0xFFFF0818717A176AULL, 0xFFFFFFF37117A76AULL,
    0xFFFF768981861621ULL, 0xF937390976192962ULL, 0xFFFFFFF206607087ULL, 0xFFFFFFFFFF276237ULL,
    0xFFFF76898A86AB32ULL, 0xF7A9A76790B72702ULL, 0xFB32A767A1871081ULL, 0xFFFF17616A71B12BULL,
    0xF63136B619768698ULL, 0xFFFFFFFFFF76B190ULL, 0xFFFF06B0B3607087ULL, 0xFFFFFFFFFFFFF6B7ULL,
    0xFFFFFFFFFFFFFB67ULL, 0xFFFFFFFFFF67B803ULL, 0xFFFFFFFFFF67B910ULL, 0xFFFFFFF67B138918ULL,
    0xFFFFFFFFFF7B621AULL, 0xFFFFFFF7B6803A21ULL, 0xFFFFFFF7B69A2092ULL, 0xFFFF89A38A3A27B6ULL,
    0xFFFFFFFFFF726327ULL, 0xFFFFFFF026067807ULL, 0xFFFFFFF910732672ULL, 0xFFFF678891681261ULL,
    0xFFFFFFF73171A67AULL, 0xFFFF801781A7167AULL, 0xFFFF7A69A0A70730ULL, 0xFFFFFFF9A88A7A67ULL,
    0xFFFFFFFFFF68B486ULL, 0xFFFFFFF640603B63ULL, 0xFFFFFFF109648B68ULL, 0xFFFF63B139369649ULL,
    0xFFFFFFF1A28B6486ULL, 0xFFFF640B60B03A21ULL, 0xFFFF9A2920B648B4ULL, 0xF36463B34923A39AULL,
    0xFFFFFFF264248328ULL, 0xFFFFFFFFFF264240ULL, 0xFFFF834642432091ULL, 0xFFFFFFF642241491ULL,
    0xFFFF1A6648168318ULL, 0xFFFFFFF40660A01AULL, 0xF39A9303A6834364ULL, 0xFFFFFFFFFF4A649AULL,
    0xFFFFFFFFFFB67594ULL, 0xFFFFFFF67B594380ULL, 0xFFFFFFFB67045105ULL, 0xFFFF51345343867BULL,
    0xFFFFFFFB6721A459ULL, 0xFFFF594380A217B6ULL, 0xFFFF204A24A45B67ULL, 0xF67B25A523453843ULL,
    0xFFFFFFF945267327ULL, 0xFFFF786260680459ULL, 0xFFFF045051673263ULL, 0xF851584812786826ULL,
    0xFFFF73167161A459ULL, 0xF459078701671A61ULL, 0xFA737A6A305A4A04ULL, 0xFFFFA84A458A7A67ULL,
    0xFFFFFFF98B9B6596ULL, 0xFFFF590650360B63ULL, 0xFFFFB65510B508B0ULL, 0xFFFFFFF1355363B6ULL,
    0xFFFF65B8B9B59A21ULL, 0xFA21965690B603B0ULL, 0xF52025A50865B58BULL, 0xFFFF35A3A25363B6ULL,
    0xFFFF283265825985ULL, 0xFFFFFFF260069659ULL, 0xF826283865081851ULL, 0xFFFFFFFFFF612651ULL,
    0xF698965683A61631ULL, 0xFFFF06505960A01AULL, 0xFFFFFFFFFFA65830ULL, 0xFFFFFFFFFFFFF65AULL,
    0xFFFFFFFFFFB57A5BULL, 0xFFFFFFF03857BA5BULL, 0xFFFFFFF091BA57B5ULL, 0xFFFF1381897BA57AULL,
    0xFFFFFFF15717B21BULL, 0xFFFFB27571721380ULL, 0xFFFF7B2209729579ULL, 0xF289823295B27257ULL,
    0xFFFFFFF573532A52ULL, 0xFFFF52A578258028ULL, 0xFFFF2A37353A5109ULL, 0xF25752A278129289ULL,
    0xFFFFFFFFFF573531ULL, 0xFFFFFFF571170780ULL, 0xFFFFFFF735539309ULL, 0xFFFFFFFFFF795789ULL,
    0xFFFFFFF8BA8A5485ULL, 0xFFFF03BBA50B5405ULL, 0xFFFF54ABA8A48910ULL, 0xF41314943B54A4BAULL,
    0xFFFF8548B2582152ULL, 0xFB151B2B543B0B40ULL, 0xF58B8545B2950520ULL, 0xFFFFFFFFFF3B2549ULL,
    0xFFFF483543253A52ULL, 0xFFFFFFF0244252A5ULL, 0xF910854583A532A3ULL, 0xFFFF2492914252A5ULL,
    0xFFFFFFF153358548ULL, 0xFFFFFFFFFF501540ULL, 0xFFFF530509358548ULL, 0xFFFFFFFFFFFFF549ULL,
    0xFFFFFFFBA9B947B4ULL, 0xFFFFBA97B9794380ULL, 0xFFFFB470414B1BA1ULL, 0xF4BAB474A1843413ULL,
    0xFFFF219B294B97B4ULL, 0xF3801B2B197B9479ULL, 0xFFFFFFF04224B47BULL, 0xFFFF42343824B47BULL,
    0xFFFF947732972A92ULL, 0xF70207872A4797A9ULL, 0xFA040A1A472A3A73ULL, 0xFFFFFFFFFF4782A1ULL,
    0xFFFFFFF317714194ULL, 0xFFFF178180714194ULL, 0xFFFFFFFFFF347304ULL, 0xFFFFFFFFFFFFF784ULL,
    0xFFFFFFFFFF8BA8A9ULL, 0xFFFFFFFA9BB93903ULL, 0xFFFFFFFBA88A0A10ULL, 0xFFFFFFFFFFA3BA13ULL,
    0xFFFFFFF8B99B1B21ULL, 0xFFFF9B2921B93903ULL, 0xFFFFFFFFFFB08B20ULL, 0xFFFFFFFFFFFFFB23ULL,
    0xFFFFFFF98AA82832ULL, 0xFFFFFFFFFF2902A9ULL, 0xFFFF8A1810A82832ULL, 0xFFFFFFFFFFFFF2A1ULL,
    0xFFFFFFFFFF819831ULL, 0xFFFFFFFFFFFFF190ULL, 0xFFFFFFFFFFFFF830ULL, 0xFFFFFFFFFFFFFFFFULL,
};

struct EdgeCache {
    std::vector<std::int32_t> x;
    std::vector<std::int32_t> y;
    std::vector<std::int32_t> z;

    EdgeCache(std::size_t nx, std::size_t ny, std::size_t nz)
        : x((nx - 1) * ny * nz, -1),
          y(nx * (ny - 1) * nz, -1),
          z(nx * ny * (nz - 1), -1) {}

    static std::size_t x_index(std::size_t nx, std::size_t ny,
                                std::size_t x, std::size_t y, std::size_t z) {
        return (z * ny + y) * (nx - 1) + x;
    }

    static std::size_t y_index(std::size_t nx, std::size_t ny,
                                std::size_t x, std::size_t y, std::size_t z) {
        return (z * (ny - 1) + y) * nx + x;
    }

    static std::size_t z_index(std::size_t nx, std::size_t ny,
                                std::size_t x, std::size_t y, std::size_t z) {
        return (z * ny + y) * nx + x;
    }
};

inline std::size_t grid_index(std::size_t nx, std::size_t ny,
                              std::size_t x, std::size_t y, std::size_t z) {
    return (z * ny + y) * nx + x;
}

inline Point lerp(const Point& a, const Point& b, double t) {
    return {
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2])
    };
}

// Returns the grid-edge origin and axis for the requested MC edge.
// axis 0 = x edge, 1 = y edge, 2 = z edge.
inline std::pair<std::array<std::size_t,3>, int>
edge_origin_and_axis(std::size_t x, std::size_t y, std::size_t z, int edge) {
    switch (edge) {
        case 0:  return {{x,   y,   z  }, 0};
        case 1:  return {{x+1, y,   z  }, 2};
        case 2:  return {{x,   y,   z+1}, 0};
        case 3:  return {{x,   y,   z  }, 2};
        case 4:  return {{x,   y+1, z  }, 0};
        case 5:  return {{x+1, y+1, z  }, 2};
        case 6:  return {{x,   y+1, z+1}, 0};
        case 7:  return {{x,   y+1, z  }, 2};
        case 8:  return {{x,   y,   z  }, 1};
        case 9:  return {{x+1, y,   z  }, 1};
        case 10: return {{x+1, y,   z+1}, 1};
        case 11: return {{x,   y,   z+1}, 1};
        default: throw std::out_of_range("Marching Cubes edge must be 0..11");
    }
}

} // namespace detail

// Sample and polygonize a scalar field f(x,y,z).
//
// f must be callable as:
//     double value = f(x, y, z);
//
// The grid dimensions are numbers of SAMPLE POINTS, not numbers of cubes.
// Therefore an N x N x N grid contains (N-1)^3 cubes.
template<class ScalarField>
Mesh extract(ScalarField&& f, const Bounds& bounds, const Options& options = {}) {
    const std::size_t nx = options.nx;
    const std::size_t ny = options.ny;
    const std::size_t nz = options.nz;

    if (nx < 2 || ny < 2 || nz < 2)
        throw std::invalid_argument("Marching Cubes requires nx, ny, nz >= 2");

    if (!(bounds.min[0] < bounds.max[0] &&
          bounds.min[1] < bounds.max[1] &&
          bounds.min[2] < bounds.max[2])) {
        throw std::invalid_argument("Invalid Marching Cubes bounds");
    }

    const std::size_t sample_count = nx * ny * nz;
    std::vector<double> values(sample_count);

    const double dx = (bounds.max[0] - bounds.min[0]) / static_cast<double>(nx - 1);
    const double dy = (bounds.max[1] - bounds.min[1]) / static_cast<double>(ny - 1);
    const double dz = (bounds.max[2] - bounds.min[2]) / static_cast<double>(nz - 1);

    for (std::size_t z = 0; z < nz; ++z) {
        const double pz = bounds.min[2] + static_cast<double>(z) * dz;
        for (std::size_t y = 0; y < ny; ++y) {
            const double py = bounds.min[1] + static_cast<double>(y) * dy;
            for (std::size_t x = 0; x < nx; ++x) {
                const double px = bounds.min[0] + static_cast<double>(x) * dx;
                values[detail::grid_index(nx, ny, x, y, z)] = f(px, py, pz);
            }
        }
    }

    Mesh mesh;
    detail::EdgeCache cache(nx, ny, nz);

    // A practical upper-bound estimate for reserve; most cells contain no
    // surface, so this is intentionally conservative.
    mesh.vertices.reserve(std::min<std::size_t>(sample_count, sample_count / 2 + 1024));
    mesh.triangles.reserve(std::min<std::size_t>((nx - 1) * (ny - 1) * (nz - 1) * 2, sample_count * 2));

    auto point_at = [&](std::size_t x, std::size_t y, std::size_t z) -> Point {
        return {
            bounds.min[0] + static_cast<double>(x) * dx,
            bounds.min[1] + static_cast<double>(y) * dy,
            bounds.min[2] + static_cast<double>(z) * dz
        };
    };

    auto get_cached_vertex = [&](std::size_t cx, std::size_t cy, std::size_t cz,
                                 int edge) -> Index {
        const auto [origin, axis] = detail::edge_origin_and_axis(cx, cy, cz, edge);
        const auto ox = origin[0], oy = origin[1], oz = origin[2];

        std::int32_t* slot = nullptr;
        if (axis == 0)
            slot = &cache.x[detail::EdgeCache::x_index(nx, ny, ox, oy, oz)];
        else if (axis == 1)
            slot = &cache.y[detail::EdgeCache::y_index(nx, ny, ox, oy, oz)];
        else
            slot = &cache.z[detail::EdgeCache::z_index(nx, ny, ox, oy, oz)];

        if (*slot >= 0)
            return static_cast<Index>(*slot);

        const int c0 = detail::edge_corners[edge][0];
        const int c1 = detail::edge_corners[edge][1];

        const std::size_t x0 = cx + static_cast<std::size_t>(detail::corner_offset[c0][0]);
        const std::size_t y0 = cy + static_cast<std::size_t>(detail::corner_offset[c0][1]);
        const std::size_t z0 = cz + static_cast<std::size_t>(detail::corner_offset[c0][2]);
        const std::size_t x1 = cx + static_cast<std::size_t>(detail::corner_offset[c1][0]);
        const std::size_t y1 = cy + static_cast<std::size_t>(detail::corner_offset[c1][1]);
        const std::size_t z1 = cz + static_cast<std::size_t>(detail::corner_offset[c1][2]);

        const double v0 = values[detail::grid_index(nx, ny, x0, y0, z0)];
        const double v1 = values[detail::grid_index(nx, ny, x1, y1, z1)];

        double t;
        const double denom = v1 - v0;
        if (std::abs(denom) <= std::numeric_limits<double>::epsilon()) {
            t = 0.5;
        } else {
            t = (options.isovalue - v0) / denom;
        }
        // Roundoff can produce tiny excursions outside [0,1].
        t = std::clamp(t, 0.0, 1.0);

        const Point p0 = point_at(x0, y0, z0);
        const Point p1 = point_at(x1, y1, z1);

        const Index index = static_cast<Index>(mesh.vertices.size());
        mesh.vertices.push_back(detail::lerp(p0, p1, t));
        *slot = static_cast<std::int32_t>(index);
        return index;
    };

    for (std::size_t z = 0; z + 1 < nz; ++z) {
        for (std::size_t y = 0; y + 1 < ny; ++y) {
            for (std::size_t x = 0; x + 1 < nx; ++x) {
                double v[8];
                for (int c = 0; c < 8; ++c) {
                    const std::size_t gx =
                        x + static_cast<std::size_t>(detail::corner_offset[c][0]);
                    const std::size_t gy =
                        y + static_cast<std::size_t>(detail::corner_offset[c][1]);
                    const std::size_t gz =
                        z + static_cast<std::size_t>(detail::corner_offset[c][2]);
                    v[c] = values[detail::grid_index(nx, ny, gx, gy, gz)];
                }

                int cube_index = 0;
                for (int c = 0; c < 8; ++c) {
                    if (v[c] < options.isovalue)
                        cube_index |= (1 << c);
                }

                const std::uint64_t config = detail::tri_table[cube_index];
                if (config == 0xFFFFFFFFFFFFFFFFULL)
                    continue;

                Index edge_vertices[12]{};
                bool used[12]{};

                for (int i = 0; i < 16; ++i) {
                    const int edge = static_cast<int>((config >> (4 * i)) & 0xFULL);
                    if (edge == 0xF)
                        break;
                    if (!used[edge]) {
                        edge_vertices[edge] = get_cached_vertex(x, y, z, edge);
                        used[edge] = true;
                    }
                }

                for (int i = 0; i < 16; i += 3) {
                    const int e0 = static_cast<int>((config >> (4 * i)) & 0xFULL);
                    if (e0 == 0xF)
                        break;
                    const int e1 = static_cast<int>((config >> (4 * (i + 1))) & 0xFULL);
                    const int e2 = static_cast<int>((config >> (4 * (i + 2))) & 0xFULL);

                    const Index a = edge_vertices[e0];
                    const Index b = edge_vertices[e1];
                    const Index c = edge_vertices[e2];

                    // Degenerate triangles can occur in numerically pathological
                    // cells. They are useless to downstream mesh processing.
                    if (a == b || b == c || c == a)
                        continue;

                    mesh.triangles.push_back({a, b, c});
                }
            }
        }
    }

    return mesh;
}

} // namespace iso::mc
