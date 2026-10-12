#pragma once
#include "ImplicitSurface.hpp"
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

inline std::unique_ptr<ImplicitSurface> make_implicit_surface(
    const std::string& name, double c)
{
    if (name == "cone") return std::make_unique<ParameterizedConeQuadric>(c);
    if (name == "elliptic-cone") return std::make_unique<EllipticConeSurface>(c);
    if (name == "rotated-cone") return std::make_unique<RotatedEllipticConeSurface>(c);
    if (name == "x-cone") return std::make_unique<XAxisEllipticConeSurface>(c);
    if (name == "y-cone") return std::make_unique<YAxisEllipticConeSurface>(c);
    if (name == "tilted-cone") return std::make_unique<TiltedEllipticConeSurface>(c);
    if (name == "crossing-planes") return std::make_unique<IntersectingPlanesSurface>(c);
    if (name == "quartic") return std::make_unique<QuarticSaddleSurface>(c, 0.05);
    if (name == "sphere") return std::make_unique<SphereImplicitSurface>(1.0);
    throw std::invalid_argument("Unknown surface: " + name +
        ". Use cone, elliptic-cone, rotated-cone, x-cone, y-cone, tilted-cone, crossing-planes, quartic, sphere, or all.");
}

inline std::vector<std::string> implicit_surface_names()
{
    return {"cone", "elliptic-cone", "rotated-cone", "x-cone", "y-cone",
            "tilted-cone", "crossing-planes", "quartic", "sphere"};
}
