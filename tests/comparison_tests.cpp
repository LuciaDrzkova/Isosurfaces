#include "DualContouring.hpp"
#include "ImplicitSurface.hpp"
#include "MarchingCubes.hpp"
#include "MeshGeometryCheck.hpp"
#include "SurfaceFactory.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

namespace {
int failures=0;
#define CHECK(expr) do { if(!(expr)){ std::cerr<<__FILE__<<':'<<__LINE__<<" CHECK failed: "<<#expr<<'\n'; ++failures; return; } } while(false)
void valid_mesh(const iso::mc::Mesh& mesh) {
    for(const auto& t:mesh.triangles) {
        CHECK(t[0]<mesh.vertices.size()); CHECK(t[1]<mesh.vertices.size()); CHECK(t[2]<mesh.vertices.size());
        CHECK(t[0]!=t[1]&&t[1]!=t[2]&&t[2]!=t[0]);
    }
}
void test_surface_factory_and_singularity_values() {
    const auto names=implicit_surface_names();
    CHECK(names.size()==9);
    for(const std::string& name:names) {
        auto surface=make_implicit_surface(name,0.0);
        CHECK(std::isfinite(surface->eval(0.0,0.0,0.0)));
    }
    for(const std::string name:{"cone","elliptic-cone","rotated-cone","x-cone","y-cone","tilted-cone","crossing-planes","quartic"}) {
        auto surface=make_implicit_surface(name,0.0);
        CHECK(std::abs(surface->eval(0.0,0.0,0.0))<1e-12);
    }
}
void test_generic_projection_distance() {
    SphereImplicitSurface sphere(1.0);
    const auto d=iso::geometry_check::distance_to_surface({1.2,0.0,0.0},sphere);
    CHECK(d.converged);
    CHECK(std::abs(d.distance-0.2)<1e-6);
    IntersectingPlanesSurface planes(0.0);
    const auto p=iso::geometry_check::distance_to_surface({0.25,0.0,0.4},planes);
    CHECK(p.converged);
    CHECK(std::abs(p.distance-0.25)<1e-5);
    const auto on_surface=iso::geometry_check::distance_to_surface({1.0,0.0,0.0},sphere);
    CHECK(on_surface.converged&&on_surface.distance==0.0);

    // A tiny |f| is not itself a distance criterion for a quadratic cone:
    // f is quadratic in position near its apex. This point is off the cone and
    // must not silently receive a zero-distance score.
    ParameterizedConeQuadric cone(0.0);
    const auto near_apex=iso::geometry_check::distance_to_surface({1e-6,0.0,0.0},cone);
    CHECK(near_apex.converged);
    CHECK(near_apex.distance>0.0);
}
void test_classic_mc_and_dual_contouring_make_meshes() {
    SphereImplicitSurface sphere(1.0);
    const iso::mc::Bounds bounds{{-1.5,-1.5,-1.5},{1.5,1.5,1.5}};
    iso::mc::Options mo; mo.nx=mo.ny=mo.nz=25; mo.isovalue=0.0;
    const auto mc=iso::mc::extract([&](double x,double y,double z){return sphere.eval(x,y,z);},bounds,mo);
    iso::dc::Options dopt; dopt.nx=dopt.ny=dopt.nz=25; dopt.isovalue=0.0;
    const auto dc=iso::dc::extract([&](double x,double y,double z){return sphere.eval(x,y,z);},bounds,dopt);
    CHECK(!mc.vertices.empty()); CHECK(!mc.triangles.empty());
    CHECK(!dc.vertices.empty()); CHECK(!dc.triangles.empty());
    valid_mesh(mc); valid_mesh(dc);
    const auto report=iso::geometry_check::measure(dc,sphere,{0.0,0.0,0.0},1.0,1);
    CHECK(report.all_vertices.attempted>0);
    CHECK(report.all_vertices.count>0);
    CHECK(report.all_vertices.mean<0.2);
}
void test_surface_family_produces_isosurfaces() {
    const iso::mc::Bounds bounds{{-1.5,-1.5,-1.5},{1.5,1.5,1.5}};
    iso::mc::Options options; options.nx=options.ny=options.nz=21; options.isovalue=0.0;
    for(const std::string name:{"cone","elliptic-cone","rotated-cone","x-cone","y-cone","tilted-cone","crossing-planes","quartic"}) {
        const auto surface=make_implicit_surface(name,0.0);
        const auto mesh=iso::mc::extract([&](double x,double y,double z){return surface->eval(x,y,z);},bounds,options);
        CHECK(!mesh.vertices.empty());
        CHECK(!mesh.triangles.empty());
        valid_mesh(mesh);
    }
}
}
int main() {
    test_surface_factory_and_singularity_values();
    test_generic_projection_distance();
    test_classic_mc_and_dual_contouring_make_meshes();
    test_surface_family_produces_isosurfaces();
    if(failures) {std::cerr<<failures<<" comparison test(s) failed\n";return 1;}
    std::cout<<"All comparison tests passed\n";
    return 0;
}
