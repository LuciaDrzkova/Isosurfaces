#pragma once
#include "MarchingCubes.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

// Uniform-grid Dual Contouring comparison baseline.
// Each active cell gets one vertex by solving a regularized QEF from edge
// intersections and numerically estimated normals. Adjacent active-cell
// vertices are connected around sign-changing grid edges. This is a compact
// reference implementation, not a topology-guaranteed MC33/manifold variant:
// one vertex per cell can merge multiple sheets in ambiguous/singular cells.
namespace iso::dc {
using Point = iso::mc::Point;
using Mesh = iso::mc::Mesh;
using Triangle = iso::mc::Triangle;
using Index = iso::mc::Index;

struct Options {
    std::size_t nx=96, ny=96, nz=96; // sample points; cell count is one less per axis
    double isovalue=0.0;
    double normal_step_scale=0.25;  // fraction of smallest grid spacing
    double qef_regularization=1e-8;
};
namespace detail {
inline constexpr int corner_offset[8][3] = {
    {0,0,0},{1,0,0},{1,1,0},{0,1,0},
    {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
inline constexpr int edge_corners[12][2] = {
    {0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},
    {0,4},{1,5},{2,6},{3,7}};
inline std::size_t grid_index(std::size_t nx,std::size_t ny,
                              std::size_t x,std::size_t y,std::size_t z) {
    return (z*ny+y)*nx+x;
}
inline std::size_t cell_index(std::size_t nx,std::size_t ny,
                              std::size_t x,std::size_t y,std::size_t z) {
    return (z*(ny-1)+y)*(nx-1)+x;
}
inline double norm(const Point& p) {
    return std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}
inline Point add(const Point& a,const Point& b) {
    return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};
}
inline Point sub(const Point& a,const Point& b) {
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}
inline Point mul(const Point& a,double s) { return {a[0]*s,a[1]*s,a[2]*s}; }
inline double dot(const Point& a,const Point& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
inline Point cross(const Point& a,const Point& b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
inline Point lerp(const Point& a,const Point& b,double t) {
    return {a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1]),a[2]+t*(b[2]-a[2])};
}
inline Point clamp_to_cell(Point p,const Point& lo,const Point& hi) {
    for(int i=0;i<3;++i) p[i]=std::clamp(p[i],lo[i],hi[i]);
    return p;
}
inline bool solve_3x3(double A[3][3],double b[3],Point& x) {
    double M[3][4];
    for(int r=0;r<3;++r){for(int c=0;c<3;++c) M[r][c]=A[r][c]; M[r][3]=b[r];}
    for(int c=0;c<3;++c){
        int pivot=c;
        for(int r=c+1;r<3;++r) if(std::abs(M[r][c])>std::abs(M[pivot][c])) pivot=r;
        if(std::abs(M[pivot][c])<1e-20) return false;
        if(pivot!=c) for(int j=c;j<4;++j) std::swap(M[pivot][j],M[c][j]);
        const double d=M[c][c]; for(int j=c;j<4;++j) M[c][j]/=d;
        for(int r=0;r<3;++r) if(r!=c){const double f=M[r][c]; for(int j=c;j<4;++j) M[r][j]-=f*M[c][j];}
    }
    x={M[0][3],M[1][3],M[2][3]};
    return std::isfinite(x[0])&&std::isfinite(x[1])&&std::isfinite(x[2]);
}
template<class ScalarField>
Point gradient(const ScalarField& field,const Point& p,double h) {
    Point g{};
    for(int a=0;a<3;++a){Point lo=p,hi=p;lo[a]-=h;hi[a]+=h;
        g[a]=(field(hi[0],hi[1],hi[2])-field(lo[0],lo[1],lo[2]))/(2.0*h);}
    const double length=norm(g);
    if(length>1e-14&&std::isfinite(length)) return mul(g,1.0/length);
    return {0.0,0.0,0.0};
}
inline bool append_triangle(Mesh& mesh,Index a,Index b,Index c) {
    if(a==b||b==c||c==a) return false;
    const Point ab=sub(mesh.vertices[b],mesh.vertices[a]);
    const Point ac=sub(mesh.vertices[c],mesh.vertices[a]);
    if(norm(cross(ab,ac))<1e-14) return false;
    mesh.triangles.push_back({a,b,c});
    return true;
}
template<class ScalarField>
void emit_quad(Mesh& mesh,std::array<std::int64_t,4> ids,
               const ScalarField& field,double iso,double normal_h) {
    std::vector<Index> unique;
    for(std::int64_t id:ids) if(id>=0){
        const Index v=static_cast<Index>(id);
        if(std::find(unique.begin(),unique.end(),v)==unique.end()) unique.push_back(v);
    }
    if(unique.size()<3) return;
    const auto orient=[&](Index a,Index b,Index c){
        const Point p=mul(add(add(mesh.vertices[a],mesh.vertices[b]),mesh.vertices[c]),1.0/3.0);
        const Point n=cross(sub(mesh.vertices[b],mesh.vertices[a]),sub(mesh.vertices[c],mesh.vertices[a]));
        const Point g=gradient(field,p,normal_h);
        if(dot(n,g)<0.0) std::swap(b,c);
        (void)iso;
        append_triangle(mesh,a,b,c);
    };
    if(unique.size()==3){orient(unique[0],unique[1],unique[2]);return;}
    // Fixed diagonal is deterministic; this is a baseline, not an ambiguity-resolving variant.
    orient(unique[0],unique[1],unique[2]);
    orient(unique[0],unique[2],unique[3]);
}
} // namespace detail

template<class ScalarField>
Mesh extract(ScalarField&& scalar_field,const iso::mc::Bounds& bounds,
             const Options& options={}) {
    const auto& f=scalar_field;
    const std::size_t nx=options.nx,ny=options.ny,nz=options.nz;
    if(nx<2||ny<2||nz<2) throw std::invalid_argument("Dual Contouring requires nx, ny, nz >= 2");
    if(!(bounds.min[0]<bounds.max[0]&&bounds.min[1]<bounds.max[1]&&bounds.min[2]<bounds.max[2]))
        throw std::invalid_argument("Invalid Dual Contouring bounds");
    if(nx>std::numeric_limits<std::size_t>::max()/ny || nx*ny>std::numeric_limits<std::size_t>::max()/nz)
        throw std::invalid_argument("Dual Contouring grid is too large");
    const std::size_t sample_count=nx*ny*nz;
    const std::size_t cell_count=(nx-1)*(ny-1)*(nz-1);
    const double dx=(bounds.max[0]-bounds.min[0])/static_cast<double>(nx-1);
    const double dy=(bounds.max[1]-bounds.min[1])/static_cast<double>(ny-1);
    const double dz=(bounds.max[2]-bounds.min[2])/static_cast<double>(nz-1);
    const double normal_h=std::max(1e-8,options.normal_step_scale*std::min({dx,dy,dz}));
    const auto point_at=[&](std::size_t x,std::size_t y,std::size_t z)->Point{
        return {bounds.min[0]+static_cast<double>(x)*dx,
                bounds.min[1]+static_cast<double>(y)*dy,
                bounds.min[2]+static_cast<double>(z)*dz};};
    std::vector<double> values(sample_count);
    for(std::size_t z=0;z<nz;++z) for(std::size_t y=0;y<ny;++y) for(std::size_t x=0;x<nx;++x){
        const Point p=point_at(x,y,z); values[detail::grid_index(nx,ny,x,y,z)]=f(p[0],p[1],p[2]);
    }
    const auto value=[&](std::size_t x,std::size_t y,std::size_t z){return values[detail::grid_index(nx,ny,x,y,z)];};
    std::vector<std::int64_t> cell_vertices(cell_count,-1);
    Mesh mesh;
    mesh.vertices.reserve(cell_count/4+1);
    // Hermite samples and a regularized QEF solve, one vertex per active cell.
    for(std::size_t z=0;z<nz-1;++z) for(std::size_t y=0;y<ny-1;++y) for(std::size_t x=0;x<nx-1;++x){
        std::array<Point,8> cp{}; std::array<double,8> cv{}; bool below=false,above=false;
        for(int c=0;c<8;++c){const std::size_t xx=x+static_cast<std::size_t>(detail::corner_offset[c][0]);
            const std::size_t yy=y+static_cast<std::size_t>(detail::corner_offset[c][1]);
            const std::size_t zz=z+static_cast<std::size_t>(detail::corner_offset[c][2]);
            cp[c]=point_at(xx,yy,zz); cv[c]=value(xx,yy,zz);
            below|=cv[c]<options.isovalue; above|=cv[c]>=options.isovalue;}
        if(!(below&&above)) continue;
        double A[3][3]{}; double b[3]{}; Point mean{0,0,0}; std::size_t intersections=0,normal_count=0;
        for(int e=0;e<12;++e){const int i=detail::edge_corners[e][0],j=detail::edge_corners[e][1];
            if((cv[i]<options.isovalue)==(cv[j]<options.isovalue)) continue;
            const double denom=cv[j]-cv[i];
            const double t=std::abs(denom)<1e-30?0.5:std::clamp((options.isovalue-cv[i])/denom,0.0,1.0);
            const Point p=detail::lerp(cp[i],cp[j],t); mean=detail::add(mean,p); ++intersections;
            const Point n=detail::gradient(f,p,normal_h); if(detail::norm(n)<1e-12) continue;
            ++normal_count; const double rhs=detail::dot(n,p);
            for(int r=0;r<3;++r){b[r]+=n[r]*rhs;for(int c=0;c<3;++c)A[r][c]+=n[r]*n[c];}
        }
        if(intersections==0) continue;
        mean=detail::mul(mean,1.0/static_cast<double>(intersections));
        Point vertex=mean;
        if(normal_count>0){
            double trace=A[0][0]+A[1][1]+A[2][2];
            const double reg=options.qef_regularization*std::max(1.0,trace);
            for(int d=0;d<3;++d){A[d][d]+=reg;b[d]+=reg*mean[d];}
            Point qef{}; if(detail::solve_3x3(A,b,qef)) vertex=qef;
        }
        const Point lo=point_at(x,y,z), hi=point_at(x+1,y+1,z+1);
        vertex=detail::clamp_to_cell(vertex,lo,hi);
        if(mesh.vertices.size()>=static_cast<std::size_t>(std::numeric_limits<Index>::max()))
            throw std::overflow_error("Dual Contouring vertex index overflow");
        cell_vertices[detail::cell_index(nx,ny,x,y,z)]=static_cast<std::int64_t>(mesh.vertices.size());
        mesh.vertices.push_back(vertex);
    }
    const auto cvtx=[&](std::size_t x,std::size_t y,std::size_t z)->std::int64_t{
        return cell_vertices[detail::cell_index(nx,ny,x,y,z)];};
    // A dual face surrounds every sign-changing primal edge whose four incident
    // cells exist. Domain-boundary edges intentionally remain open.
    for(std::size_t z=1;z+1<nz;++z) for(std::size_t y=1;y+1<ny;++y) for(std::size_t x=0;x+1<nx;++x)
        if((value(x,y,z)<options.isovalue)!=(value(x+1,y,z)<options.isovalue))
            detail::emit_quad(mesh,{cvtx(x,y-1,z-1),cvtx(x,y,z-1),cvtx(x,y,z),cvtx(x,y-1,z)},f,options.isovalue,normal_h);
    for(std::size_t z=1;z+1<nz;++z) for(std::size_t y=0;y+1<ny;++y) for(std::size_t x=1;x+1<nx;++x)
        if((value(x,y,z)<options.isovalue)!=(value(x,y+1,z)<options.isovalue))
            detail::emit_quad(mesh,{cvtx(x-1,y,z-1),cvtx(x-1,y,z),cvtx(x,y,z),cvtx(x,y,z-1)},f,options.isovalue,normal_h);
    for(std::size_t z=0;z+1<nz;++z) for(std::size_t y=1;y+1<ny;++y) for(std::size_t x=1;x+1<nx;++x)
        if((value(x,y,z)<options.isovalue)!=(value(x,y,z+1)<options.isovalue))
            detail::emit_quad(mesh,{cvtx(x-1,y-1,z),cvtx(x,y-1,z),cvtx(x,y,z),cvtx(x-1,y,z)},f,options.isovalue,normal_h);
    // Remove active-cell vertices that were not referenced by a face (usually
    // boundary-only/ambiguous cells). They should not inflate component counts.
    std::vector<char> used(mesh.vertices.size(),0);
    for(const Triangle& t:mesh.triangles) used[t[0]]=used[t[1]]=used[t[2]]=1;
    std::vector<Index> remap(mesh.vertices.size(),std::numeric_limits<Index>::max());
    std::vector<Point> compact; compact.reserve(mesh.vertices.size());
    for(std::size_t i=0;i<mesh.vertices.size();++i) if(used[i]) {
        remap[i]=static_cast<Index>(compact.size()); compact.push_back(mesh.vertices[i]);
    }
    for(Triangle& t:mesh.triangles) {t[0]=remap[t[0]];t[1]=remap[t[1]];t[2]=remap[t[2]];}
    mesh.vertices=std::move(compact);
    return mesh;
}
} // namespace iso::dc
