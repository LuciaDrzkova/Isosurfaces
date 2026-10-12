#pragma once
#include "ImplicitSurface.hpp"
#include "MarchingCubes.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Geometric-error diagnostics. The legacy overload preserves the original
// exact distance-to-double-cone benchmark. The generic overload estimates a
// distance to any implicit zero set by damped normal projection. The latter is
// explicitly an estimate (not a global closest-point guarantee), and reports
// projection failures rather than silently treating them as zero error.
namespace iso::geometry_check {
using Point = iso::mc::Point;
using Mesh = iso::mc::Mesh;
using Triangle = iso::mc::Triangle;

struct Statistics {
    std::size_t count = 0;       // successful finite samples
    std::size_t attempted = 0;
    std::size_t failures = 0;
    double mean = std::numeric_limits<double>::quiet_NaN();
    double rms = std::numeric_limits<double>::quiet_NaN();
    double p95 = std::numeric_limits<double>::quiet_NaN();
    double p99 = std::numeric_limits<double>::quiet_NaN();
    double maximum = std::numeric_limits<double>::quiet_NaN();
};
struct Report {
    Statistics all_vertices;
    Statistics all_triangle_samples;
    Statistics local_vertices;
    Statistics local_triangle_samples;
};
struct ProjectionOptions {
    std::size_t max_iterations = 32;
    std::size_t max_backtracking = 12;
    double function_tolerance = 1e-10;
    double distance_tolerance = 1e-10;
    double gradient_squared_tolerance = 1e-24;
    double finite_difference_scale = 1e-5;
    double max_step_scale = 0.5;
};
struct DistanceResult {
    bool converged = false;
    double distance = std::numeric_limits<double>::quiet_NaN();
    std::size_t iterations = 0;
    double residual = std::numeric_limits<double>::infinity();
};

inline double distance(const Point& a, const Point& b) {
    return std::sqrt((a[0]-b[0])*(a[0]-b[0]) +
                     (a[1]-b[1])*(a[1]-b[1]) +
                     (a[2]-b[2])*(a[2]-b[2]));
}
inline double squared_distance(const Point& a, const Point& b) {
    const double dx=a[0]-b[0], dy=a[1]-b[1], dz=a[2]-b[2];
    return dx*dx+dy*dy+dz*dz;
}
inline double distance_to_double_cone(const Point& p) {
    const double radial=std::hypot(p[0],p[1]);
    return std::abs(radial-std::abs(p[2]))/std::sqrt(2.0);
}

inline Point numeric_gradient(const ImplicitSurface& surface, const Point& p,
                              double scale=1e-5) {
    const double norm=std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    const double h=scale*std::max(1.0,norm);
    Point g{};
    for (std::size_t axis=0; axis<3; ++axis) {
        Point lo=p, hi=p;
        lo[axis]-=h; hi[axis]+=h;
        const double flo=surface.eval(lo[0],lo[1],lo[2]);
        const double fhi=surface.eval(hi[0],hi[1],hi[2]);
        g[axis]=(fhi-flo)/(2.0*h);
    }
    return g;
}
inline DistanceResult distance_to_surface(const Point& source,
                                          const ImplicitSurface& surface,
                                          const ProjectionOptions& options = {}) {
    DistanceResult result;
    Point p=source;
    const auto norm_point=[](const Point& q) {
        return std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
    };
    double f=surface.eval(p[0],p[1],p[2]);
    if (!std::isfinite(f)) return result;
    result.residual=std::abs(f);
    // Exact zeros are safe to report as zero distance. Do not use only an
    // absolute |f| threshold here: for quadratic cones, f is O(distance^2),
    // so nearby off-surface points can have tiny residuals while still having
    // a meaningful geometric distance.
    if (f == 0.0) {
        result.converged=true; result.distance=0.0; return result;
    }
    const double distance_tolerance = options.distance_tolerance *
        std::max(1.0, norm_point(source));
    for (std::size_t iter=0; iter<options.max_iterations; ++iter) {
        result.iterations=iter+1;
        const Point g=numeric_gradient(surface,p,options.finite_difference_scale);
        const double g2=g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
        if (!std::isfinite(g2) || g2<=options.gradient_squared_tolerance)
            return result;
        double step[3]={f*g[0]/g2, f*g[1]/g2, f*g[2]/g2};
        const double step_norm=std::sqrt(step[0]*step[0]+step[1]*step[1]+step[2]*step[2]);
        if (!std::isfinite(step_norm)) return result;
        if (step_norm <= distance_tolerance) {
            result.converged=true;
            result.distance=distance(source,p);
            return result;
        }
        const double limit=options.max_step_scale*std::max(1.0,norm_point(p));
        if (step_norm>limit && step_norm>0.0) {
            const double factor=limit/step_norm;
            for (double& value : step) value*=factor;
        }
        bool accepted=false;
        double alpha=1.0;
        Point candidate{};
        double fc=std::numeric_limits<double>::infinity();
        for (std::size_t bt=0; bt<options.max_backtracking; ++bt) {
            for (std::size_t axis=0; axis<3; ++axis)
                candidate[axis]=p[axis]-alpha*step[axis];
            fc=surface.eval(candidate[0],candidate[1],candidate[2]);
            if (std::isfinite(fc) &&
                (std::abs(fc)<std::abs(f) || std::abs(fc)<=options.function_tolerance)) {
                accepted=true; break;
            }
            alpha*=0.5;
        }
        if (!accepted) return result;
        p=candidate; f=fc; result.residual=std::abs(f);
    }
    return result;
}

inline Statistics summarize(std::vector<double> values, std::size_t attempted=0) {
    Statistics out;
    out.attempted=attempted ? attempted : values.size();
    values.erase(std::remove_if(values.begin(),values.end(),[](double v){return !std::isfinite(v);}),values.end());
    out.count=values.size();
    out.failures=out.attempted>=out.count ? out.attempted-out.count : 0;
    if (values.empty()) return out;
    double sum=0.0, sumsq=0.0;
    for (double v:values) { sum+=v; sumsq+=v*v; }
    out.mean=sum/static_cast<double>(values.size());
    out.rms=std::sqrt(sumsq/static_cast<double>(values.size()));
    out.maximum=*std::max_element(values.begin(),values.end());
    std::sort(values.begin(),values.end());
    const auto percentile=[&values](double q) {
        const auto idx=static_cast<std::size_t>(std::ceil(q*static_cast<double>(values.size())))-1;
        return values[std::min(idx,values.size()-1)];
    };
    out.p95=percentile(0.95); out.p99=percentile(0.99);
    return out;
}

inline Report measure(const Mesh& mesh, const Point& center, double local_radius,
                      std::size_t triangle_subdivisions=4) {
    Report report;
    std::vector<double> all_v, local_v, all_f, local_f;
    const double r2=local_radius*local_radius;
    const auto local=[&](const Point& p){return squared_distance(p,center)<=r2;};
    for (const Point& p:mesh.vertices) {
        const double e=distance_to_double_cone(p); all_v.push_back(e);
        if (local(p)) local_v.push_back(e);
    }
    const std::size_t n=std::max<std::size_t>(1,triangle_subdivisions);
    for (const Triangle& t:mesh.triangles) {
        const Point& a=mesh.vertices[t[0]]; const Point& b=mesh.vertices[t[1]]; const Point& c=mesh.vertices[t[2]];
        for (std::size_t i=0;i<=n;++i) for (std::size_t j=0;j<=n-i;++j) {
            const std::size_t k=n-i-j; const double inv=1.0/static_cast<double>(n);
            Point p{(i*a[0]+j*b[0]+k*c[0])*inv,
                    (i*a[1]+j*b[1]+k*c[1])*inv,
                    (i*a[2]+j*b[2]+k*c[2])*inv};
            const double e=distance_to_double_cone(p); all_f.push_back(e);
            if(local(p)) local_f.push_back(e);
        }
    }
    report.all_vertices=summarize(std::move(all_v));
    report.local_vertices=summarize(std::move(local_v));
    report.all_triangle_samples=summarize(std::move(all_f));
    report.local_triangle_samples=summarize(std::move(local_f));
    return report;
}

inline Report measure(const Mesh& mesh, const ImplicitSurface& surface,
                      const Point& center, double local_radius,
                      std::size_t triangle_samples_per_face=1,
                      const ProjectionOptions& projection={}) {
    Report report;
    std::vector<double> all_v, local_v, all_f, local_f;
    std::size_t attempted_v=0, attempted_lv=0, attempted_f=0, attempted_lf=0;
    const double r2=local_radius*local_radius;
    const auto local=[&](const Point& p){return squared_distance(p,center)<=r2;};
    const auto append=[&](const Point& p,std::vector<double>& all,std::vector<double>& in,
                          std::size_t& nall,std::size_t& nin) {
        ++nall; const bool islocal=local(p); if(islocal) ++nin;
        const auto d=distance_to_surface(p,surface,projection);
        const double error=d.converged?d.distance:std::numeric_limits<double>::quiet_NaN();
        all.push_back(error); if(islocal) in.push_back(error);
    };
    for(const Point& p:mesh.vertices) append(p,all_v,local_v,attempted_v,attempted_lv);
    // Use barycentric samples to capture faceting error that vertex-only tests miss.
    const std::size_t n=std::max<std::size_t>(1,triangle_samples_per_face);
    for(const Triangle& t:mesh.triangles) {
        const Point& a=mesh.vertices[t[0]]; const Point& b=mesh.vertices[t[1]]; const Point& c=mesh.vertices[t[2]];
        Point p{(a[0]+b[0]+c[0])/3.0,(a[1]+b[1]+c[1])/3.0,(a[2]+b[2]+c[2])/3.0};
        append(p,all_f,local_f,attempted_f,attempted_lf);
        if(n>1) {
            const Point mids[3]={{(a[0]+b[0])*0.5,(a[1]+b[1])*0.5,(a[2]+b[2])*0.5},
                                 {(b[0]+c[0])*0.5,(b[1]+c[1])*0.5,(b[2]+c[2])*0.5},
                                 {(c[0]+a[0])*0.5,(c[1]+a[1])*0.5,(c[2]+a[2])*0.5}};
            for(const Point& q:mids) append(q,all_f,local_f,attempted_f,attempted_lf);
        }
    }
    report.all_vertices=summarize(std::move(all_v),attempted_v);
    report.local_vertices=summarize(std::move(local_v),attempted_lv);
    report.all_triangle_samples=summarize(std::move(all_f),attempted_f);
    report.local_triangle_samples=summarize(std::move(local_f),attempted_lf);
    return report;
}

inline void print_statistics(const std::string& label,const Statistics& s) {
    std::cout<<"  "<<label<<" ("<<s.count<<"/"<<s.attempted<<" successful samples; "<<s.failures<<" projection failures)\n"
             <<"    mean distance : "<<s.mean<<'\n'<<"    RMS distance  : "<<s.rms<<'\n'
             <<"    95th percentile: "<<s.p95<<'\n'<<"    99th percentile: "<<s.p99<<'\n'
             <<"    maximum       : "<<s.maximum<<'\n';
}
inline void print_report(const std::string& label,const Report& r) {
    std::cout<<"\nGEOMETRIC ERROR: "<<label<<'\n';
    std::cout<<"  Generic values are damped normal-projection estimates, not guaranteed global nearest-point distances.\n";
    print_statistics("all mesh vertices",r.all_vertices);
    print_statistics("triangle samples",r.all_triangle_samples);
    print_statistics("vertices inside local radius",r.local_vertices);
    print_statistics("triangle samples inside local radius",r.local_triangle_samples);
}
} // namespace iso::geometry_check
