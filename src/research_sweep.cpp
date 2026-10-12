#include "DualContouring.hpp"
#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "LocalVoxelStitcher.hpp"
#include "MarchingCubes.hpp"
#include "MeshGeometryCheck.hpp"
#include "SingularityDetector.hpp"
#include "SurfaceFactory.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Config {
    double c_min=-0.10, c_max=0.10, c_step=0.05, extent=2.0;
    std::vector<std::size_t> resolutions{33,65};
    std::vector<double> radii{0.20,0.30};
    std::vector<std::size_t> levels{1,2};
    std::string surface="all";
    std::filesystem::path output="outputs/parameter_sweep.csv";
};
struct Candidate {
    bool found=false;
    iso::mc::Point point{0.0,0.0,0.0};
    double gradient_norm=std::numeric_limits<double>::infinity();
    SingularityType type=SingularityType::Regular;
};

std::string number(double x) {
    if(!std::isfinite(x)) return "";
    std::ostringstream out; out<<std::setprecision(17)<<x; return out.str();
}
std::string type_name(SingularityType t) {
    if(t==SingularityType::NonDegenerateSingular) return "NonDegenerateSingular";
    if(t==SingularityType::DegenerateSingular) return "DegenerateSingular";
    return "Regular";
}
std::vector<double> c_values(const Config& cfg,const std::string& surface) {
    if(surface=="sphere") return {0.0}; // sphere radius is a separate fixed control parameter
    std::vector<double> values;
    const auto count=static_cast<long long>(std::floor((cfg.c_max-cfg.c_min)/cfg.c_step+1e-9));
    for(long long i=0;i<=count;++i) values.push_back(cfg.c_min+static_cast<double>(i)*cfg.c_step);
    if(values.empty() || values.back()<cfg.c_max-1e-9) values.push_back(cfg.c_max);
    return values;
}
Candidate find_singularity(const iso::mc::Mesh& mesh,const SingularityDetector& detector) {
    Candidate candidate;
    // Probe the analytic critical point before the sampled mesh. This also
    // works at even resolutions, whose lattice may not include the origin.
    const auto origin=detector.classifyPoint(pmp::Point(0.0f,0.0f,0.0f));
    if(origin.type!=SingularityType::Regular && std::abs(origin.function_value)<=1e-12) {
        candidate.found=true; candidate.point={0.0,0.0,0.0};
        candidate.gradient_norm=origin.gradient_norm; candidate.type=origin.type;
        return candidate;
    }
    for(const auto& p:mesh.vertices) {
        const auto r=detector.classifyPoint(pmp::Point(static_cast<float>(p[0]),
                                                       static_cast<float>(p[1]),
                                                       static_cast<float>(p[2])));
        if(r.type==SingularityType::Regular || std::abs(r.function_value)>1e-12) continue;
        if(!candidate.found || r.gradient_norm<candidate.gradient_norm) {
            candidate.found=true; candidate.point=p; candidate.gradient_norm=r.gradient_norm; candidate.type=r.type;
        }
    }
    return candidate;
}
struct MeshStats {
    iso::stitch::TopologyStatistics topology;
    iso::geometry_check::Report geometry;
};
MeshStats measure_mesh(const iso::mc::Mesh& mesh,const ImplicitSurface& surface,
                       const iso::mc::Point& center,double local_radius) {
    MeshStats result;
    result.topology=iso::stitch::compute_topology(mesh);
    result.geometry=iso::geometry_check::measure(mesh,surface,center,local_radius,1);
    return result;
}
void write_header(std::ofstream& csv) {
    csv << "surface,c,resolution,extent,radius,level,singularity_found,singularity_type,method,"
           "vertices,triangles,components,boundary_edges,nonmanifold_edges,"
           "vertex_error_mean,vertex_error_rms,vertex_error_p95,vertex_error_max,"
           "triangle_error_mean,triangle_error_rms,triangle_error_p95,triangle_error_max,"
           "vertex_projection_failures,triangle_projection_failures,"
           "local_vertex_error_mean,local_vertex_error_rms,local_vertex_error_p95,local_vertex_error_max,"
           "local_triangle_error_mean,local_triangle_error_rms,local_triangle_error_p95,local_triangle_error_max,"
           "local_vertex_projection_failures,local_triangle_projection_failures,"
           "selection_changed,sphere_selected_faces,topological_selected_faces,"
           "stitch_success,seam_loops,max_seam_vertex_distance\n";
}
void write_row(std::ofstream& csv,const std::string& name,double c,std::size_t resolution,
    double extent,double radius,std::size_t level,const Candidate& singularity,
    const std::string& method,const iso::mc::Mesh& mesh,const MeshStats* stats=nullptr,
    const LocalUnfolder::LevelStatistics* unfolder_level=nullptr,
    const iso::stitch::StitchResult* stitch=nullptr) {
    iso::stitch::TopologyStatistics topo;
    if(stats) topo=stats->topology; else topo=iso::stitch::compute_topology(mesh);
    csv << name << ',' << number(c) << ',' << resolution << ',' << number(extent) << ','
        << number(radius) << ',' << level << ',' << (singularity.found?1:0) << ','
        << type_name(singularity.type) << ',' << method << ','
        << topo.vertices << ',' << topo.triangles << ',' << topo.connected_components << ','
        << topo.boundary_edges << ',' << topo.nonmanifold_edges << ',';
    if(stats) {
        const auto& g=stats->geometry;
        const auto append=[&csv](const iso::geometry_check::Statistics& x) {
            csv << number(x.mean) << ',' << number(x.rms) << ',' << number(x.p95) << ','
                << number(x.maximum) << ',';
        };
        append(g.all_vertices); append(g.all_triangle_samples);
        csv << g.all_vertices.failures << ',' << g.all_triangle_samples.failures << ',';
        append(g.local_vertices); append(g.local_triangle_samples);
        csv << g.local_vertices.failures << ',' << g.local_triangle_samples.failures << ',';
    } else {
        // Twenty geometry/error-count fields. Keep failed/empty stats blank.
        for(int i=0;i<20;++i) csv << ',';
    }
    if(unfolder_level)
        csv << (unfolder_level->selection_changed?1:0) << ','
            << unfolder_level->selected_faces_sphere << ','
            << unfolder_level->selected_faces_topological << ',';
    else csv << ",,,";
    if(stitch)
        csv << (stitch->stitched?1:0) << ',' << stitch->seam_loops << ','
            << number(stitch->max_seam_vertex_distance) << '\n';
    else csv << ",,\n";
}
std::string require_value(int& i,int argc,char** argv,const std::string& option) {
    if(i+1>=argc) throw std::invalid_argument("Missing value for "+option);
    return argv[++i];
}
template<class T,class Parse> std::vector<T> parse_list(const std::string& s,Parse parse) {
    std::vector<T> out; std::stringstream ss(s); std::string item;
    while(std::getline(ss,item,',')) if(!item.empty()) out.push_back(parse(item));
    return out;
}
}

int main(int argc,char** argv) {
    std::cout << std::unitbuf;
    try {
        Config cfg;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--surface") cfg.surface=require_value(i,argc,argv,arg);
            else if(arg=="--c-min") cfg.c_min=std::stod(require_value(i,argc,argv,arg));
            else if(arg=="--c-max") cfg.c_max=std::stod(require_value(i,argc,argv,arg));
            else if(arg=="--c-step") cfg.c_step=std::stod(require_value(i,argc,argv,arg));
            else if(arg=="--extent") cfg.extent=std::stod(require_value(i,argc,argv,arg));
            else if(arg=="--output") cfg.output=require_value(i,argc,argv,arg);
            else if(arg=="--resolutions") cfg.resolutions=parse_list<std::size_t>(require_value(i,argc,argv,arg),[](const std::string& s){return static_cast<std::size_t>(std::stoul(s));});
            else if(arg=="--radii") cfg.radii=parse_list<double>(require_value(i,argc,argv,arg),[](const std::string& s){return std::stod(s);});
            else if(arg=="--levels") cfg.levels=parse_list<std::size_t>(require_value(i,argc,argv,arg),[](const std::string& s){return static_cast<std::size_t>(std::stoul(s));});
            else if(arg=="--quick") {cfg.c_min=cfg.c_max=0.0;cfg.c_step=1.0;cfg.resolutions={33};cfg.radii={0.25};cfg.levels={1};}
            else if(arg=="--help") {
                std::cout << "IsosurfacesSweep options\n"
                    << "  --surface all|cone|elliptic-cone|rotated-cone|x-cone|y-cone|tilted-cone|crossing-planes|quartic|sphere\n"
                    << "  --c-min V --c-max V --c-step V --extent V\n"
                    << "  --resolutions N,N,... --radii R,R,... --levels L,L,...\n"
                    << "  --output FILE --quick\n"
                    << "Methods: classic_mc, dual_contouring, local_unfolder_secondary, local_voxel_patch, local_voxel_stitched.\n";
                return 0;
            } else throw std::invalid_argument("Unknown argument: "+arg);
        }
        if(!(cfg.c_step>0.0) || !(cfg.extent>0.0) || !std::isfinite(cfg.c_min) ||
           !std::isfinite(cfg.c_max) || cfg.c_min>cfg.c_max)
            throw std::invalid_argument("Invalid c range, c-step, or extent");
        if(cfg.resolutions.empty()||cfg.radii.empty()||cfg.levels.empty())
            throw std::invalid_argument("Resolution/radius/level lists must not be empty");
        for(auto n:cfg.resolutions) if(n<2) throw std::invalid_argument("Resolution values must be >= 2");
        for(auto r:cfg.radii) if(!(r>0.0)||!std::isfinite(r)) throw std::invalid_argument("Radius values must be finite and positive");
        for(auto l:cfg.levels) if(l==0||l>8) throw std::invalid_argument("Levels must be in [1,8]");
        const std::vector<std::string> surfaces=(cfg.surface=="all")?implicit_surface_names():std::vector<std::string>{cfg.surface};
        for(const auto& name:surfaces) (void)make_implicit_surface(name,0.0); // validate names early
        const auto parent=cfg.output.parent_path(); if(!parent.empty()) std::filesystem::create_directories(parent);
        std::ofstream csv(cfg.output); if(!csv) throw std::runtime_error("Could not open output CSV: "+cfg.output.string());
        csv << std::setprecision(17); write_header(csv);
        const iso::mc::Bounds bounds{{-cfg.extent,-cfg.extent,-cfg.extent},{cfg.extent,cfg.extent,cfg.extent}};
        std::size_t groups=0, rows=0;
        for(const auto& name:surfaces) for(double c:c_values(cfg,name)) for(std::size_t resolution:cfg.resolutions) {
            const auto surface=make_implicit_surface(name,c);
            iso::mc::Options mcopt; mcopt.nx=mcopt.ny=mcopt.nz=resolution; mcopt.isovalue=0.0;
            std::cout<<"[Start] surface"<< name << " c=" << c << " resolution=" << resolution << " stage=MarchingCubes" << std::endl;
            const auto mc=iso::mc::extract([&](double x,double y,double z){return surface->eval(x,y,z);},bounds,mcopt);
            iso::dc::Options dcopt; dcopt.nx=dcopt.ny=dcopt.nz=resolution; dcopt.isovalue=0.0;
            std::cout<<"[Start] surface"<< name << " c=" << c << " resolution=" << resolution << " stage=DualContouring" << std::endl;
            const auto dc=iso::dc::extract([&](double x,double y,double z){return surface->eval(x,y,z);},bounds,dcopt);
            const SingularityDetector detector(*surface,1e-2,1e-3);
            const Candidate singularity=find_singularity(mc,detector);
            const iso::mc::Point center=singularity.found?singularity.point:iso::mc::Point{0.0,0.0,0.0};
            if(!singularity.found) {
                // Still include regular controls and non-singular c-levels as global baselines.
                const auto mcs=measure_mesh(mc,*surface,center,0.0);
                const auto dcs=measure_mesh(dc,*surface,center,0.0);
                write_row(csv,name,c,resolution,cfg.extent,0.0,0,singularity,"classic_mc",mc,&mcs); ++rows;
                write_row(csv,name,c,resolution,cfg.extent,0.0,0,singularity,"dual_contouring",dc,&dcs); ++rows;
                std::cout << '[' << name << "] c=" << c << " N=" << resolution
                          << " MC V/T=" << mc.vertices.size() << '/' << mc.triangles.size()
                          << " DC V/T=" << dc.vertices.size() << '/' << dc.triangles.size()
                          << " singularity=none\n";
                ++groups; continue;
            }
            std::cout << '[' << name << "] c=" << c << " N=" << resolution
                      << " MC V/T=" << mc.vertices.size() << '/' << mc.triangles.size()
                      << " DC V/T=" << dc.vertices.size() << '/' << dc.triangles.size()
                      << " singularity=" << type_name(singularity.type) << '\n';
            for(double radius:cfg.radii) for(std::size_t level:cfg.levels) {
                // Re-measure MC/DC for this radius so the CSV contains a valid
                // local-neighborhood baseline for the same (radius, level).
                const auto mcstats=measure_mesh(mc,*surface,center,radius);
                const auto dcstats=measure_mesh(dc,*surface,center,radius);
                write_row(csv,name,c,resolution,cfg.extent,radius,level,singularity,"classic_mc",mc,&mcstats); ++rows;
                write_row(csv,name,c,resolution,cfg.extent,radius,level,singularity,"dual_contouring",dc,&dcstats); ++rows;

                LocalUnfolder::Options uopt; uopt.radius=radius;uopt.levels=level;
                uopt.region_mode=LocalUnfolder::RegionMode::TopologicalBfs;uopt.sampling_bounds=bounds;uopt.base_resolution=resolution;uopt.isovalue=0.0;
                LocalUnfolder unfolder; const auto baseline=unfolder.refine(mc,*surface,center,uopt);
                const auto baseline_stats=measure_mesh(baseline.mesh,*surface,center,radius);
                const auto* level_stats=baseline.levels.empty()?nullptr:&baseline.levels.back();
                write_row(csv,name,c,resolution,cfg.extent,radius,level,singularity,"local_unfolder_secondary",baseline.mesh,&baseline_stats,level_stats); ++rows;

                LocalMarchingCubes::Options popt; popt.radius=radius;popt.level=level;popt.padding_cells=1;
                popt.sampling_bounds=bounds;popt.base_resolution=resolution;popt.isovalue=0.0;
                LocalMarchingCubes local_mc; const auto patch=local_mc.extract(*surface,center,popt);
                const auto patch_stats=measure_mesh(patch.mesh,*surface,center,radius);
                write_row(csv,name,c,resolution,cfg.extent,radius,level,singularity,"local_voxel_patch",patch.mesh,&patch_stats); ++rows;
                iso::stitch::LocalVoxelStitcher::Options sopt;sopt.radius=radius;sopt.use_sphere_seam=true;
                iso::stitch::LocalVoxelStitcher stitcher; const auto stitched=stitcher.stitch(mc,patch.mesh,center,sopt);
                const auto stitched_stats=measure_mesh(stitched.mesh,*surface,center,radius);
                write_row(csv,name,c,resolution,cfg.extent,radius,level,singularity,"local_voxel_stitched",stitched.mesh,&stitched_stats,nullptr,&stitched); ++rows;
                std::cout << "  radius=" << radius << " level=" << level
                          << " patch V/T=" << patch.mesh.vertices.size() << '/' << patch.mesh.triangles.size()
                          << " stitched=" << (stitched.stitched?"YES":"NO")
                          << " seam loops=" << stitched.seam_loops << '\n';
            }
            ++groups;
        }
        csv.close();
        std::cout << "Wrote " << cfg.output << " (" << groups << " surface/parameter/resolution groups, " << rows << " rows).\n"
                  << "Primary contribution: LocalMarchingCubes + LocalVoxelStitcher. LocalUnfolder is a secondary comparison.\n";
        return 0;
    } catch(const std::exception& e) { std::cerr << "IsosurfacesSweep error: " << e.what() << '\n'; return 2; }
}
