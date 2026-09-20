// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include <pmp/visualization/mesh_viewer.h>
#include <pmp/algorithms/smoothing.h>
#include <pmp/algorithms/remeshing.h>
#include <pmp/algorithms/differential_geometry.h>
#include <pmp/algorithms/normals.h>
#include <pmp/algorithms/utilities.h>
#include <pmp/algorithms/distance_point_triangle.h>
#include <pmp/io/read_obj.h>
#include <pmp/io/read_off.h>
#include <pmp/io/io.h>
#include <pmp/io/write_obj.h>
#include <pmp/io/write_off.h>
#include <pmp/bounding_box.h>
#include <functional>

#include <pmp/mat_vec.h>
#include <pmp/surface_mesh.h> // Include Point Mesh Processing library for SurfaceMesh
//#include <pmp/algorithms/remeshing_2.h>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <string>
#include <limits>
#include <unordered_set> // Include for using unordered_set
#include <Eigen/Dense>
#include <Eigen/Core> // Include Eigen
#include <unsupported/Eigen/AutoDiff>
#include <numeric> // For std::accumulate
#include <chrono>  // For runtime measurement

//Autodiff library
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/dual/dual.hpp>

class MyViewer : public pmp::MeshViewer
{
public:
    //! constructor
    MyViewer(const char* title, int width, int height)
        : MeshViewer(title, width, height)
    {
        set_draw_mode("Smooth Shading");
    }

protected:
    //! this function handles keyboard events
    void keyboard(int key, int code, int action, int mod) override;
};
