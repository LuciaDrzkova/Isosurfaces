#pragma once

#include <pmp/visualization/mesh_viewer.h>

namespace iso {

class MyViewer : public pmp::MeshViewer
{
public:
    MyViewer(const char* title, int w, int h) : MeshViewer(title, w, h) {}
    void load_mesh(const char* path) { MeshViewer::load_mesh(path); }
    int run() { return MeshViewer::run(); }
};

} // namespace iso
