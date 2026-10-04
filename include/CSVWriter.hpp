#pragma once

#include "MetricsEngine.hpp"
#include <fstream>
#include <string>
#include <iostream>

class CSVResultWriter {
private:
    std::ofstream file_;

public:
    explicit CSVResultWriter(const std::string& filename) {
        file_.open(filename);
        if (file_.is_open()) {
            file_ << "c_parameter,num_vertices,num_faces,connected_components,euler_characteristic,"
                  << "mean_imp_error,max_imp_error,mean_euc_dist,max_euc_dist,skinny_triangles\n";
        }
    }

    ~CSVResultWriter() {
        if (file_.is_open()) file_.close();
    }

    void writeRow(double c, const MeshMetrics& m) {
        if (!file_.is_open()) return;
        file_ << c << ","
              << m.num_vertices << ","
              << m.num_faces << ","
              << m.num_connected_components << ","
              << m.euler_characteristic << ","
              << m.mean_implicit_error << ","
              << m.max_implicit_error << ","
              << m.mean_euclidean_dist << ","
              << m.max_euclidean_dist << ","
              << m.skinny_triangles_count << "\n";
    }
};
