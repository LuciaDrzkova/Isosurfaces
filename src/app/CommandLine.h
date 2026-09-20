// Reads the options typed after the program name, e.g.
//   Isosurfaces assets/sphere.obj --function sphere --divisor 1.0

#pragma once

#include "method/Job.h"

#include <iosfwd>

namespace iso {

struct CommandLine
{
    JobSettings job;
    bool window = true;          //!< false with --no-gui
    bool input_given = false;    //!< a mesh file was named
    bool function_given = false; //!< --function was used
};

//! Fills \p cl from the arguments; values not mentioned are left as they were.
//! Returns false if only information was printed (--help, --list-functions,
//! --list-experiments) and the program should exit.
//! \throws std::invalid_argument for wrong usage.
bool parse_command_line(int argc, char** argv, CommandLine& cl);

void print_usage(std::ostream& out);

} // namespace iso
