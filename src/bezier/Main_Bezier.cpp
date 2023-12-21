//=============================================================================
//
//   Exercise code for the lecture "Computer Graphics"
//     by Prof. Mario Botsch, TU Dortmund
//
//   Copyright (C)  Computer Graphics Group, TU Dortmund
//
//=============================================================================

#include "BezierViewer.h"

//=============================================================================

int main(int argc, char **argv)
{
    BezierViewer viewer("BezierViewer", 800, 600);
    if (argc > 1)
    {
        viewer.loadMesh(argv[1]);
    }
    return viewer.run();
}

//=============================================================================
