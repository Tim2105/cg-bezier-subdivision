//=============================================================================
//
//   Exercise code for the lecture "Computer Graphics"
//     by Prof. Mario Botsch, TU Dortmund
//
//   Copyright (C)  Computer Graphics Group, TU Dortmund
//
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/visualization/GL.h>
#include <pmp/visualization/Shader.h>
#include <pmp/MatVec.h>
#include <pmp/visualization/SurfaceMeshGL.h>

//=============================================================================

/// Class for subdividable surface mesh
class SubdivisionMesh : public pmp::SurfaceMeshGL
{
public:
    /// Constructor
    SubdivisionMesh();

    /// subdivides the mesh using the Catmull-Clark scheme
    void subdivide();
};
//=============================================================================
