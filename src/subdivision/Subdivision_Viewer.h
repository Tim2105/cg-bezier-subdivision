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

#include "Mesh.h"
#include <pmp/visualization/MeshViewer.h>

//=============================================================================

/// Simple viewer for a SurfaceMesh

class Subdivision_Viewer : public pmp::MeshViewer
{
public:
    /// constructor
    Subdivision_Viewer(const char* title, int width, int height,
                       bool showgui = true);

    /// load a mesh from file \c filename
    virtual bool load_mesh(const char* filename) override;

    /// draw the scene in different draw modes
    virtual void draw(const std::string& _draw_mode) override;

    /// handle ImGUI interface
    virtual void process_imgui() override;

    /// this function handles keyboard events
    virtual void keyboard(int key, int code, int action, int mod) override;

protected:
    /// the subdivided mesh
    SubdivisionMesh surface_mesh_;

    /// indicates if control mesh should be drawn
    bool draw_control_mesh_;

    /// index of currently loaded mesh
    int mesh_index_;
};
//=============================================================================
