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

#include <pmp/visualization/TrackballViewer.h>
#include <pmp/visualization/Shader.h>

#include "bezier_surface.h"

//=============================================================================

//! Simple viewer for a SurfaceMesh
class BezierViewer : public pmp::TrackballViewer
{
    //! constructor
public:
    /** Constructor, which gets the filename of the Bezier object as its
          parameter. The contructor loads the Bezier object,
          triggers its tessellation, and sets up the model matrix to transform
          the object into the unit sphere. */
    BezierViewer(const char *_title, int _width, int _height,
                 bool showgui = true);

    void loadMesh(const char *filename);

private:
    /// render/handle GUI
    virtual void process_imgui() override;

    /// draw the scene in different draw modes
    virtual void draw(const std::string &_draw_mode) override;

    /// this function handles mouse button events
    virtual void mouse(int button, int action, int mods) override;

    /// this function handles mouse motion (passive/active position)
    virtual void motion(double x, double y) override;

private:
    /// the Bezier object
    Bezier_surface bezier_;

    int tesselation_resolution_;

    /// Phong shader
    pmp::Shader phong_shader_;

    /// should we render the control mesh?
    bool render_control_mesh_;

    int meshIndex_;
};

//=============================================================================
