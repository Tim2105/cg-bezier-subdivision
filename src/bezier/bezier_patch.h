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

#include <pmp/MatVec.h>
#include <pmp/visualization/GL.h>

#include <vector>

//=============================================================================

/// Bi-cubic Bezier patch.
/** This class represents a bicubic tensor-product Bezier patch with a
    control polygon of 4x4 control points. The class Bezier_surface
    stores a set of Bezier patches to represent more complex surfaces.
    \sa Bezier_surface
*/
class Bezier_patch
{
public:
    /// default constructor
    Bezier_patch();

    /// destructor
    ~Bezier_patch();

    /// compute bounding box, return min/max position in _bbmin and _bbmax.
    void bounding_box(pmp::vec3 &_bbmin, pmp::vec3 &_bbmax) const;

    /// tessellate Bezier patch into a triangle mesh with a prescribed
    /// resolution.
    void tessellate(unsigned int _resolution);

    /// render the control polygon
    void draw_control_polygon();

    /// render the tessellated surface
    void draw_surface(std::string drawmode, bool upload = false);

    /// pick control point from mouse input
    float pick(const pmp::vec2 &coord2d, const pmp::mat4 &mvp);

    /// getter for currently selected control point
    pmp::vec3 get_selected_control_point();

    /// setter for currently selected control point
    void set_selected_control_point(const pmp::vec3 &p);

private:
    /// compute position `_p` and normal `_n` of Bezier patch at parameter (_u,_v)
    void position_normal(float _u, float _v, pmp::vec3 &_p,
                         pmp::vec3 &_n) const;

    /// upload data to OpenGL buffers for control polgyon and tessellated mesh
    void upload_opengl_buffers();

    /// toggle bezier evaluation via de Casteljau or Bernstein polynomials
    void toggle_de_Casteljau();

    /// get current evaluation status
    bool de_Casteljau_used() { return use_de_Casteljau_; }

private:
    // Bezier_surface has to access control_points_ during file load
    friend class Bezier_surface;

    /// control polygon of 4x4 points
    pmp::vec3 control_points_[4][4];
    /// array of edge indices for the control polygon (two indices per edge)
    std::vector<GLuint> control_edges_;

    /// selected control point (used for mouse input)
    int selected_control_point_;

    /// array of vertex positions for the tessellated surface
    std::vector<pmp::vec3> surface_vertices_;
    /// array of normal vectors for the tessellated surface
    std::vector<pmp::vec3> surface_normals_;
    /// array of triangles for tessellated surface (three indices per triangle)
    std::vector<GLuint> surface_triangles_;

    /// OpenGL vertex array object for control polygon
    GLuint cpoly_vertex_array_;
    /// OpenGL buffer object for control points
    GLuint cpoly_vertex_buffer_;
    /// OpenGL buffer object for edge indices of control polygon
    GLuint cpoly_index_buffer_;

    /// OpenGL vertex array object for surface mesh
    GLuint surf_vertex_array_;
    /// OpenGL buffer object for surface vertices
    GLuint surf_vertex_buffer_;
    /// OpenGL buffer object for surface normals
    GLuint surf_normal_buffer_;
    /// OpenGL buffer object for surface triangle indices
    GLuint surf_index_buffer_;
    /// boolean to determine, if evaluation is done via de Casteljau or
    /// Bernstein polynomials
    bool use_de_Casteljau_;
};

//=============================================================================
