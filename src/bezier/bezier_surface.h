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

#include "bezier_patch.h"

//=============================================================================

/// A surface represented by a collection of Bezier patches.
/** This class stores a surface that is represented by a collection of
    bi-cubic tensor-product Bezier patches.
    \sa Bezier_patch
*/
class Bezier_surface
{
public:
    /// constructor, get name of the file to be loaded
    Bezier_surface(const char *_filename = NULL);

    /// destructor
    ~Bezier_surface();

    /// is the object empty, i.e., has nothing been loaded yet?
    bool empty() const { return patches_.empty(); }

    /// load Bezier object from a *.bez file
    bool load_file(const char *_filename);

    /// compute bounding box, store min/max position in _bbmin and _bbmax.
    bool bounding_box(pmp::vec3 &_bbmin, pmp::vec3 &_bbmax) const;

    /// tessellate Bezier surface into triangles with a prescribed resolution.
    void tessellate(unsigned int _resolution);

    /// draw the control polygon for all Bezier patches.
    void draw_control_polygon();

    /// draw the tessellated surface of all Bezier patches.
    void draw_surface(std::string drawmode, bool upload = false);

    /// toggles between de casteljau and bernstein bezier patch evaluation
    void toggle_de_Casteljau();

    /// return current mode (de casteljau or bernstein)
    bool get_bezier_mode();

    /// selects nearest control point to mouse cursor
    void pick(const pmp::vec2 &coord2d, const pmp::mat4 &mvp);

    /// getter for currently selected control point
    pmp::vec3 get_selected_control_point();

    /// setter for currently selected control point
    void set_selected_control_point(const pmp::vec3 &p);

    /// time needed for last tesselation
    float tesselation_time_;

private:
    /// array of all Bezier patches
    std::vector<Bezier_patch> patches_;

    /// currently picked Bezier patch
    Bezier_patch *picked_patch_;
};
//=============================================================================
