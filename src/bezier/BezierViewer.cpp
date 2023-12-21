//=============================================================================
//
//   Exercise code for the lecture "Computer Graphics"
//     by Prof. Mario Botsch, TU Dortmund
//
//   Copyright (C)  Computer Graphics Group, TU Dortmund
//
//=============================================================================

#include "BezierViewer.h"

#include <cfloat>
#include <imgui.h>
#include <iostream>
#include <sstream>

#include <pmp/visualization/PhongShader.h>

//=============================================================================

using namespace pmp;

BezierViewer::BezierViewer(const char *title, int width, int height,
                           bool showgui)
    : TrackballViewer(title, width, height, showgui),
      tesselation_resolution_(20)
{
    render_control_mesh_ = true;
    meshIndex_ = 0;

    clear_draw_modes();
    add_draw_mode("Hidden Line");
    add_draw_mode("Smooth Shading");
    add_draw_mode("Points");
    set_draw_mode("Hidden Line");

    // add imgui help items
    add_help_item("CTRL + LMB", "Move Control Points", 0);
}

//-----------------------------------------------------------------------------

void BezierViewer::loadMesh(const char *filename)
{
    bezier_.load_file(filename);

    // store model dimensions
    vec3 bbmin, bbmax;
    bezier_.bounding_box(bbmin, bbmax);
    vec3 center = 0.5 * (bbmin + bbmax);
    float radius = 0.5 * distance(bbmin, bbmax);
    set_scene(center, radius);
    // tessellate Bezier surface with resolution
    bezier_.tessellate(tesselation_resolution_);
}

//-----------------------------------------------------------------------------

void BezierViewer::process_imgui()
{
    if (ImGui::CollapsingHeader("Load Mesh", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::PushItemWidth(100);
        const char *listbox_items[] = {"simple.bez", "heart.bez", "teacup.bez",
                                       "teapot.bez", "car.bez"};
        int listbox_item_current = meshIndex_;
        ImGui::ListBox(" ", &listbox_item_current, listbox_items,
                       IM_ARRAYSIZE(listbox_items), 5);
        if (listbox_item_current != meshIndex_ || bezier_.empty())
        {
            std::stringstream ss;
            ss << DATA_PATH << listbox_items[listbox_item_current];

            loadMesh(ss.str().c_str());
            meshIndex_ = listbox_item_current;
        }
        ImGui::PopItemWidth();
        ImGui::Spacing();
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Drawmode", ImGuiTreeNodeFlags_DefaultOpen))
    {
        int mode = draw_mode_;
        ImGui::PushItemWidth(100);
        for (int i = 0; i < draw_mode_names_.size(); ++i)
        {
            ImGui::RadioButton(draw_mode_names_[i].c_str(), &mode, i);
        }
        ImGui::PopItemWidth();
        draw_mode_ = mode;

        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Checkbox("Draw Control Mesh", &render_control_mesh_);
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Bezier", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Text("Tesselation Resolution");
        ImGui::PushItemWidth(120);
        int tesres = tesselation_resolution_;
        ImGui::SliderInt("  ", &tesres, 2, 30);
        ImGui::PopItemWidth();
        if (tesres != tesselation_resolution_)
        {
            tesselation_resolution_ = tesres;
            bezier_.tessellate(tesselation_resolution_);
        }

        static int last_tesselation_res = tesselation_resolution_;
        if (ImGui::Button("Tesselate"))
        {
            bezier_.tessellate(tesselation_resolution_);
            last_tesselation_res = tesselation_resolution_;
        }
        ImGui::Spacing();
        ImGui::Spacing();

        int t_mode = (int)bezier_.get_bezier_mode();
        ImGui::RadioButton("Bernstein Evaluation", &t_mode, 0);
        ImGui::RadioButton("Casteljau Evaluation", &t_mode, 1);
        if (t_mode != (int)bezier_.get_bezier_mode())
        {
            bezier_.toggle_de_Casteljau();
            bezier_.tessellate(tesselation_resolution_);
        }

        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Text("Tesselation time:\n%.2fms", bezier_.tesselation_time_);
    }
}

//-----------------------------------------------------------------------------

void BezierViewer::draw(const std::string &_draw_mode)
{
    // load shader?
    if (!phong_shader_.is_valid())
    {
        if (!phong_shader_.source(phong_vshader, phong_fshader))
            exit(1);
    }

    // clear framebuffer and depth buffer first
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // setup matrices
    mat4 mv_matrix = modelview_matrix_;
    mat4 mvp_matrix = projection_matrix_ * mv_matrix;
    mat3 n_matrix = inverse(transpose(linear_part(mv_matrix)));

    // render filled surface triangles with Phong lighting
    phong_shader_.use();
    phong_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
    phong_shader_.set_uniform("modelview_matrix", mv_matrix);
    phong_shader_.set_uniform("normal_matrix", n_matrix);
    phong_shader_.set_uniform("light1", normalize(vec3(1.0, 1.0, 1.0)));
    phong_shader_.set_uniform("light2", normalize(vec3(-1.0, 1.0, 1.0)));
    phong_shader_.set_uniform("front_color", vec3(0.9, 0.0, 0.0));
    phong_shader_.set_uniform("back_color", vec3(0.5, 0.0, 0.0));
    phong_shader_.set_uniform("ambient", 0.1f);
    phong_shader_.set_uniform("diffuse", 0.8f);
    phong_shader_.set_uniform("specular", 0.6f);
    phong_shader_.set_uniform("shininess", 100.0f);
    phong_shader_.set_uniform("alpha", 1.0f);
    phong_shader_.set_uniform("use_lighting", true);
    phong_shader_.set_uniform("use_texture", false);
    phong_shader_.set_uniform("use_srgb", false);
    phong_shader_.set_uniform("show_texture_layout", false);

    glDepthRange(0.01, 1.0); // reduce depth range for solid rendering,
    // such that wireframe will be slightly in front
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    bezier_.draw_surface(_draw_mode);
    glDepthRange(0.0, 1.0); // restore default depth range

    // render surface edges with lighting, but a bit darker
    if (_draw_mode == "Hidden Line")
    {
        phong_shader_.set_uniform("front_color", vec3(0.2, 0.2, 0.2));
        phong_shader_.set_uniform("back_color", vec3(0.2, 0.0, 0.0));
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        bezier_.draw_surface(_draw_mode);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    // render control mesh in grey
    if (render_control_mesh_)
    {
        phong_shader_.use();
        phong_shader_.set_uniform("use_lighting", false);
        phong_shader_.set_uniform("front_color", vec3(0.4));
        phong_shader_.set_uniform("back_color", vec3(0.4));
        bezier_.draw_control_polygon();
    }

    glCheckError();
}

//-----------------------------------------------------------------------------

void BezierViewer::mouse(int button, int action, int mods)
{
    // CTRL not pressed -> rotate model
    if (!ctrl_pressed())
    {
        TrackballViewer::mouse(button, action, mods);
        return;
    }

    // CTRL pressed -> drag control points around
    if (action == GLFW_PRESS)
    {
        // select joint
        if (button == GLFW_MOUSE_BUTTON_LEFT)
        {
            // get point under mouse cursor
            double x, y;
            cursor_pos(x, y);

            // get viewport paramters
            vec4 vp(0.0f);
            glGetFloatv(GL_VIEWPORT, &vp[0]);

            // in OpenGL y=0 is at the 'bottom'
            y = vp[3] - y;

            // invert viewport mapping
            vec2 p2d;
            p2d[0] = ((float)x - (float)vp[0]) / ((float)vp[2]) * 2.0f - 1.0f;
            p2d[1] = ((float)y - (float)vp[1]) / ((float)vp[3]) * 2.0f - 1.0f;

            // find closest joint
            mat4 mvp = projection_matrix_ * modelview_matrix_;
            bezier_.pick(vec2(p2d[0], p2d[1]), mvp);
        }
    }
}

//-----------------------------------------------------------------------------

void BezierViewer::motion(double _x, double _y)
{
    // CTRL not pressed -> rotate model
    if (!ctrl_pressed())
    {
        TrackballViewer::motion(_x, _y);
        return;
    }

    if (left_mouse_pressed())
    {
        vec3 p = bezier_.get_selected_control_point();

        // project joint position to 2D
        mat4 mvp = projection_matrix_ * modelview_matrix_;
        vec4 pp = mvp * vec4(p[0], p[1], p[2], 1.0f);
        pp /= pp[3];

        // get viewport data
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        int x = (int)_x;
        int y = (int)_y;

        // in OpenGL y=0 is at the 'bottom'
        y = viewport[3] - y;

        // unproject mouse position to 3D
        vec4 pn;
        pn[0] = ((float)x - (float)viewport[0]) / ((float)viewport[2]) * 2.0f -
                1.0f;
        pn[1] = ((float)y - (float)viewport[1]) / ((float)viewport[3]) * 2.0f -
                1.0f;
        // use z of projected joint position for correct depth
        pn[2] = pp[2];
        pn[3] = 1.0f;

        // unproject
        pn = inverse(mvp) * pn;
        pn /= pn[3];

        // set target position for effector and solve IK
        vec3 target(pn[0], pn[1], pn[2]);
        bezier_.set_selected_control_point(target);
        bezier_.tessellate(tesselation_resolution_);
    }
}

//=============================================================================
