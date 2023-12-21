//=============================================================================
//
//   Exercise code for the lecture "Computer Graphics"
//     by Prof. Mario Botsch, TU Dortmund
//
//   Copyright (C)  Computer Graphics Group, TU Dortmund
//
//=============================================================================

#include <imgui.h>
#include "Subdivision_Viewer.h"
#include <cfloat>
#include <iostream>
#include <sstream>

//=============================================================================
Subdivision_Viewer::Subdivision_Viewer(const char* title, int width, int height,
                                       bool showgui)
    : pmp::MeshViewer(title, width, height, showgui)
{
    // setup draw modes
    clear_draw_modes();
    add_draw_mode("Hidden Line");
    add_draw_mode("Smooth Shading");
    set_draw_mode("Hidden Line");
    crease_angle_ = 30.0;
    mesh_index_ = 0;
    draw_control_mesh_ = false;

    // add imgui help items
    add_help_item("S", "Subdivide", 0);
    add_help_item("C", "Toggle drawing the original mesh", 1);
    add_help_item("O", "Write mesh to an .off file", 2);
}

//-----------------------------------------------------------------------------

bool Subdivision_Viewer::load_mesh(const char* filename)
{
    bool success = MeshViewer::load_mesh(filename);

    if (success)
    {
        surface_mesh_.read(filename);

        // update scene center and bounds
        pmp::BoundingBox bb = mesh_.bounds();
        set_scene(bb.center(), 0.7 * bb.size());

        mesh_.set_crease_angle(crease_angle_);
        surface_mesh_.set_crease_angle(crease_angle_);

        // compute face & vertex normals, update face indices
        surface_mesh_.update_opengl_buffers();

        std::cout << std::endl;
        return true;
    }

    std::cerr << "Failed to read mesh from " << filename << " !" << std::endl;
    return false;
}

//-----------------------------------------------------------------------------

void Subdivision_Viewer::process_imgui()
{
    if (ImGui::CollapsingHeader("Load Mesh", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::PushItemWidth(100);
        const char* listbox_items[] = {"cube.off", "kissmouth.obj",
                                       "suzanne.obj", "teacup.off",
                                       "teapot.off"};
        int listbox_item_current = mesh_index_;
        ImGui::ListBox(" ", &listbox_item_current, listbox_items,
                       IM_ARRAYSIZE(listbox_items), 5);
        if (listbox_item_current != mesh_index_ || mesh_.n_vertices() == 0)
        {
            std::stringstream ss;
            ss << DATA_PATH << listbox_items[listbox_item_current];

            load_mesh(ss.str().c_str());

            mesh_index_ = listbox_item_current;
        }
        ImGui::PopItemWidth();
        ImGui::Spacing();
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Mesh Info", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // output mesh statistics
        ImGui::BulletText("%d vertices", (int)surface_mesh_.n_vertices());
        ImGui::BulletText("%d edges", (int)surface_mesh_.n_edges());
        ImGui::BulletText("%d faces", (int)surface_mesh_.n_faces());

        ImGui::Spacing();
        ImGui::Spacing();

        // control crease angle
        ImGui::PushItemWidth(80);
        ImGui::SliderFloat("Crease Angle", &crease_angle_, 0.0f, 180.0f,
                           "%.0f");
        ImGui::PopItemWidth();
        if (crease_angle_ != surface_mesh_.crease_angle())
        {
            surface_mesh_.set_crease_angle(crease_angle_);
        }
        ImGui::Spacing();
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Drawmode", ImGuiTreeNodeFlags_DefaultOpen))
    {
        int mode = draw_mode_;
        ImGui::PushItemWidth(100);
        for (int i = 0; i < (int)draw_mode_names_.size(); ++i)
        {
            ImGui::RadioButton(draw_mode_names_[i].c_str(), &mode, i);
        }
        ImGui::PopItemWidth();
        draw_mode_ = mode;

        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Checkbox("Draw Control Mesh", &draw_control_mesh_);
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Subdivision", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Spacing();
        if (ImGui::Button("Subdivide"))
        {
            surface_mesh_.subdivide();
        }
    }
}

//-----------------------------------------------------------------------------

void Subdivision_Viewer::draw(const std::string& drawMode)
{
    // draw mesh
    surface_mesh_.draw(projection_matrix_, modelview_matrix_, drawMode);
    if (draw_control_mesh_)
    {
        mesh_.draw(projection_matrix_, modelview_matrix_, "Edges");
    }
}

//-----------------------------------------------------------------------------

void Subdivision_Viewer::keyboard(int key, int scancode, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
        return;

    switch (key)
    {
        case GLFW_KEY_S: // subdivide model
        {
            surface_mesh_.subdivide();
            break;
        }

        case GLFW_KEY_C: // toggle control mesh
        {
            draw_control_mesh_ = !draw_control_mesh_;
            break;
        }

        case GLFW_KEY_O: // write mesh
        {
            surface_mesh_.write("output.off");
            break;
        }

        default:
        {
            MeshViewer::keyboard(key, scancode, action, mods);
            break;
        }
    }
}

//=============================================================================
