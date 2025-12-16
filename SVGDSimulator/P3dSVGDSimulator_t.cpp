/**
 * P3dSVGDSimulator_t.cpp
 * * Standalone 3D SVGD Simulation with two windows.
 * * Window 1: 2D XY Projection (Non-interactive monitor).
 * * Window 2: Interactive 3D Particle View
 * * MOUSE: LClick Drag=Orbit, MClick=Select/Move Peak, RClick=Add Peak, LClick=Remove Peak.
 * * KEYBOARD: When a peak is selected (MClick), use ARROW KEYS (Up/Down/Left/Right)
 * * and W/S/A/D to move the peak center.
 */

#include <GLFW/glfw3.h>
#include <GL/glu.h> // REQUIRED for gluProject and gluUnProject
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <numeric>
#include <functional>

// Include your headers (Assuming these exist and define Vec3, Vec2, P3dGenerator, and SVGDCaculator)
#include "P3dGenerator.hpp"
#include "P2dGenerator.hpp"
#include "SVGDCaculator.hpp"

using Vec = Vec3;
using Generator = P3dGenerator;

// ==========================================
// 1. Graphics & Palette
// ==========================================
struct Color {
    float r, g, b, a;
    void apply() const { glColor4f(r, g, b, a); }
};

namespace Palette {
    const Color BackgroundSim(0.1f, 0.1f, 0.2f, 1.0f);
    const Color BackgroundProj(0.0f, 0.0f, 0.0f, 1.0f);
    const Color Particle(1.0f, 1.0f, 1.0f, 1.0f);
    const Color GridLines(0.3f, 0.3f, 0.3f, 1.0f);
    const Color PeakHandle(1.0f, 0.0f, 1.0f, 1.0f); // MAGENTA for 3D peaks
    const Color ActivePeakHighlight(1.0f, 1.0f, 0.0f, 1.0f); // YELLOW for selection
}

// ==========================================
// 2. Logic Classes (3D Distributions)
// ==========================================
class GaussianMixture3D {
public:
    struct Peak { Vec3 pos; double sigma; double weight; };
    std::vector<Peak> peaks;

    GaussianMixture3D() {
        peaks.push_back({Vec3(5, 5, 5), 1.5, 1.0});
        peaks.push_back({Vec3(-5, -5, -5), 1.5, 1.0});
    }

    double prob(const Vec3& x) const {
        double total = 0.0;
        for (const auto& p : peaks) {
            double sq_dist = (x - p.pos).normSq();
            total += p.weight * std::exp(-sq_dist / (2.0 * p.sigma * p.sigma));
        }
        return total;
    }

    Vec3 score(const Vec3& x) const {
        double den = 0.0; Vec3 grad_num;
        for (const auto& p : peaks) {
            Vec3 diff = x - p.pos; double var = p.sigma * p.sigma;
            double val = p.weight * std::exp(-diff.normSq() / (2.0 * var));
            den += val;
            grad_num = grad_num + (diff * (-1.0 / var)) * val;
        }
        if (den < 1e-12) return Vec3();
        return grad_num * (1.0 / den);
    }
};

class Simulation {
public:
    std::unique_ptr<Generator> generator;
    std::unique_ptr<GaussianMixture3D> target;

    Simulation() {
        generator = std::make_unique<Generator>(-10, 10, -10, 10, -10, 10);
        generator->generate_uniform(500);
        target = std::make_unique<GaussianMixture3D>();
    }

    void update(double dt) {
        generator->swim(dt);
        auto pos = generator->get_positions();
        SVGDCaculator<Vec3>::step(pos, [&](const Vec3& x){ return target->score(x); }, 0.5);
        generator->set_positions(pos);
    }
};

// ==========================================
// 3. 3D Camera & View Setup
// ==========================================
struct OrbitCamera {
    float theta = 0.5f;
    float phi = 0.5f;
    float zoom = 30.0f;

    void drag(float dx, float dy) {
        theta += dx * 0.01f;
        phi   += dy * 0.01f;
        phi = std::max(0.01f, std::min(3.14f - 0.01f, phi));
    }

    void scroll(float dy) {
        zoom = std::max(5.0f, zoom + dy * -1.0f);
    }

    void apply_view() const {
        // Step 1: Translate camera away
        glTranslatef(0.0f, 0.0f, -zoom);
        // Step 2: Apply rotation (phi for pitch, theta for yaw)
        glRotatef(phi * 180.0f / 3.14159f, 1.0f, 0.0f, 0.0f);
        glRotatef(theta * 180.0f / 3.14159f, 0.0f, 1.0f, 0.0f);
    }
};

void set_perspective(int width, int height) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double aspect = (double)width / height;
    double fov = 45.0;
    double zNear = 0.1;
    double zFar = 100.0;
    double fH = tan(fov / 360 * 3.14159) * zNear;
    double fW = fH * aspect;
    glFrustum(-fW, fW, -fH, fH, zNear, zFar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// ==========================================
// 4. Application
// ==========================================
class App {
public:
    GLFWwindow* win_projection;
    GLFWwindow* win_3d;
    Simulation sim;

    // Interaction State
    OrbitCamera camera;
    int width = 600, height = 600;

    // Peak Editing State
    bool dragging_peak = false;
    int active_peak = -1;
    double active_peak_depth = 0.0;
    const double KEY_MOVE_STEP = 0.1; // Step size for keyboard movement

    App() { init_windows(); }
    ~App() { glfwTerminate(); }

    void init_windows() {
        if(!glfwInit()) exit(-1);

        win_projection = glfwCreateWindow(width, height, "Window 1: 2D XY Projection", NULL, NULL);
        glfwSetWindowPos(win_projection, 100, 100);

        win_3d = glfwCreateWindow(width, height, "Window 2: Interactive 3D (LClick=Orbit, MClick=Select, Arrows=Move)", NULL, NULL);
        glfwSetWindowPos(win_3d, 100 + width + 20, 100);
        glfwSetWindowUserPointer(win_3d, this);

        // Input Callbacks for 3D Window
        glfwSetMouseButtonCallback(win_3d, [](GLFWwindow* w, int b, int a, int m) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->on_click(b, a);
        });
        glfwSetCursorPosCallback(win_3d, [](GLFWwindow* w, double x, double y) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->on_drag(x, y);
        });
        glfwSetScrollCallback(win_3d, [](GLFWwindow* w, double dx, double dy) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->camera.scroll((float)dy);
        });
        // Keyboard Callback
        glfwSetKeyCallback(win_3d, [](GLFWwindow* w, int key, int scancode, int action, int mods) {
            if (action == GLFW_PRESS || action == GLFW_REPEAT) {
                static_cast<App*>(glfwGetWindowUserPointer(w))->on_key_press(key);
            }
        });
    }

    void run() {
        while(!glfwWindowShouldClose(win_projection) && !glfwWindowShouldClose(win_3d)) {
            sim.update(0.001);

            // --- Window 1: 2D Projection ---
            glfwMakeContextCurrent(win_projection);
            Palette::BackgroundProj.apply();
            glClear(GL_COLOR_BUFFER_BIT);
            glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(-10,10,-10,10,-1,1);
            glMatrixMode(GL_MODELVIEW); glLoadIdentity();
            render_2d_projection();
            glfwSwapBuffers(win_projection);

            // --- Window 2: Interactive 3D View ---
            glfwMakeContextCurrent(win_3d);
            Palette::BackgroundSim.apply();
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glEnable(GL_DEPTH_TEST);

            set_perspective(width, height);

            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            camera.apply_view();

            render_3d_domain();
            render_3d_particles();

            glDisable(GL_DEPTH_TEST);
            glfwSwapBuffers(win_3d);
            glfwPollEvents();
        }
    }

private:

    // --- Utility: 3D to 2D Projection ---
    Vec2 project_3d(const Vec3& p) {
        GLdouble modelview[16]; GLdouble projection[16]; GLint viewport[4]; GLdouble winX, winY, winZ;
        glGetDoublev(GL_MODELVIEW_MATRIX, modelview); glGetDoublev(GL_PROJECTION_MATRIX, projection);
        glGetIntegerv(GL_VIEWPORT, viewport);
        gluProject(p.x, p.y, p.z, modelview, projection, viewport, &winX, &winY, &winZ);
        return Vec2((double)winX, (double)height - winY);
    }

    // --- Utility: 2D Screen to 3D World (on a fixed Z-plane) ---
    Vec3 unproject_2d(double mx, double my, double world_z) {
        GLdouble modelview[16]; GLdouble projection[16]; GLint viewport[4]; GLdouble worldX, worldY, worldZ;
        glGetDoublev(GL_MODELVIEW_MATRIX, modelview); glGetDoublev(GL_PROJECTION_MATRIX, projection);
        glGetIntegerv(GL_VIEWPORT, viewport);
        double opengl_y = height - my;
        gluUnProject(mx, opengl_y, 0.5, modelview, projection, viewport, &worldX, &worldY, &worldZ);
        double final_z = dragging_peak ? active_peak_depth : world_z;
        return Vec3((double)worldX, (double)worldY, final_z);
    }


    // --- Interaction ---
    void on_click(int btn, int action) {
        double mx, my; glfwGetCursorPos(win_3d, &mx, &my);
        const double PEAK_TOLERANCE_SQ = 15.0 * 15.0;

        if(action == GLFW_PRESS) {
            // Modified code (Always adds center at (0, 0, 0))
            if(btn == GLFW_MOUSE_BUTTON_RIGHT) {
                // RClick: ADD a new Gaussian peak at the origin (0, 0, 0)
                // No need to use mouse coordinates (mx, my) or unproject_2d
                Vec3 world_pos(0.0, 0.0, 0.0); // <-- HARDCODED ORIGIN

                sim.target->peaks.push_back({world_pos, 1.5, 1.0});
                active_peak = sim.target->peaks.size() - 1;
            }
            else if(btn == GLFW_MOUSE_BUTTON_MIDDLE) {
                // MClick: Select a peak (Start drag/selection)
                double min_d_sq = PEAK_TOLERANCE_SQ;
                int candidate_peak = -1;

                for(size_t i=0; i<sim.target->peaks.size(); ++i) {
                    Vec2 peak_screen_pos = project_3d(sim.target->peaks[i].pos);
                    double d_sq = (peak_screen_pos - Vec2((double)mx, (double)my)).normSq();

                    if(d_sq < min_d_sq) { min_d_sq = d_sq; candidate_peak = i; }
                }

                if(candidate_peak != -1) {
                    active_peak = candidate_peak; // Set selection
                    dragging_peak = true;        // Start mouse dragging (for immediate movement)
                    active_peak_depth = sim.target->peaks[active_peak].pos.z;
                } else {
                    active_peak = -1; // Deselect if click missed
                }
            }
            else if (btn == GLFW_MOUSE_BUTTON_LEFT) {
                // LClick: REMOVE the closest Gaussian peak
                double min_d_sq = PEAK_TOLERANCE_SQ;
                int peak_to_remove = -1;

                for(size_t i=0; i<sim.target->peaks.size(); ++i) {
                    Vec2 peak_screen_pos = project_3d(sim.target->peaks[i].pos);
                    double d_sq = (peak_screen_pos - Vec2((double)mx, (double)my)).normSq();
                    if(d_sq < min_d_sq) { min_d_sq = d_sq; peak_to_remove = i; }
                }

                if(peak_to_remove != -1 && sim.target->peaks.size() > 1) {
                    sim.target->peaks.erase(sim.target->peaks.begin() + peak_to_remove);
                    if (active_peak == peak_to_remove) active_peak = -1;
                    else if (active_peak > peak_to_remove) active_peak--;
                } else {
                    active_peak = -1; // Deselect if click missed the removal threshold
                }
            }
        } else if (action == GLFW_RELEASE) {
            if(btn == GLFW_MOUSE_BUTTON_MIDDLE) {
                dragging_peak = false;
                // Keep active_peak set for keyboard control
            }
        }
    }

    void on_drag(double x, double y) {
        static double last_x = 0, last_y = 0;
        double dx = x - last_x;
        double dy = y - last_y;
        last_x = x; last_y = y;

        // Middle Drag: Move the active peak
        if(dragging_peak && active_peak != -1) {
            Vec3 world_pos = unproject_2d(x, y, active_peak_depth);
            if(active_peak < sim.target->peaks.size())
                sim.target->peaks[active_peak].pos = world_pos;
        }
        // Left Drag: Orbital Rotation
        else if (glfwGetMouseButton(win_3d, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            camera.drag((float)dx, (float)dy);
        }
    }

    // --- New: Keyboard Control for Active Peak ---
    void on_key_press(int key) {
        if (active_peak == -1 || active_peak >= sim.target->peaks.size()) return;

        double& px = sim.target->peaks[active_peak].pos.x;
        double& py = sim.target->peaks[active_peak].pos.y;
        double& pz = sim.target->peaks[active_peak].pos.z;

        // X-Y Plane Movement (Horizontal/Vertical)
        if (key == GLFW_KEY_LEFT) px -= KEY_MOVE_STEP;
        else if (key == GLFW_KEY_RIGHT) px += KEY_MOVE_STEP;
        else if (key == GLFW_KEY_UP) py += KEY_MOVE_STEP;
        else if (key == GLFW_KEY_DOWN) py -= KEY_MOVE_STEP;

        // Z-Depth Movement (W/S keys)
        else if (key == GLFW_KEY_W) pz += KEY_MOVE_STEP;
        else if (key == GLFW_KEY_S) pz -= KEY_MOVE_STEP;
    }


    // --- Rendering Functions ---

    void render_2d_projection() {
        glPointSize(3.0); Palette::Particle.apply(); glBegin(GL_POINTS);
        sim.generator->render([](Vec3 p){ glVertex2f(p.x, p.y); });
        glEnd();
    }

    void render_3d_domain() {
        Palette::GridLines.apply();
        glLineWidth(1.0f);

        float min = -10.0f, max = 10.0f;

        // Draw the main X-Z grid (Y=0 plane)
        glBegin(GL_LINES);
        glColor3f(0.5f, 0.5f, 0.5f);
        for(int i=-10; i<=10; i+=2) {
             // Lines parallel to Z-axis
             glVertex3f((float)i, 0.0f, min); glVertex3f((float)i, 0.0f, max);
             // Lines parallel to X-axis
             glVertex3f(min, 0.0f, (float)i); glVertex3f(max, 0.0f, (float)i);
        }
        glEnd();

        // Draw Bounding Box (Edges of the domain)
        glBegin(GL_LINES);
        glColor3f(0.3f, 0.3f, 0.3f);
        // X-parallel lines
        glVertex3f(min, min, min); glVertex3f(max, min, min);
        glVertex3f(min, max, min); glVertex3f(max, max, min);
        glVertex3f(min, min, max); glVertex3f(max, min, max);
        glVertex3f(min, max, max); glVertex3f(max, max, max);
        // Y-parallel lines
        glVertex3f(min, min, min); glVertex3f(min, max, min);
        glVertex3f(max, min, min); glVertex3f(max, max, min);
        glVertex3f(min, min, max); glVertex3f(min, max, max);
        glVertex3f(max, min, max); glVertex3f(max, max, max);
        // Z-parallel lines
        glVertex3f(min, min, min); glVertex3f(min, min, max);
        glVertex3f(max, min, min); glVertex3f(max, min, max);
        glVertex3f(min, max, min); glVertex3f(min, max, max);
        glVertex3f(max, max, min); glVertex3f(max, max, max);
        glEnd();

        // Draw Axes (Colored)
        glLineWidth(3.0f);
        glBegin(GL_LINES);
        // X-axis (Red)
        glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(0, 0, 0); glVertex3f(12, 0, 0);
        // Y-axis (Green)
        glColor3f(0.0f, 1.0f, 0.0f); glVertex3f(0, 0, 0); glVertex3f(0, 12, 0);
        // Z-axis (Blue)
        glColor3f(0.0f, 0.0f, 1.0f); glVertex3f(0, 0, 0); glVertex3f(0, 0, 12);
        glEnd();
    }

    void render_3d_particles() {
        // Draws the actual particle cloud
        glPointSize(4.0); Palette::Particle.apply(); glBegin(GL_POINTS);
        sim.generator->render([](Vec3 p){ glVertex3f(p.x, p.y, p.z); });
        glEnd();

        // Draw Target Peaks (Magenta)
        glPointSize(10.0); Palette::PeakHandle.apply(); glBegin(GL_POINTS);
        for(size_t i=0; i<sim.target->peaks.size(); ++i) {
            glVertex3f(sim.target->peaks[i].pos.x, sim.target->peaks[i].pos.y, sim.target->peaks[i].pos.z);
        }
        glEnd();

        // HIGHLIGHT the active peak (Yellow)
        if (active_peak != -1 && active_peak < sim.target->peaks.size()) {
            Palette::ActivePeakHighlight.apply();
            glPointSize(18.0); // Make the highlight larger than the peak dot
            glBegin(GL_POINTS);
            const auto& p = sim.target->peaks[active_peak].pos;
            glVertex3f(p.x, p.y, p.z);
            glEnd();
        }
    }
};

// ==========================================
// 5. Entry Point
// ==========================================
int main() {
    App app;
    app.run();
    return 0;
}