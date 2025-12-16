/**
 * P2dSVGDSimulator_t.cpp
 * * * Standalone 2D SVGD Simulation, 2D only.
 * * FINAL MOUSE CONTROLS: Middle=Move Held Peak, Right=Add Peak, Left=Remove Peak.
 * * Peak handles changed to CYAN. Panning is removed.
 * * Contains ALL class definitions to fix "no member named" errors.
 */

#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <numeric>
#include <functional>

// Include your headers (Assuming these exist and define Vec2, P2dGenerator, and SVGDCaculator)
#include "P2dGenerator.hpp"
#include "SVGDCaculator.hpp"

using Vec = Vec2;
using Generator = P2dGenerator;

// ==========================================
// 1. Graphics & Palette
// ==========================================
struct Color {
    float r, g, b, a;
    void apply() const { glColor4f(r, g, b, a); }
};

namespace Palette {
    const Color BackgroundPart(0.0f, 0.0f, 0.0f, 1.0f);
    const Color BackgroundDist(0.05f, 0.05f, 0.1f, 1.0f);
    const Color Particle(1.0f, 1.0f, 1.0f, 1.0f);
    const Color PeakHandle(0.0f, 1.0f, 1.0f, 1.0f); // CYAN

    Color get_contour_color(double val) {
        val = std::max(0.0, std::min(1.0, val));
        float r=0,g=0,b=0;
        if(val < 0.5) { float t=val/0.5f; r=t; g=1.0f; b=0.0f; }
        else { float t=(val-0.5f)/0.5f; r=1.0f; g=1.0f-t; b=0.0f; }
        return Color(r, g, b, 1.0f);
    }

    Color get_heatmap_color(double val) {
        val = std::max(0.0, std::min(1.0, val));
        float r=0,g=0,b=0;
        if (val < 0.25) { float t=val/0.25f; b=1.0f; g=t; }
        else if (val < 0.5) { float t=(val-0.25f)/0.25f; b=1.0f-t; g=1.0f; }
        else if (val < 0.75) { float t=(val-0.5f)/0.25f; g=1.0f; r=t; }
        else { float t=(val-0.75f)/0.25f; g=1.0f-t; r=1.0f; }
        return Color(r, g, b, 1.0f);
    }
}

// ==========================================
// 2. Logic Classes
// ==========================================

// --- Target Distribution (2D Gaussian Mixture) ---
class GaussianMixture2D {
public:
    struct Peak { Vec2 pos; double sigma; double weight; };
    std::vector<Peak> peaks;

    GaussianMixture2D() {
        peaks.push_back({Vec2(0,0), 1.5, 1.0});
    }

    double prob(const Vec2& x) const {
        double total = 0.0;
        for (const auto& p : peaks) {
            double sq_dist = (x - p.pos).normSq();
            total += p.weight * std::exp(-sq_dist / (2.0 * p.sigma * p.sigma));
        }
        return total;
    }

    Vec2 score(const Vec2& x) const {
        double den = 0.0; Vec2 grad_num;
        for (const auto& p : peaks) {
            Vec2 diff = x - p.pos; double var = p.sigma * p.sigma;
            double val = p.weight * std::exp(-diff.normSq() / (2.0 * var));
            den += val;
            grad_num = grad_num + (diff * (-1.0 / var)) * val;
        }
        if (den < 1e-12) return Vec2();
        return grad_num * (1.0 / den);
    }
};

// --- Particle Density Grid (2D) ---
class ParticleDensityGrid {
    int res; double min_xy, max_xy;
    std::vector<double> grid;
    int idx(int x, int y) const { return std::max(0, std::min(res-1, y))*res + std::max(0, std::min(res-1, x)); }
public:
    ParticleDensityGrid(int r, double min, double max) : res(r), min_xy(min), max_xy(max) { grid.resize(r*r, 0.0); }

    void compute(const std::vector<Vec2>& particles) {
        std::fill(grid.begin(), grid.end(), 0.0);
        double range = max_xy - min_xy;
        double max_val = 0.0;

        for (const auto& p : particles) {
            int ix = (int)((p.x - min_xy) / range * res);
            int iy = (int)((p.y - min_xy) / range * res);
            if (ix >= 0 && ix < res && iy >= 0 && iy < res) grid[idx(ix, iy)] += 1.0;
        }
        // Smooth
        auto smooth = [&]() {
            std::vector<double> temp = grid;
            for (int y=1; y<res-1; ++y) for (int x=1; x<res-1; ++x) {
                double s = 0;
                for(int dy=-1; dy<=1; ++dy) for(int dx=-1; dx<=1; ++dx) s += temp[idx(x+dx, y+dy)];
                grid[idx(x, y)] = s / 9.0;
            }
        };
        smooth(); smooth();

        for(double v : grid) max_val = std::max(max_val, v);
        if(max_val > 0) for(double& v : grid) v /= max_val;
    }

    double get(double wx, double wy) const {
        double range = max_xy - min_xy;
        int ix = (int)((wx - min_xy) / range * res);
        int iy = (int)((wy - min_xy) / range * res);
        if (ix < 0 || ix >= res || iy < 0 || iy >= res) return 0.0;
        return grid[idx(ix, iy)];
    }
};

// --- Simulation Engine (2D) ---
class Simulation {
public:
    // *** THESE ARE THE MEMBERS THE COMPILER WAS COMPLAINING ABOUT! ***
    std::unique_ptr<Generator> generator;
    std::unique_ptr<GaussianMixture2D> target;
    std::unique_ptr<ParticleDensityGrid> density;

    Simulation() {
        generator = std::make_unique<Generator>(-10, 10, -10, 10);
        generator->generate_uniform(300);
        target = std::make_unique<GaussianMixture2D>();
        density = std::make_unique<ParticleDensityGrid>(60, -10.0, 10.0);
    }

    void update(double dt) {
        generator->swim(dt);
        auto pos = generator->get_positions();
        SVGDCaculator<Vec2>::step(pos, [&](const Vec2& x){ return target->score(x); }, 0.5);
        generator->set_positions(pos);
        density->compute(pos);
    }
};

// ==========================================
// 3. Application
// ==========================================
class App {
public:
    GLFWwindow* win_particles;
    GLFWwindow* win_target;
    // The member 'sim' is correctly defined here:
    Simulation sim;

    // Interaction State
    bool dragging_peak = false;
    int active_peak = -1;
    int width = 600, height = 600;

    App() { init_windows(); }
    ~App() { glfwTerminate(); }

    void init_windows() {
        if(!glfwInit()) exit(-1);

        win_particles = glfwCreateWindow(width, height, "Window 1: 2D Particle Density & Contours", NULL, NULL);
        glfwSetWindowPos(win_particles, 100, 100);

        win_target = glfwCreateWindow(width, height, "Window 2: Target Heatmap (Middle=Move Peak, Right=Add, Left=Remove)", NULL, NULL);
        glfwSetWindowPos(win_target, 100 + width + 20, 100);
        glfwSetWindowUserPointer(win_target, this);

        // Input Callbacks
        glfwSetMouseButtonCallback(win_target, [](GLFWwindow* w, int b, int a, int m) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->on_click(b, a);
        });
        glfwSetCursorPosCallback(win_target, [](GLFWwindow* w, double x, double y) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->on_drag(x, y);
        });
    }

    void run() {
        while(!glfwWindowShouldClose(win_particles) && !glfwWindowShouldClose(win_target)) {
            // Error was here: sim.update is now correctly defined in class Simulation above
            sim.update(0.001);

            // --- Window 1: 2D Particles (Always 2D) ---
            glfwMakeContextCurrent(win_particles);
            Palette::BackgroundPart.apply();
            glClear(GL_COLOR_BUFFER_BIT);
            glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(-10,10,-10,10,-1,1);
            glMatrixMode(GL_MODELVIEW); glLoadIdentity();
            render_particle_view_2d();
            glfwSwapBuffers(win_particles);

            // --- Window 2: 2D Heatmap (Target Distribution) ---
            glfwMakeContextCurrent(win_target);
            Palette::BackgroundDist.apply();
            glClear(GL_COLOR_BUFFER_BIT);

            glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(-10,10,-10,10,-1,1);
            glMatrixMode(GL_MODELVIEW); glLoadIdentity();

            render_2d_target_heatmap();

            glfwSwapBuffers(win_target);
            glfwPollEvents();
        }
    }

private:

    Vec2 screen_to_world(double mx, double my) const {
        double range = 20.0;
        double wx = (mx / width) * range - 10.0;
        double wy = -((my / height) * range - 10.0);
        return Vec2(wx, wy);
    }

    // --- Interaction ---
    void on_click(int btn, int action) {
        double mx, my; glfwGetCursorPos(win_target, &mx, &my);
        Vec2 world_pos = screen_to_world(mx, my);

        const double PEAK_TOLERANCE_SQ = 0.5 * 0.5;

        if(action == GLFW_PRESS) {
            if(btn == GLFW_MOUSE_BUTTON_RIGHT) {
                // Right Click: ADD a new Gaussian peak
                sim.target->peaks.push_back({world_pos, 1.5, 1.0}); // Error fixed
            }
            else if(btn == GLFW_MOUSE_BUTTON_LEFT) {
                // Left Click: REMOVE the closest Gaussian peak
                double min_d_sq = PEAK_TOLERANCE_SQ;
                int peak_to_remove = -1;

                for(size_t i=0; i<sim.target->peaks.size(); ++i) { // Error fixed
                    double d_sq = (sim.target->peaks[i].pos - world_pos).normSq(); // Error fixed
                    if(d_sq < min_d_sq) {
                        min_d_sq = d_sq;
                        peak_to_remove = i;
                    }
                }

                if(peak_to_remove != -1 && sim.target->peaks.size() > 1) {
                    sim.target->peaks.erase(sim.target->peaks.begin() + peak_to_remove); // Error fixed
                }
            }
            else if(btn == GLFW_MOUSE_BUTTON_MIDDLE) {
                // Middle Click: Select a peak to move (Start drag)
                double min_d_sq = PEAK_TOLERANCE_SQ;
                active_peak = -1;
                for(size_t i=0; i<sim.target->peaks.size(); ++i) { // Error fixed
                    double d_sq = (sim.target->peaks[i].pos - world_pos).normSq(); // Error fixed
                    if(d_sq < min_d_sq) {
                        min_d_sq = d_sq;
                        active_peak = i;
                    }
                }
                if(active_peak != -1) dragging_peak = true;
            }
        } else if (action == GLFW_RELEASE) {
            if(btn == GLFW_MOUSE_BUTTON_MIDDLE) {
                // Middle Release: Stop moving peak
                dragging_peak = false;
                active_peak = -1;
            }
        }
    }

    void on_drag(double x, double y) {
        if(dragging_peak && active_peak != -1) {
            // Middle Drag: Move the active peak
            Vec2 world_pos = screen_to_world(x, y);
            if(active_peak < sim.target->peaks.size()) // Error fixed
                sim.target->peaks[active_peak].pos = world_pos; // Error fixed
        }
    }

    // --- Rendering Functions ---

    void render_particle_view_2d() {
        // Draw Contours (Empirical Density - Window 1)
        glLineWidth(1.5f); glBegin(GL_LINES);
        int res = 60; double step = 20.0/res;
        std::vector<double> levels = {0.1, 0.3, 0.5, 0.7, 0.9};

        for(int i=0; i<res-1; ++i) {
            for(int j=0; j<res-1; ++j) {
                double x=-10.0 + i*step, y=-10.0 + j*step;
                // Error fixed: sim.density is now defined
                double v0 = sim.density->get(x, y);
                double v1 = sim.density->get(x+step, y); double v2 = sim.density->get(x+step, y+step); double v3 = sim.density->get(x, y+step);

                for(double lvl : levels) {
                    if(lvl > std::max({v0,v1,v2,v3}) || lvl < std::min({v0,v1,v2,v3})) continue;
                    Palette::get_contour_color(lvl).apply();
                    auto t = [&](double a, double b){ return (lvl-a)/(b-a+1e-5); };

                    double p[8]; int c=0;
                    if((v0<lvl)!=(v1<lvl)) { p[c]=x+t(v0,v1)*step; p[c+1]=y; c+=2; }
                    if((v1<lvl)!=(v2<lvl)) { p[c]=x+step; p[c+1]=y+t(v1,v2)*step; c+=2; }
                    if((v2<lvl)!=(v3<lvl)) { p[c]=x+(1-t(v3,v2))*step; p[c+1]=y+step; c+=2; }
                    if((v3<lvl)!=(v0<lvl)) { p[c]=x; p[c+1]=y+t(v0,v3)*step; c+=2; }
                    if(c>=4) { glVertex2f(p[0], p[1]); glVertex2f(p[2], p[3]); }
                }
            }
        }
        glEnd();

        // Draw Particles
        glPointSize(3.0); Palette::Particle.apply(); glBegin(GL_POINTS);
        sim.generator->render([](Vec2 p){ glVertex2f(p.x, p.y); }); // Error fixed
        glEnd();
    }

    void render_2d_target_heatmap() {
        // Draw 2D Heatmap (Target Distribution - Window 2)
        glShadeModel(GL_SMOOTH); glBegin(GL_QUADS);
        int res = 40; double step = 20.0/res;

        for(int i=0; i<res; ++i) {
            for(int j=0; j < res; ++j) {
                double x=-10+i*step, y=-10+j*step;
                // Error fixed: sim.target is now defined
                Palette::get_heatmap_color(sim.target->prob({x,y})).apply(); glVertex2f(x,y);
                Palette::get_heatmap_color(sim.target->prob({x+step,y})).apply(); glVertex2f(x+step,y);
                Palette::get_heatmap_color(sim.target->prob({x+step,y+step})).apply(); glVertex2f(x+step,y+step);
                Palette::get_heatmap_color(sim.target->prob({x,y+step})).apply(); glVertex2f(x,y+step);
            }
        }
        glEnd();

        // Draw Handles (CYAN)
        glPointSize(8.0); Palette::PeakHandle.apply(); glBegin(GL_POINTS);
        for(const auto& p : sim.target->peaks) glVertex2f(p.pos.x, p.pos.y); // Error fixed
        glEnd();
    }
};

// ==========================================
// 4. Entry Point
// ==========================================
int main() {
    App app;
    app.run();
    return 0;
}