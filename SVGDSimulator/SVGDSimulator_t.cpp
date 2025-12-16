/**
 * main.cpp
 *
 * Refactored Dual-Window SVGD Simulation.
 *
 * Architecture:
 * 1. Palette: Color definitions.
 * 2. Math/Logic: GaussianMixture (Target) & ParticleDensityGrid (Empirical).
 * 3. Simulation: Manages the physics step and data state.
 * 4. App: Manages GLFW windows, input, and the main loop.
 */

#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>

// Include your headers
#include "P2dGenerator.hpp"
// #include "P3dGenerator.hpp" // Not used in this 2D dual-window setup
#include "SVGD.hpp"

// ==========================================
// 1. Graphics & Palette
// ==========================================
struct Color {
    float r, g, b, a;
    void apply() const { glColor4f(r, g, b, a); }
};

namespace Palette {
    const Color BackgroundPart(0.0f, 0.0f, 0.0f, 1.0f);   // Window 1 BG
    const Color BackgroundDist(0.05f, 0.05f, 0.1f, 1.0f); // Window 2 BG
    const Color Particle(1.0f, 1.0f, 1.0f, 1.0f);

    // Gradient: Green -> Yellow -> Red
    Color get_contour_color(double val) {
        val = std::max(0.0, std::min(1.0, val));
        float r=0,g=0,b=0;
        if(val < 0.5) { float t=val/0.5f; r=t; g=1.0f; b=0.0f; }
        else { float t=(val-0.5f)/0.5f; r=1.0f; g=1.0f-t; b=0.0f; }
        return Color(r, g, b, 1.0f);
    }

    // Gradient: Blue -> Green -> Red
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
// 2. Logic Classes (Distributions)
// ==========================================

// The "Target" Math Distribution
class GaussianMixture {
public:
    struct Peak { Vec2 pos; double sigma; double weight; };
    std::vector<Peak> peaks;

    GaussianMixture() {
        // Start with one peak at center
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
        double den = 0.0;
        Vec2 grad;
        for (const auto& p : peaks) {
            Vec2 diff = x - p.pos;
            double var = p.sigma * p.sigma;
            double val = p.weight * std::exp(-diff.normSq() / (2.0 * var));
            den += val;
            grad = grad + (diff * (-1.0 / var)) * val;
        }
        if (den < 1e-12) return Vec2();
        return grad * (1.0 / den);
    }
};

// The "Empirical" Density Grid (Histogram + Blur)
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

        // 1. Binning
        for (const auto& p : particles) {
            int ix = (int)((p.x - min_xy) / range * res);
            int iy = (int)((p.y - min_xy) / range * res);
            if (ix >= 0 && ix < res && iy >= 0 && iy < res) grid[idx(ix, iy)] += 1.0;
        }
        // 2. Blur (Simple 3x3 box blur x2)
        auto smooth = [&]() {
            std::vector<double> temp = grid;
            for (int y=1; y<res-1; ++y) for (int x=1; x<res-1; ++x) {
                double s = 0;
                for(int dy=-1; dy<=1; ++dy) for(int dx=-1; dx<=1; ++dx) s += temp[idx(x+dx, y+dy)];
                grid[idx(x, y)] = s / 9.0;
            }
        };
        smooth(); smooth();

        // 3. Normalize
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

// ==========================================
// 3. Simulation State & Physics
// ==========================================
class Simulation {
public:
    std::unique_ptr<P2dGenerator> generator;
    std::unique_ptr<GaussianMixture> target;
    std::unique_ptr<ParticleDensityGrid> density;

    Simulation() {
        generator = std::make_unique<P2dGenerator>(-10, 10, -10, 10);
        generator->generate_uniform(300);

        target = std::make_unique<GaussianMixture>();

        // 60x60 grid covering -10 to 10
        density = std::make_unique<ParticleDensityGrid>(60, -10.0, 10.0);
    }

    void update(double dt) {
        // 1. Swim
        generator->swim(dt);

        // 2. SVGD Step
        auto pos = generator->get_positions();
        SVGD<Vec2>::step(pos, [&](const Vec2& x){ return target->score(x); }, 0.5); // 0.5 is LR
        generator->set_positions(pos);

        // 3. Update Density
        density->compute(pos);
    }
};

// ==========================================
// 4. Application (Window & Input Manager)
// ==========================================
class App {
public:
    GLFWwindow* win_particles; // Left Window
    GLFWwindow* win_target;    // Right Window
    Simulation sim;

    // Interaction State
    bool dragging = false;
    int active_peak = -1;
    int width = 600, height = 600;

    App() {
        init_windows();
    }

    ~App() {
        glfwTerminate();
    }

    void init_windows() {
        if(!glfwInit()) exit(-1);

        win_particles = glfwCreateWindow(width, height, "Window 1: Particles & Contours", NULL, NULL);
        glfwSetWindowPos(win_particles, 100, 100);

        win_target = glfwCreateWindow(width, height, "Window 2: Target Heatmap (Right Click Add)", NULL, NULL);
        glfwSetWindowPos(win_target, 100 + width + 20, 100);

        // Set User Pointer to access 'this' inside callbacks
        glfwSetWindowUserPointer(win_target, this);

        // Callbacks for Target Window only
        glfwSetMouseButtonCallback(win_target, [](GLFWwindow* w, int b, int a, int m) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->on_click(b, a);
        });
        glfwSetCursorPosCallback(win_target, [](GLFWwindow* w, double x, double y) {
            static_cast<App*>(glfwGetWindowUserPointer(w))->on_drag(x, y);
        });
    }

    void run() {
        while(!glfwWindowShouldClose(win_particles) && !glfwWindowShouldClose(win_target)) {
            sim.update(0.1); // Update Physics

            // --- Render Window 1 ---
            glfwMakeContextCurrent(win_particles);
            Palette::BackgroundPart.apply();
            glClear(GL_COLOR_BUFFER_BIT);
            setup_viewport();
            render_particle_view();
            glfwSwapBuffers(win_particles);

            // --- Render Window 2 ---
            glfwMakeContextCurrent(win_target);
            Palette::BackgroundDist.apply();
            glClear(GL_COLOR_BUFFER_BIT);
            setup_viewport();
            render_target_view();
            glfwSwapBuffers(win_target);

            glfwPollEvents();
        }
    }

private:
    void setup_viewport() {
        glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(-10,10,-10,10,-1,1);
        glMatrixMode(GL_MODELVIEW); glLoadIdentity();
    }

    // --- Interaction Logic ---
    void on_click(int btn, int action) {
        double mx, my; glfwGetCursorPos(win_target, &mx, &my);
        double wx = (mx/width)*20.0 - 10.0;
        double wy = -((my/height)*20.0 - 10.0);

        if(btn == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
            sim.target->peaks.push_back({Vec2(wx, wy), 1.5, 1.0});
        }
        else if(btn == GLFW_MOUSE_BUTTON_LEFT) {
            if(action == GLFW_PRESS) {
                // Find closest peak
                double min_d = 4.0;
                active_peak = -1;
                for(size_t i=0; i<sim.target->peaks.size(); ++i) {
                    double d = (sim.target->peaks[i].pos - Vec2(wx, wy)).normSq();
                    if(d < min_d) { min_d = d; active_peak = i; }
                }
                if(active_peak != -1) dragging = true;
            } else {
                dragging = false;
            }
        }
    }

    void on_drag(double x, double y) {
        if(dragging && active_peak != -1) {
            double wx = (x/width)*20.0 - 10.0;
            double wy = -((y/height)*20.0 - 10.0);
            if(active_peak < sim.target->peaks.size())
                sim.target->peaks[active_peak].pos = Vec2(wx, wy);
        }
    }

    // --- Drawing Functions ---

    void render_particle_view() {
        // 1. Draw Contours (Marching Squares Lite)
        glLineWidth(1.5f); glBegin(GL_LINES);
        int res = 60; double step = 20.0/res;
        std::vector<double> levels = {0.1, 0.3, 0.5, 0.7, 0.9};

        for(int i=0; i<res-1; ++i) {
            for(int j=0; j<res-1; ++j) {
                double x=-10.0 + i*step, y=-10.0 + j*step;
                double v0 = sim.density->get(x, y);
                double v1 = sim.density->get(x+step, y);
                double v2 = sim.density->get(x+step, y+step);
                double v3 = sim.density->get(x, y+step);

                for(double lvl : levels) {
                    if(lvl > std::max({v0,v1,v2,v3}) || lvl < std::min({v0,v1,v2,v3})) continue;
                    Palette::get_contour_color(lvl).apply();
                    auto t = [&](double a, double b){ return (lvl-a)/(b-a+1e-5); };

                    std::vector<Vec2> p;
                    if((v0<lvl)!=(v1<lvl)) p.push_back(Vec2(x+t(v0,v1)*step, y));
                    if((v1<lvl)!=(v2<lvl)) p.push_back(Vec2(x+step, y+t(v1,v2)*step));
                    if((v2<lvl)!=(v3<lvl)) p.push_back(Vec2(x+(1-t(v3,v2))*step, y+step));
                    if((v3<lvl)!=(v0<lvl)) p.push_back(Vec2(x, y+t(v0,v3)*step));
                    if(p.size()>=2) { glVertex2f(p[0].x, p[0].y); glVertex2f(p[1].x, p[1].y); }
                }
            }
        }
        glEnd();

        // 2. Draw Particles
        glPointSize(3.0); Palette::Particle.apply(); glBegin(GL_POINTS);
        sim.generator->render([](Vec2 p){ glVertex2f(p.x, p.y); });
        glEnd();
    }

    void render_target_view() {
        // 1. Draw Heatmap
        glShadeModel(GL_SMOOTH); glBegin(GL_QUADS);
        int res = 40; double step = 20.0/res;
        for(int i=0; i<res; ++i) for(int j=0; j<res; ++j) {
            double x=-10+i*step, y=-10+j*step;
            Palette::get_heatmap_color(sim.target->prob(Vec2(x,y))).apply(); glVertex2f(x,y);
            Palette::get_heatmap_color(sim.target->prob(Vec2(x+step,y))).apply(); glVertex2f(x+step,y);
            Palette::get_heatmap_color(sim.target->prob(Vec2(x+step,y+step))).apply(); glVertex2f(x+step,y+step);
            Palette::get_heatmap_color(sim.target->prob(Vec2(x,y+step))).apply(); glVertex2f(x,y+step);
        }
        glEnd();

        // 2. Draw Peak Handles
        glPointSize(8.0); glColor3f(1,1,0); glBegin(GL_POINTS);
        for(const auto& p : sim.target->peaks) glVertex2f(p.pos.x, p.pos.y);
        glEnd();
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