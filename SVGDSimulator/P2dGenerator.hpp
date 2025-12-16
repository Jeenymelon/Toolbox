#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <functional>

struct Vec2 {
    double x, y;
    Vec2() : x(0), y(0) {}
    Vec2(double _x, double _y) : x(_x), y(_y) {}
    Vec2 operator+(const Vec2& other) const { return Vec2(x + other.x, y + other.y); }
    Vec2 operator-(const Vec2& other) const { return Vec2(x - other.x, y - other.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    Vec2& operator+=(const Vec2& other) { x += other.x; y += other.y; return *this; }
    double normSq() const { return x * x + y * y; }
    double norm() const { return std::sqrt(normSq()); }
};

class P2dGenerator {
public:
    struct Particle { Vec2 pos, vel; };
    std::vector<Particle> particles;
    double x_min, x_max, y_min, y_max;

    P2dGenerator(double xmin, double xmax, double ymin, double ymax)
        : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax) {}

    void generate_uniform(int N) {
        particles.resize(N);
        std::random_device rd; std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(x_min, x_max);
        std::uniform_real_distribution<> disY(y_min, y_max);
        std::uniform_real_distribution<> disVel(-1.0, 1.0);
        for (int i = 0; i < N; ++i) {
            particles[i].pos = Vec2(disX(gen), disY(gen));
            particles[i].vel = Vec2(disVel(gen), disVel(gen));
        }
    }

    void swim(double dt) {
        static std::random_device rd; static std::mt19937 gen(rd());
        std::uniform_real_distribution<> noise(-0.1, 0.1);
        for (auto& p : particles) {
            p.pos = p.pos + p.vel * dt;
            p.vel.x += noise(gen); p.vel.y += noise(gen);
            p.vel = p.vel * 0.99; // Damping
            // Bounce
            if (p.pos.x < x_min) { p.pos.x = x_min; p.vel.x *= -1; }
            if (p.pos.x > x_max) { p.pos.x = x_max; p.vel.x *= -1; }
            if (p.pos.y < y_min) { p.pos.y = y_min; p.vel.y *= -1; }
            if (p.pos.y > y_max) { p.pos.y = y_max; p.vel.y *= -1; }
        }
    }

    std::vector<Vec2> get_positions() const {
        std::vector<Vec2> pos; pos.reserve(particles.size());
        for(const auto& p : particles) pos.push_back(p.pos);
        return pos;
    }

    void set_positions(const std::vector<Vec2>& new_pos) {
        if(new_pos.size() != particles.size()) return;
        for(size_t i=0; i<particles.size(); ++i) particles[i].pos = new_pos[i];
    }

    void render(const std::function<void(Vec2)>& draw_func) const {
        for (const auto& p : particles) draw_func(p.pos);
    }
};