#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <functional>

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    Vec3 operator+(const Vec3& other) const { return Vec3(x + other.x, y + other.y, z + other.z); }
    Vec3 operator-(const Vec3& other) const { return Vec3(x - other.x, y - other.y, z - other.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3& operator+=(const Vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    double normSq() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(normSq()); }
};

class P3dGenerator {
public:
    struct Particle { Vec3 pos, vel; };
    std::vector<Particle> particles;
    double x_min, x_max, y_min, y_max, z_min, z_max;

    P3dGenerator(double xm, double xM, double ym, double yM, double zm, double zM)
        : x_min(xm), x_max(xM), y_min(ym), y_max(yM), z_min(zm), z_max(zM) {}

    void generate_uniform(int N) {
        particles.resize(N);
        std::random_device rd; std::mt19937 gen(rd());
        std::uniform_real_distribution<> dX(x_min, x_max);
        std::uniform_real_distribution<> dY(y_min, y_max);
        std::uniform_real_distribution<> dZ(z_min, z_max);
        std::uniform_real_distribution<> dV(-1.0, 1.0);
        for (int i = 0; i < N; ++i) {
            particles[i].pos = Vec3(dX(gen), dY(gen), dZ(gen));
            particles[i].vel = Vec3(dV(gen), dV(gen), dV(gen));
        }
    }

    void swim(double dt) {
        static std::random_device rd; static std::mt19937 gen(rd());
        std::uniform_real_distribution<> noise(-0.1, 0.1);
        for (auto& p : particles) {
            p.pos = p.pos + p.vel * dt;
            p.vel.x += noise(gen); p.vel.y += noise(gen); p.vel.z += noise(gen);
            p.vel = p.vel * 0.99;
            if (p.pos.x < x_min) { p.pos.x=x_min; p.vel.x*=-1; } else if (p.pos.x > x_max) { p.pos.x=x_max; p.vel.x*=-1; }
            if (p.pos.y < y_min) { p.pos.y=y_min; p.vel.y*=-1; } else if (p.pos.y > y_max) { p.pos.y=y_max; p.vel.y*=-1; }
            if (p.pos.z < z_min) { p.pos.z=z_min; p.vel.z*=-1; } else if (p.pos.z > z_max) { p.pos.z=z_max; p.vel.z*=-1; }
        }
    }

    std::vector<Vec3> get_positions() const {
        std::vector<Vec3> pos; pos.reserve(particles.size());
        for(const auto& p : particles) pos.push_back(p.pos);
        return pos;
    }

    void set_positions(const std::vector<Vec3>& new_pos) {
        if(new_pos.size() != particles.size()) return;
        for(size_t i=0; i<particles.size(); ++i) particles[i].pos = new_pos[i];
    }

    void render(const std::function<void(Vec3)>& draw_func) const {
        for (const auto& p : particles) draw_func(p.pos);
    }
};