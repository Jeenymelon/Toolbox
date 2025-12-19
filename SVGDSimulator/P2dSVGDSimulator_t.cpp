#include <GLFW/glfw3.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>

struct Vec2 {
    double x, y;
    Vec2(double _x = 0, double _y = 0) : x(_x), y(_y) {}
    Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    double distSq(Vec2 v) const { return (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y); }
    void clamp(double b) { x = std::max(-b, std::min(b, x)); y = std::max(-b, std::min(b, y)); }
};

class SVGDApp
{
    GLFWwindow* window;
    std::vector<Vec2> particles, peaks, sample_pool;

    bool is_median_mode = true;
    double sigma = 0.6;          // 紫色条: 采样噪声 (UP/DN)
    double step_size = 0.06;     // 红色条: 仿真步长 (W/S)
    double alpha = 1.0;          // 蓝色条: 带宽倍率 (A/D)
    double h_eff = 1.0, h_base = 1.0;

    int p_count = 200;           // 粒子数量 (O/P)
    int memory_depth = 150;      // 绿色条: 记忆深度 (N/M)
    float grid[81][81];

    // 交互状态
    int dragged_idx = -1;

public:
    SVGDApp() {
        if (!glfwInit()) exit(-1);
        window = glfwCreateWindow(1200, 800, "SVGD Perfect Layout Lab", NULL, NULL);
        glfwMakeContextCurrent(window);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        peaks.push_back({0, 0});
        reset_particles();
        print_ui();
    }

    void reset_particles() {
        particles.clear();
        for(int i=0; i<p_count; i++)
            particles.push_back({(double)(rand()%16-8), (double)(rand()%16-8)});
    }

    void print_ui() {
#ifdef _WIN32
        system("cls");
#else
        std::cout << "\033[2J\033[1;1H";
#endif
        std::cout << "========================================================\n";
        std::cout << "             SVGD STRATEGY LAB - UI OPTIMIZED           \n";
        std::cout << "========================================================\n";
        std::cout << "[MODE]: " << (is_median_mode ? "AUTO-MEDIAN" : "FIXED-MANUAL") << "\n";
        std::cout << "--------------------------------------------------------\n";
        std::cout << "[TWEAK PANEL MAPPING]:\n";
        std::cout << "  - Blue  (A/D)  : Bandwidth Alpha  -> " << alpha << "\n";
        std::cout << "  - Red   (W/S)  : Sim Step Size    -> " << step_size << "\n";
        std::cout << "  - Green (N/M)  : Memory Depth     -> " << memory_depth << "\n";
        std::cout << "  - Purple(UP/DN): Sampling Sigma   -> " << sigma << "\n";
        std::cout << "  - Yellow(READ) : Effective H      -> " << h_eff << "\n";
        std::cout << "--------------------------------------------------------\n";
        std::cout << " [TAB]: Toggle Mode | [O/P]: Particles | [L-Click]: Peak\n";
        std::cout << "========================================================\n" << std::endl;
    }

    void handle_input() {
        static bool l_lock = false, r_lock = false, m_lock = false;
        double mx, my; glfwGetCursorPos(window, &mx, &my);
        int w, h; glfwGetWindowSize(window, &w, &h);

        // 坐标映射：屏幕像素 -> 仿真空间 [-10, 10]
        Vec2 m_pos(((mx - 200.0) / (w - 200.0)) * 20.0 - 10.0, -((my / h) * 20.0 - 10.0));
        m_pos.clamp(9.8);

        bool in_canvas = (mx > 200);

        // --- 鼠标交互逻辑 ---

        // 左键：添加
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            if (!l_lock && in_canvas) {
                peaks.push_back(m_pos);
                l_lock = true;
            }
        } else l_lock = false;

        // 右键：消除最近的点
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
            if (!r_lock && in_canvas && !peaks.empty()) {
                int best_idx = -1; double min_d = 2.0; // 交互半径
                for(int i=0; i<peaks.size(); i++) {
                    double d = m_pos.distSq(peaks[i]);
                    if(d < min_d) { min_d = d; best_idx = i; }
                }
                if(best_idx != -1) peaks.erase(peaks.begin() + best_idx);
                r_lock = true;
            }
        } else r_lock = false;

        // 中键：移动拖拽
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS) {
            if (in_canvas) {
                if (!m_lock) { // 初始按下，寻找最近点
                    double min_d = 2.0;
                    for(int i=0; i<peaks.size(); i++) {
                        double d = m_pos.distSq(peaks[i]);
                        if(d < min_d) { min_d = d; dragged_idx = i; }
                    }
                    m_lock = true;
                }
                if (dragged_idx != -1) {
                    peaks[dragged_idx] = m_pos; // 持续跟随
                }
            }
        } else {
            m_lock = false;
            dragged_idx = -1;
        }

        // --- 键盘参数控制 ---
        if(glfwGetKey(window, GLFW_KEY_TAB)==GLFW_PRESS) { /* 模式切换逻辑...略，同前 */ }
        if(glfwGetKey(window, GLFW_KEY_D)==GLFW_PRESS) alpha += 0.02;
        if(glfwGetKey(window, GLFW_KEY_A)==GLFW_PRESS) alpha = std::max(0.01, alpha - 0.02);
        if(glfwGetKey(window, GLFW_KEY_W)==GLFW_PRESS) step_size += 0.002;
        if(glfwGetKey(window, GLFW_KEY_S)==GLFW_PRESS) step_size = std::max(0.001, step_size - 0.002);
    }

    void update() {
        static std::mt19937 gen(42); std::normal_distribution<> nd(0, sigma);
        for(int i=0; i<10; i++){
            if(peaks.empty())break;
            Vec2 s = peaks[rand()%peaks.size()] + Vec2(nd(gen), nd(gen));
            if(sample_pool.size() < (size_t)memory_depth) sample_pool.push_back(s);
            else { static int r=0; sample_pool[r%memory_depth]=s; r++; }
        }
        if(sample_pool.size() > (size_t)memory_depth) sample_pool.resize(memory_depth);

        size_t n = particles.size();
        if(is_median_mode && n > 1) {
            std::vector<double> d; int lim=std::min((int)n,30);
            for(int i=0;i<lim;i++) for(int j=i+1;j<lim;j++) d.push_back(particles[i].distSq(particles[j]));
            std::sort(d.begin(), d.end());
            h_base = (d.empty()||d[d.size()/2]<1e-4) ? 0.5 : d[d.size()/2]/std::log(n+1.0);
            h_eff = alpha * h_base;
        } else h_eff = alpha;

        std::vector<Vec2> phi(n, Vec2(0,0));
        for(size_t i=0; i<n; i++) {
            Vec2 grad_logp(0,0);
            for(auto& s : sample_pool) grad_logp = grad_logp + (s - particles[i]) * (1.0/(sigma*sigma+0.1));
            grad_logp = grad_logp * (1.0/std::max(1.0,(double)sample_pool.size())) + particles[i] * -0.015;
            Vec2 force(0,0);
            for(size_t j=0; j<n; j++) {
                Vec2 diff = particles[i] - particles[j];
                double k = std::exp(-diff.distSq({0,0}) / h_eff);
                force = force + grad_logp * k + diff * (k * 2.0 / h_eff);
            }
            phi[i] = force * (1.0/n);
        }
        for(size_t i=0; i<n; i++) { particles[i] = particles[i] + phi[i]*step_size; particles[i].clamp(9.8); }
    }

    void draw() {
        int w, h; glfwGetFramebufferSize(window, &w, &h);
        glClearColor(0.01, 0.01, 0.04, 1.0); glClear(GL_COLOR_BUFFER_BIT);

        // --- 1. 左侧 Tweak Bars (保持不变) ---
        glViewport(0, 0, (int)(200.0 * w / 1200.0), h);
        glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(0, 200, 0, 800, -1, 1);
        glColor3f(0.05, 0.05, 0.1); glBegin(GL_QUADS); glVertex2f(0, 0); glVertex2f(200, 0); glVertex2f(200, 800); glVertex2f(0, 800); glEnd();

        auto draw_bar = [&](int i, float val, float max_v, float r, float g, float b) {
            float y_base = 720 - i * 65;
            glColor3f(0.15, 0.15, 0.25);
            glBegin(GL_QUADS); glVertex2f(20, y_base); glVertex2f(180, y_base); glVertex2f(180, y_base + 10); glVertex2f(20, y_base + 10); glEnd();
            float fill = std::min(1.0f, val / max_v);
            glColor3f(r, g, b);
            glBegin(GL_QUADS); glVertex2f(20, y_base); glVertex2f(20 + fill * 160, y_base); glVertex2f(20 + fill * 160, y_base + 10); glVertex2f(20, y_base + 10); glEnd();
        };
        draw_bar(0, (float)alpha, 4.0f, 0.2, 0.6, 1.0);
        draw_bar(1, (float)step_size, 0.25f, 0.8, 0.2, 0.2);
        draw_bar(2, (float)memory_depth, 1500.0f, 0.2, 0.8, 0.4);
        draw_bar(3, (float)sigma, 2.0f, 0.6, 0.3, 0.9);
        draw_bar(4, (float)h_eff, 8.0f, 0.9, 0.8, 0.2);

        // --- 2. Main View (优化热力图) ---
        glViewport((int)(200.0 * w / 1200.0), 0, (int)(1000.0 * w / 1200.0), h);
        glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(-10, 10, -10, 10, -1, 1);

        // 计算高精度网格密度 (80x80)
        for (int i = 0; i < 80; i++) for (int j = 0; j < 80; j++) grid[i][j] = 0;
        for (auto& p : particles) {
            int gx = (int)((p.x + 10) / 20 * 80), gy = (int)((p.y + 10) / 20 * 80);
            if (gx >= 0 && gx < 80 && gy >= 0 && gy < 80) grid[gx][gy] += 1.2f;
        }

        // 绘制热力图格子
        float step = 20.0f / 80.0f;
        for (int i = 0; i < 80; i++) {
            for (int j = 0; j < 80; j++) {
                float v = grid[i][j];
                if (v < 0.05f) continue;

                float x = i * step - 10, y = j * step - 10;
                float alpha_val = std::min(0.6f, v * 0.15f); // 降低透明度

                // 填充色：根据密度平滑过渡
                glBegin(GL_QUADS);
                if (v < 2.0f) glColor4f(0.1, 0.3, 0.6, alpha_val);
                else if (v < 5.0f) glColor4f(0.2, 0.6, 0.5, alpha_val);
                else glColor4f(0.9, 0.4, 0.1, alpha_val);

                glVertex2f(x, y); glVertex2f(x + step, y);
                glVertex2f(x + step, y + step); glVertex2f(x, y + step);
                glEnd();

                // 强化边缘：绘制格子线（仅在有密度的区域）
                if (v > 0.5f) {
                    glLineWidth(0.5f);
                    glColor4f(1.0, 1.0, 1.0, 0.1f); // 极细的白色半透描边
                    glBegin(GL_LINE_LOOP);
                    glVertex2f(x, y); glVertex2f(x + step, y);
                    glVertex2f(x + step, y + step); glVertex2f(x, y + step);
                    glEnd();
                }
            }
        }

        // 采样云 (青色点)
        glPointSize(3.0f); glBegin(GL_POINTS);
        for(auto& s : sample_pool) { glColor4f(0.0, 0.9, 0.9, 0.15); glVertex2f(s.x, s.y); }
        glEnd();

        // 粒子 (白色点 - 核心焦点)
        glPointSize(4.0f); glColor4f(1.0, 1.0, 1.0, 0.9); glBegin(GL_POINTS);
        for(auto& p : particles) glVertex2f(p.x, p.y);
        glEnd();

        // 目标红十字
        glLineWidth(2.0f);
        for(auto& pk:peaks){
            glColor4f(1.0, 0.2, 0.2, 0.8);
            glBegin(GL_LINES);
            glVertex2f(pk.x-0.3, pk.y); glVertex2f(pk.x+0.3, pk.y);
            glVertex2f(pk.x, pk.y-0.3); glVertex2f(pk.x, pk.y+0.3);
            glEnd();
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    void run() { while(!glfwWindowShouldClose(window)){ handle_input(); update(); draw(); } }
};

int main() { SVGDApp().run(); return 0; }