#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/xfeatures2d/nonfree.hpp>

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

// --- 配置常量 ---
const string VIDEO_PATH = "../llama.mp4";
const string FILE_PATH = "../points.yml";
const int MAX_TOTAL_FEATURES = 500;
const int GRID_ROWS = 6; // 网格行数
const int GRID_COLS = 6; // 网格列数
const Scalar KEYPOINT_COLOR = Scalar(0, 255, 255); // 黄色 (BGR)
const Scalar GRID_COLOR = Scalar(255, 0, 0);       // 蓝色 (BGR)
const Scalar TEXT_COLOR = Scalar(255, 0, 255);     // 紫色 (BGR)
// --------------------

void write_keypoints_all_features_to_text(const string& filename, const vector<KeyPoint>& keypoints) {
    FileStorage fs(FILE_PATH, FileStorage::APPEND);
    write(fs, "keypoints", keypoints);
    fs.release();
}

void test_orb_standard() {
     // 1. Initialize VideoCapture and Check
    VideoCapture cap(VIDEO_PATH);

    if (!cap.isOpened()) {
        cerr << "FATAL: Could not open video file at: " << VIDEO_PATH << endl;
        return;
    }

    // 2. Initialize the ORB Detector and Descriptor
    // This resolves the patent error (no need for xfeatures2d::SURF)
    Ptr<ORB> detector = ORB::create(MAX_TOTAL_FEATURES);

    cout << "*** ORB Feature Extractor Initialized (Patent-Free) ***" << endl;
    cout << "Max Features: " << MAX_TOTAL_FEATURES << endl;

    Mat img;
    int frame_counter = 0;

    // 3. Main Video Processing Loop
    while (true) {
        // Read a new frame from the video stream
        cap >> img;

        // Exit loop if frame read fails
        if (img.empty()) {
            cout << "INFO: End of video reached." << endl;
            break;
        }

        frame_counter++;

        // --- Feature Extraction (Detection & Computation) ---

        // Feature detection often performs better on grayscale images
        Mat gray_img;
        cvtColor(img, gray_img, COLOR_BGR2GRAY);

        vector<KeyPoint> keypoints;
        Mat descriptors;

        // Detect keypoints and compute binary descriptors in one call
        detector->detectAndCompute(gray_img, noArray(), keypoints, descriptors);

        // --- Render Features ---
        Mat img_keypoints;

        // Draw ORB features onto the colored image
        drawKeypoints(img,
                      keypoints,
                      img_keypoints,
                      KEYPOINT_COLOR,
                      // Use DEFAULT flag for ORB; it does not draw orientation like RICH
                      DrawMatchesFlags::DEFAULT);

        // Overlay statistics
        putText(img_keypoints,
                "ORB Keypoints: " + to_string(keypoints.size()),
                Point(10, 30),
                FONT_HERSHEY_SIMPLEX,
                0.8,
                TEXT_COLOR, // Yellow text for contrast
                2);

        // --- write down ---
        write_keypoints_all_features_to_text(FILE_PATH, keypoints);

        // --- Display ---
        imshow("ORB Feature Extraction", img_keypoints);

        // Exit on 'q' press
        if (waitKey(1) == 'q') {
            cout << "INFO: User termination requested." << endl;
            break;
        }
    }

    // Clean up
    destroyAllWindows();
    cout << "Processed " << frame_counter << " total frames. Execution finished." << endl;
}

void test_orb_in_grid() {
    // 1. 初始化视频捕获并检查
    VideoCapture cap(VIDEO_PATH);

    if (!cap.isOpened()) {
        cerr << "FATAL: 无法打开视频文件: " << VIDEO_PATH << endl;
        return;
    }

    // 2. 初始化核心参数
    cout << "*** ORB 网格特征提取器启动 (" << GRID_ROWS << "x" << GRID_COLS << " 网格) ***" << endl;

    // 确定每个网格允许的最大特征数
    const int features_per_grid = MAX_TOTAL_FEATURES / (GRID_ROWS * GRID_COLS);

    Mat frame;
    int frame_counter = 0;

    // 3. 视频处理主循环
    while (true) {
        cap >> frame;
        if (frame.empty()) {
            cout << "INFO: 视频播放结束。" << endl;
            break;
        }

        frame_counter++;

        // 提取灰度图
        Mat gray_img;
        cvtColor(frame, gray_img, COLOR_BGR2GRAY);

        vector<KeyPoint> all_keypoints;
        Mat all_descriptors;

        // 计算单个网格单元的尺寸
        int grid_width = frame.cols / GRID_COLS;
        int grid_height = frame.rows / GRID_ROWS;

        // --- 4. 核心网格划分与局部检测循环 ---
        for (int r = 0; r < GRID_ROWS; ++r) {
            for (int c = 0; c < GRID_COLS; ++c) {

                // 定义当前网格单元的 ROI (Region of Interest)
                Rect roi(c * grid_width,
                         r * grid_height,
                         grid_width,
                         grid_height);

                // 确保最后一个单元覆盖可能存在的剩余像素
                if (c == GRID_COLS - 1) roi.width = frame.cols - roi.x;
                if (r == GRID_ROWS - 1) roi.height = frame.rows - roi.y;

                Mat roi_img = gray_img(roi);

                // 局部 ORB 检测器
                // 确保局部检测器使用相同的尺度和层级参数，以保证描述子的兼容性
                Ptr<ORB> local_detector = ORB::create(features_per_grid,
                                                       1.2f, // scaleFactor
                                                       8);   // nlevels

                vector<KeyPoint> local_keypoints;
                Mat local_descriptors;

                // 在 ROI 上运行 ORB
                local_detector->detectAndCompute(roi_img, noArray(), local_keypoints, local_descriptors);

                // 5. 坐标校正与合并

                // 坐标校正：将局部坐标平移到全局坐标系
                for (auto& kp : local_keypoints) {
                    kp.pt.x += roi.x;
                    kp.pt.y += roi.y;
                    all_keypoints.push_back(kp);
                }

                // 描述子合并
                if (!local_descriptors.empty()) {
                    if (all_descriptors.empty()) {
                        all_descriptors = local_descriptors;
                    } else {
                        vconcat(all_descriptors, local_descriptors, all_descriptors); // 垂直拼接描述子矩阵
                    }
                }
            }
        }

        // --- 6. 渲染输出 ---
        Mat img_output = frame.clone();

        // --- write down ---
        write_keypoints_all_features_to_text(FILE_PATH, all_keypoints);

        // 绘制关键点
        drawKeypoints(img_output,
                      all_keypoints,
                      img_output,
                      KEYPOINT_COLOR,
                      DrawMatchesFlags::DEFAULT);

        // 绘制网格线
        for (int i = 1; i < GRID_COLS; ++i) {
            line(img_output, Point(i * grid_width, 0), Point(i * grid_width, frame.rows), GRID_COLOR, 1);
        }
        for (int i = 1; i < GRID_ROWS; ++i) {
            line(img_output, Point(0, i * grid_height), Point(frame.cols, i * grid_height), GRID_COLOR, 1);
        }

        // 统计信息
        putText(img_output,
                "Grid KPs: " + to_string(all_keypoints.size()) + " / Descriptors: " + to_string(all_descriptors.rows),
                Point(10, 30),
                FONT_HERSHEY_SIMPLEX,
                0.8,
                TEXT_COLOR,
                2);

        // --- 7. 显示与控制 ---
        imshow("ORB 网格特征提取 (test_orb_in_grid)", img_output);

        if (waitKey(1) == 'q') {
            cout << "INFO: 用户请求终止。" << endl;
            break;
        }
    }

    // 8. 清理资源
    destroyAllWindows();
    cout << "总共处理了 " << frame_counter << " 帧。程序结束。" << endl;
}

void test_sift_standard() {
    VideoCapture cap(VIDEO_PATH);

    if (!cap.isOpened()) {
        cerr << "FATAL: 无法打开视频文件: " << VIDEO_PATH << endl;
        return;
    }

    // 1. 初始化 SIFT 检测器
    // nfeatures: 要保留的最佳特征数量
    // nOctaveLayers: 每层金字塔的层数
    Ptr<SIFT> detector = SIFT::create(MAX_TOTAL_FEATURES, 3, 0.04, 10, 1.6);

    cout << "*** SIFT 标准特征提取器启动 ***" << endl;
    cout << "目标特征数: " << MAX_TOTAL_FEATURES << endl;

    Mat frame;
    int frame_counter = 0;

    // 2. 视频处理主循环
    while (true) {
        cap >> frame;
        if (frame.empty()) {
            cout << "INFO: 视频播放结束。" << endl;
            break;
        }

        frame_counter++;

        // 提取灰度图
        Mat gray_img;
        cvtColor(frame, gray_img, COLOR_BGR2GRAY);

        vector<KeyPoint> keypoints;
        Mat descriptors;

        // 3. 核心 SIFT 操作：检测关键点并计算描述子
        detector->detectAndCompute(gray_img, noArray(), keypoints, descriptors);

        // --- 4. 渲染输出 ---
        Mat img_output = frame.clone();

        // --- write down ---
        write_keypoints_all_features_to_text(FILE_PATH, keypoints);

        // 绘制关键点 (使用 DRAW_RICH_KEYPOINTS 来显示方向和大小)
        drawKeypoints(img_output,
                      keypoints,
                      img_output,
                      KEYPOINT_COLOR,
                      DrawMatchesFlags::DRAW_RICH_KEYPOINTS);

        // 统计信息
        putText(img_output,
                "SIFT KPs: " + to_string(keypoints.size()) + " / Descriptors: " + to_string(descriptors.rows),
                Point(10, 30),
                FONT_HERSHEY_SIMPLEX,
                0.8,
                TEXT_COLOR, // 青色文字
                2);

        // --- 5. 显示与控制 ---
        imshow("SIFT 标准特征提取", img_output);

        if (waitKey(1) == 'q') {
            cout << "INFO: 用户请求终止。" << endl;
            break;
        }
    }

    // 6. 清理资源
    destroyAllWindows();
    cout << "总共处理了 " << frame_counter << " 帧。程序结束。" << endl;

}

void test_sift_in_grid() {
    // 1. 初始化视频捕获并检查
    VideoCapture cap(VIDEO_PATH);

    if (!cap.isOpened()) {
        cerr << "FATAL: 无法打开视频文件: " << VIDEO_PATH << endl;
        return;
    }

    // 2. 初始化核心参数
    cout << "*** SIFT 网格特征提取器启动 (" << GRID_ROWS << "x" << GRID_COLS << " 网格) ***" << endl;

    // 确定每个网格允许的最大特征数
    const int features_per_grid = MAX_TOTAL_FEATURES / (GRID_ROWS * GRID_COLS);

    Mat frame;
    int frame_counter = 0;

    // 3. 视频处理主循环
    while (true) {
        cap >> frame;
        if (frame.empty()) {
            cout << "INFO: 视频播放结束。" << endl;
            break;
        }

        frame_counter++;

        // 提取灰度图
        Mat gray_img;
        cvtColor(frame, gray_img, COLOR_BGR2GRAY);

        vector<KeyPoint> all_keypoints;
        Mat all_descriptors;

        // 计算单个网格单元的尺寸
        int grid_width = frame.cols / GRID_COLS;
        int grid_height = frame.rows / GRID_ROWS;

        // --- 4. 核心网格划分与局部检测循环 ---
        for (int r = 0; r < GRID_ROWS; ++r) {
            for (int c = 0; c < GRID_COLS; ++c) {

                // 定义当前网格单元的 ROI (Region of Interest)
                Rect roi(c * grid_width,
                         r * grid_height,
                         grid_width,
                         grid_height);

                // 确保最后一个单元覆盖可能存在的剩余像素
                if (c == GRID_COLS - 1) roi.width = frame.cols - roi.x;
                if (r == GRID_ROWS - 1) roi.height = frame.rows - roi.y;

                Mat roi_img = gray_img(roi);

                // 局部 SIFT 检测器 (设置特征数量限制为 features_per_grid)
                Ptr<SIFT> local_detector = SIFT::create(features_per_grid);

                vector<KeyPoint> local_keypoints;
                Mat local_descriptors;

                // 在 ROI 上运行 SIFT
                local_detector->detectAndCompute(roi_img, noArray(), local_keypoints, local_descriptors);

                // 5. 坐标校正与合并

                // 坐标校正：将局部坐标平移到全局坐标系
                for (auto& kp : local_keypoints) {
                    kp.pt.x += roi.x;
                    kp.pt.y += roi.y;
                    all_keypoints.push_back(kp);
                }

                // 描述子合并
                if (!local_descriptors.empty()) {
                    if (all_descriptors.empty()) {
                        all_descriptors = local_descriptors;
                    } else {
                        // SIFT 描述子是 CV_32F (浮点数) 类型，使用 vconcat 垂直拼接
                        vconcat(all_descriptors, local_descriptors, all_descriptors);
                    }
                }
            }
        }

        // --- 6. 渲染输出 ---
        Mat img_output = frame.clone();

        // --- write down ---
        write_keypoints_all_features_to_text(FILE_PATH, all_keypoints);

        // 绘制关键点 (使用 DRAW_RICH_KEYPOINTS 来显示方向和大小)
        drawKeypoints(img_output,
                      all_keypoints,
                      img_output,
                      KEYPOINT_COLOR,
                      DrawMatchesFlags::DRAW_RICH_KEYPOINTS);

        // 绘制网格线
        for (int i = 1; i < GRID_COLS; ++i) {
            line(img_output, Point(i * grid_width, 0), Point(i * grid_width, frame.rows), GRID_COLOR, 1);
        }
        for (int i = 1; i < GRID_ROWS; ++i) {
            line(img_output, Point(0, i * grid_height), Point(frame.cols, i * grid_height), GRID_COLOR, 1);
        }

        // 统计信息
        putText(img_output,
                "Grid SIFT KPs: " + to_string(all_keypoints.size()) + " / Descriptors: " + to_string(all_descriptors.rows),
                Point(10, 30),
                FONT_HERSHEY_SIMPLEX,
                0.8,
                TEXT_COLOR,
                2);

        // --- 7. 显示与控制 ---
        imshow("SIFT 网格特征提取", img_output);

        if (waitKey(1) == 'q') {
            cout << "INFO: 用户请求终止。" << endl;
            break;
        }
    }

    // 8. 清理资源
    destroyAllWindows();
    cout << "总共处理了 " << frame_counter << " 帧。程序结束。" << endl;
}

int main() {
    // test_orb_standard();

    // test_orb_in_grid();

    // test_sift_standard();

    // test_sift_in_grid();

    return 0;
}