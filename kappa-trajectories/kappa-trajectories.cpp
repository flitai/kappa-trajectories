#include <iostream>
#include <vector>
#include <cmath>    // 数学函数
#include <algorithm> // 用于 std::min, std::max
#include <stdexcept> // 用于运行时错误
#include <string>    // 用于 std::to_string

// 如果cmath中没有定义 M_PI (C++标准不保证，但常见)
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// 二维点/向量结构体
struct Point2D {
    double x, y;

    Point2D(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    Point2D operator+(const Point2D& other) const { return {x + other.x, y + other.y}; }
    Point2D operator-(const Point2D& other) const { return {x - other.x, y - other.y}; }
    Point2D operator*(double scalar) const { return {x * scalar, y * scalar}; }
    Point2D operator/(double scalar) const {
        if (std::abs(scalar) < 1e-12) throw std::runtime_error("二维点标量除法时发生除零错误 (分母过小)。");
        return {x / scalar, y / scalar};
    }
    double dot(const Point2D& other) const { return x * other.x + y * other.y; }
    double norm_sq() const { return x * x + y * y; }
    double norm() const { return std::sqrt(norm_sq()); }
    Point2D normalize() const {
        double n = norm();
        if (n < 1e-9) throw std::runtime_error("无法单位化接近零的向量。");
        return *this / n;
    }
    // 按角度旋转 (弧度)
    Point2D rotate(double angle_rad) const {
        double cos_a = std::cos(angle_rad);
        double sin_a = std::sin(angle_rad);
        return {x * cos_a - y * sin_a, x * sin_a + y * cos_a};
    }
    // 计算两点间距离
    static double distance(const Point2D& p1, const Point2D& p2) {
        return (p1 - p2).norm();
    }
     // 判断两点是否足够接近
    static bool are_close(const Point2D& p1, const Point2D& p2, double tol = 1e-7) {
        return distance(p1, p2) < tol;
    }
};

// 用于存储生成的轨迹部分
struct TurnGeometry {
    Point2D tangent_point_on_incoming_segment; // T1: 入航线段与第一个圆弧的切点
    Point2D tangent_point_on_outgoing_segment; // T2: 最后一个圆弧与出航线段的切点
    std::vector<Point2D> path_points;          // 包含所有圆弧上的点 (从T1开始，到T2结束)
    bool is_valid;
    std::string info; // 用于存储一些调试或状态信息
};

class KappaTrajectoryGenerator {
public:
    /**
     * @brief 在航点间生成平滑的kappa轨迹的圆弧部分。
     *
     * @param w_im1 上一个航点 (w_{i-1})
     * @param w_i 当前航点，即转弯点 (w_i)
     * @param w_ip1 下一个航点 (w_{i+1})
     * @param v_hat 期望速度 (m/s)
     * @param turn_rate_c 最大转向速率 (rad/s)
     * @param kappa kappa参数, 取值范围 [0, 1]
     * @param num_points_per_arc 每个圆弧段生成的点数 (不包括起点，除非是总路径的第一个点)
     * @return TurnGeometry 结构体，包含切点和圆弧上的点
     */
    TurnGeometry generateTurnTrajectory(
        const Point2D& w_im1, const Point2D& w_i, const Point2D& w_ip1,
        double v_hat, double turn_rate_c, double kappa,
        int num_points_per_arc = 20) {

        TurnGeometry result;
        result.is_valid = false;
        const double epsilon = 1e-7; // 用于浮点数比较的小量

        if (turn_rate_c <= epsilon) {
            result.info = "错误: 转向速率过小或为零。";
            // std::cerr << result.info << std::endl;
            return result;
        }
        if (kappa < -epsilon || kappa > 1.0 + epsilon) { // 允许一点点浮点误差
            result.info = "错误: kappa 值 (" + std::to_string(kappa) + ") 必须在 [0, 1] 范围内。";
            // std::cerr << result.info << std::endl;
            return result;
        }
        // 将kappa修正到标准范围，如果只是略微超出
        kappa = std::max(0.0, std::min(1.0, kappa));


        const double R = v_hat / turn_rate_c; // 转弯半径

        // 计算从w_i指向上一个和下一个航点的向量
        Point2D v_from_wi_to_im1 = w_im1 - w_i;
        Point2D v_from_wi_to_ip1 = w_ip1 - w_i;

        if (v_from_wi_to_im1.norm_sq() < epsilon * epsilon || v_from_wi_to_ip1.norm_sq() < epsilon * epsilon) {
            result.info = "错误: 航点过于接近或重合 (" + std::to_string(v_from_wi_to_im1.norm()) + ", " + std::to_string(v_from_wi_to_ip1.norm()) + ")，无法形成有效转弯。";
            // std::cerr << result.info << std::endl;
            return result; 
        }

        Point2D u_from_wi_to_im1 = v_from_wi_to_im1.normalize(); // 指向 w_{i-1}的单位向量
        Point2D u_from_wi_to_ip1 = v_from_wi_to_ip1.normalize(); // 指向 w_{i+1}的单位向量

        double dot_product = u_from_wi_to_im1.dot(u_from_wi_to_ip1);
        dot_product = std::max(-1.0, std::min(1.0, dot_product)); 
        double beta_corner = std::acos(dot_product); // 角 w_{i-1} w_i w_{i+1}

        // 入航线段和出航线段的单位向量 (从前一个航点指向当前航点，从当前航点指向下一个航点)
        Point2D q_in_segment_dir = (w_i - w_im1).normalize();
        Point2D q_out_segment_dir = (w_ip1 - w_i).normalize();

        // 检查是否为直线 (beta_corner接近0度) 或U型转弯 (beta_corner接近PI)
        if (beta_corner < 1e-4) { // 路径几乎是直线 (0度角)
            result.info = "提示: 路径接近直线 (beta_corner=" + std::to_string(beta_corner*180/M_PI) + "度)，无需转弯平滑。";
            result.tangent_point_on_incoming_segment = w_i;
            result.tangent_point_on_outgoing_segment = w_i;
            result.path_points.push_back(w_i);
            result.is_valid = true;
            return result;
        }
         if (std::abs(beta_corner - M_PI) < 1e-4) { // 路径几乎是U型转弯 (180度角)
            // 对于 kappa < 1, 三圆弧模型会退化，中间圆弧的 theta_xi = acos(kappa), p(kappa)=w_i.
            // 对于 kappa = 1, theta_xi = 0, p(kappa)=w_i.
            // 这种情况下，通常期望一个单一的半圆转弯。
            // 此处简化处理：返回w_i作为路径点，表示模型在此处不适用或退化。
            // 一个更完整的实现可能需要为U型转弯提供专门的几何构建方法。
            result.info = "提示: 路径接近U型转弯 (beta_corner=" + std::to_string(beta_corner*180/M_PI) + "度)。kappa轨迹在此定义下可能退化。";
            result.tangent_point_on_incoming_segment = w_i;
            result.tangent_point_on_outgoing_segment = w_i;
            result.path_points.push_back(w_i);
            result.is_valid = true; // 退化情况也认为是有效的，但路径只有一个点
            return result;
        }

        // 角平分线单位向量 (从 w_i 指向转弯外侧/对称轴)
        Point2D u_bisector = (u_from_wi_to_im1 + u_from_wi_to_ip1).normalize();

        // 转弯方向: cross_product > 0 左转 (CCW), < 0 右转 (CW)
        // 使用入航段 (w_im1 -> w_i) 和出航段 (w_i -> w_ip1) 来判断
        double turn_cross_product = q_in_segment_dir.x * q_out_segment_dir.y - q_in_segment_dir.y * q_out_segment_dir.x;
        int turn_direction_sign = (turn_cross_product > 0) ? 1 : -1; // 1 for Left (CCW), -1 for Right (CW)

        double sin_beta_half = std::sin(beta_corner / 2.0);
        if (std::abs(sin_beta_half) < epsilon) { // 避免除零, 此时beta_corner接近0或2PI, 已被前面覆盖
            result.info = "错误: sin(beta_corner/2)过小 ("+ std::to_string(sin_beta_half) +")。beta_corner = " + std::to_string(beta_corner*180/M_PI);
            // std::cerr << result.info << std::endl;
            return result;
        }

        // --- 情况1: kappa = 1 (单一圆弧转弯) ---
        // 对应论文 Fig.3 的内切圆 C_bar, 和 Fig.8 的简化DTS
        if (std::abs(kappa - 1.0) < epsilon) {
            result.info = "Kappa=1: 单圆弧转弯模式。";
            // p(1) 点，即论文中的 p_bar
            Point2D pk_at_1 = w_i + u_bisector * R * (1.0 / sin_beta_half - 1.0);
            // 单一转弯圆弧的圆心 (即论文Fig.3中 C_bar 的圆心)
            Point2D center_single_arc = pk_at_1 - u_bisector * R;

            // 计算入航线段与该圆弧的切点 T1
            Point2D normal_in_to_center = q_in_segment_dir.rotate(turn_direction_sign * M_PI / 2.0);
            result.tangent_point_on_incoming_segment = center_single_arc - normal_in_to_center * R;

            // 计算出航线段与该圆弧的切点 T2
            Point2D normal_out_to_center = q_out_segment_dir.rotate(turn_direction_sign * M_PI / 2.0);
            result.tangent_point_on_outgoing_segment = center_single_arc - normal_out_to_center * R;
            
            addArcPoints(result.path_points, center_single_arc,
                         result.tangent_point_on_incoming_segment,
                         result.tangent_point_on_outgoing_segment,
                         R, turn_direction_sign, num_points_per_arc * 3); // 单圆弧，点数多一些
            result.is_valid = true;
            return result;
        }

        // --- 情况2: kappa < 1 (包括 kappa = 0) (三圆弧模型) ---
        result.info = "Kappa<1 (" + std::to_string(kappa) + "): 三圆弧转弯模式。";
        // 计算 p(kappa) 点
        double dist_wi_pk = kappa * R * (1.0 / sin_beta_half - 1.0);
        Point2D pk = w_i + u_bisector * dist_wi_pk; // p(kappa)
        if (std::abs(kappa) < epsilon) { // kappa = 0, p(kappa) = w_i
            pk = w_i;
            result.info = "Kappa=0: 路径将通过航点 w_i。";
        }

        // 中间圆弧 C_p(kappa) (论文中Fig.9的d点是其圆心) 的圆心 c_d
        Point2D c_d = pk + u_bisector * R;

        // 计算 Xi(kappa, beta_corner) 和 theta_xi (论文中为theta)
        double xi_val = ((1.0 + kappa) + (1.0 - kappa) * sin_beta_half) / 2.0;
        xi_val = std::max(-1.0 + epsilon, std::min(1.0 - epsilon, xi_val)); // 确保 xi_val 在 acos 的有效范围内
        double theta_xi = std::acos(xi_val); // theta_xi 是论文 Eq. (39) 的 theta

        // 中间圆弧的起点 P_A 和终点 P_B (对应Fig.9中的点a, 通过对称性得到另一个)
        // 向量从c_d指向pk是 -u_bisector
        Point2D P_A = c_d + (-u_bisector).rotate(turn_direction_sign * theta_xi) * R;
        Point2D P_B = c_d + (-u_bisector).rotate(-turn_direction_sign * theta_xi) * R;

        // 第一个过渡圆弧 C_i 的圆心 c_1 (对应Fig.9中的点c)
        Point2D u_cd_pa = (P_A - c_d).normalize(); 
        Point2D c_1 = P_A + u_cd_pa * R;

        // 第二个过渡圆弧 C_{i+1} 的圆心 c_2
        Point2D u_cd_pb = (P_B - c_d).normalize(); 
        Point2D c_2 = P_B + u_cd_pb * R;

        // 计算入航线段与C_i的切点 T1
        Point2D normal_in_to_c1 = q_in_segment_dir.rotate(turn_direction_sign * M_PI / 2.0);
        result.tangent_point_on_incoming_segment = c_1 - normal_in_to_c1 * R; // T1

        // 计算出航线段与C_{i+1}的切点 T2
        Point2D normal_out_to_c2 = q_out_segment_dir.rotate(turn_direction_sign * M_PI / 2.0);
        result.tangent_point_on_outgoing_segment = c_2 - normal_out_to_c2 * R; // T2
        
        // --- 生成圆弧上的点 ---
        // 圆弧1: 从 T1 到 P_A, 圆心 c_1
        addArcPoints(result.path_points, c_1, result.tangent_point_on_incoming_segment, P_A, R, turn_direction_sign, num_points_per_arc);
        
        // 圆弧2: 从 P_A 到 P_B, 圆心 c_d
        // 只有当 P_A 和 P_B 不太接近时才添加中间圆弧 (theta_xi 较大时)
        if (!Point2D::are_close(P_A, P_B, R * 1e-3)) { // 如果P_A, P_B不重合
             addArcPoints(result.path_points, c_d, P_A, P_B, R, turn_direction_sign, num_points_per_arc);
        } else {
            //如果P_A, P_B非常接近，确保P_A (或P_B) 在路径中 (如果addArcPoints没有添加它)
            if (result.path_points.empty() || !Point2D::are_close(result.path_points.back(), P_A)) {
                result.path_points.push_back(P_A);
            }
        }


        // 圆弧3: 从 P_B 到 T2, 圆心 c_2
        addArcPoints(result.path_points, c_2, P_B, result.tangent_point_on_outgoing_segment, R, turn_direction_sign, num_points_per_arc);
        
        result.is_valid = true;
        return result;
    }

private:
    // 辅助函数：生成圆弧上的点
    void addArcPoints(std::vector<Point2D>& points_vec,
                      const Point2D& center, const Point2D& start_point, const Point2D& end_point,
                      double radius, int turn_sign, int num_points) {
        
        // 如果起点和终点非常接近，则只添加终点（如果它与vector中最后一个点不同）
        if (Point2D::are_close(start_point, end_point)) {
            if (points_vec.empty() || !Point2D::are_close(points_vec.back(), end_point)) {
                points_vec.push_back(end_point);
            }
            return;
        }

        Point2D start_vec = start_point - center;
        Point2D end_vec = end_point - center;

        double start_angle = std::atan2(start_vec.y, start_vec.x);
        double end_angle = std::atan2(end_vec.y, end_vec.x);
        
        double angle_diff = end_angle - start_angle;

        // 调整angle_diff以匹配转弯方向 (turn_sign: 1 for CCW/Left, -1 for CW/Right)
        if (turn_sign > 0) { // CCW (Left Turn), angle_diff应该是正的
            if (angle_diff < -1e-6) angle_diff += 2.0 * M_PI; // e.g. start_angle=350, end_angle=10 -> diff=-340 -> +20
            else if (angle_diff > 2.0 * M_PI - 1e-6) angle_diff -= 2.0 * M_PI; // e.g. start=10, end=350, but was intended as short way
        } else { // CW (Right Turn), angle_diff应该是负的
            if (angle_diff > 1e-6) angle_diff -= 2.0 * M_PI; // e.g. start_angle=10, end_angle=350 -> diff=340 -> -20
            else if (angle_diff < -2.0 * M_PI + 1e-6) angle_diff += 2.0 * M_PI;
        }

        // 再次确保 angle_diff 和 turn_sign 符号一致性（或 angle_diff 为0）
        // 防止因atan2的周期性导致转错方向或转大圈
        if (std::abs(angle_diff) > 1e-6 && (angle_diff * turn_sign < 0)) {
             angle_diff += turn_sign * 2.0 * M_PI;
        }
        
        // 如果路径点列表为空，或当前圆弧的起点与列表末尾的点不同，则添加起点
        if (points_vec.empty() || !Point2D::are_close(points_vec.back(), start_point)) {
            points_vec.push_back(start_point);
        }

        // 如果调整后的angle_diff非常小，说明起点终点几乎一致，不需要再插值
        if (std::abs(angle_diff) < 1e-5) {
            if (!Point2D::are_close(points_vec.back(), end_point)) { // 仅当终点与最后一个点不同时添加
                 points_vec.push_back(end_point);
            }
            return;
        }

        for (int i = 1; i <= num_points; ++i) {
            double current_angle_rad = start_angle + angle_diff * (static_cast<double>(i) / num_points);
            points_vec.push_back(center + Point2D(std::cos(current_angle_rad), std::sin(current_angle_rad)) * radius);
        }
    }
};


// --- 示例用法 ---
int main() {
    KappaTrajectoryGenerator generator;

    Point2D w0 = {0, 0}; Point2D w1 = {10, 0}; Point2D w2 = {10, 10}; Point2D w3 = {20, 10}; Point2D w4 = {20,0};

    double v_hat = 5.0;    // m/s
    double c_rate = 0.5; // rad/s --> R = 10m
    int points_per_arc = 10;

    std::cout << std::fixed; std::cout.precision(3); // 设置输出精度

    // 测试 kappa = 0.5 (三圆弧)
    double kappa_test = 0.5;
    std::cout << "## 测试 kappa = " << kappa_test << " (w0-w1-w2 右转90度)" << std::endl;
    TurnGeometry turn1 = generator.generateTurnTrajectory(w0, w1, w2, v_hat, c_rate, kappa_test, points_per_arc);
    if (turn1.is_valid) {
        std::cout << turn1.info << std::endl;
        std::cout << "  T1: (" << turn1.tangent_point_on_incoming_segment.x << ", " << turn1.tangent_point_on_incoming_segment.y << ")" << std::endl;
        std::cout << "  T2: (" << turn1.tangent_point_on_outgoing_segment.x << ", " << turn1.tangent_point_on_outgoing_segment.y << ")" << std::endl;
        std::cout << "  路径点数: " << turn1.path_points.size() << std::endl;
        // for(const auto& p : turn1.path_points) std::cout << "    (" << p.x << ", " << p.y << ")" << std::endl;
    } else { std::cout << "  生成失败: " << turn1.info << std::endl;}
    std::cout << "------------------------------------" << std::endl;

    // 测试 kappa = 0 (通过航点 w_i)
    kappa_test = 0.0;
    std::cout << "## 测试 kappa = " << kappa_test << " (w0-w1-w2 右转90度)" << std::endl;
    TurnGeometry turn_k0 = generator.generateTurnTrajectory(w0, w1, w2, v_hat, c_rate, kappa_test, points_per_arc);
    if (turn_k0.is_valid) {
        std::cout << turn_k0.info << std::endl;
        std::cout << "  T1: (" << turn_k0.tangent_point_on_incoming_segment.x << ", " << turn_k0.tangent_point_on_incoming_segment.y << ")" << std::endl;
        std::cout << "  T2: (" << turn_k0.tangent_point_on_outgoing_segment.x << ", " << turn_k0.tangent_point_on_outgoing_segment.y << ")" << std::endl;
        std::cout << "  路径点数: " << turn_k0.path_points.size() << std::endl;
        // 检查中间圆弧是否经过w1 (即p(0))
        bool passes_through_wi = false;
        for(const auto& p : turn_k0.path_points) {
            if (Point2D::distance(p, w1) < 0.1) { // 允许一些误差
                passes_through_wi = true;
                break;
            }
        }
        std::cout << "  路径是否近似通过 w_i (" << w1.x << "," << w1.y << ")? " << (passes_through_wi ? "是" : "否") << std::endl;

    } else { std::cout << "  生成失败: " << turn_k0.info << std::endl;}
    std::cout << "------------------------------------" << std::endl;

    // 测试 kappa = 1 (单一圆弧)
    kappa_test = 1.0;
    std::cout << "## 测试 kappa = " << kappa_test << " (w0-w1-w2 右转90度)" << std::endl;
    TurnGeometry turn_k1 = generator.generateTurnTrajectory(w0, w1, w2, v_hat, c_rate, kappa_test, points_per_arc);
    if (turn_k1.is_valid) {
        std::cout << turn_k1.info << std::endl;
        std::cout << "  T1: (" << turn_k1.tangent_point_on_incoming_segment.x << ", " << turn_k1.tangent_point_on_incoming_segment.y << ")" << std::endl;
        std::cout << "  T2: (" << turn_k1.tangent_point_on_outgoing_segment.x << ", " << turn_k1.tangent_point_on_outgoing_segment.y << ")" << std::endl;
        std::cout << "  路径点数: " << turn_k1.path_points.size() << std::endl;
    } else { std::cout << "  生成失败: " << turn_k1.info << std::endl; }
    std::cout << "------------------------------------" << std::endl;

    // 测试直线情况
    std::cout << "## 测试直线情况 (w1-w2-w3)" << std::endl;
    TurnGeometry turn_straight = generator.generateTurnTrajectory(w1, w2, w3, v_hat, c_rate, 0.5, points_per_arc);
    if (turn_straight.is_valid) {
        std::cout << turn_straight.info << std::endl;
        std::cout << "  T1: (" << turn_straight.tangent_point_on_incoming_segment.x << ", " << turn_straight.tangent_point_on_incoming_segment.y << ")" << std::endl;
        std::cout << "  T2: (" << turn_straight.tangent_point_on_outgoing_segment.x << ", " << turn_straight.tangent_point_on_outgoing_segment.y << ")" << std::endl;
        std::cout << "  路径点数: " << turn_straight.path_points.size() << std::endl;
    } else { std::cout << "  生成失败: " << turn_straight.info << std::endl; }
    std::cout << "------------------------------------" << std::endl;

    // 测试U型转弯
    std::cout << "## 测试U型转弯 (w2-w3-w4)" << std::endl; // (10,10) -> (20,10) -> (20,0)
    TurnGeometry turn_uturn = generator.generateTurnTrajectory(w2, w3, w4, v_hat, c_rate, 0.5, points_per_arc);
     if (turn_uturn.is_valid) {
        std::cout << turn_uturn.info << std::endl;
        std::cout << "  T1: (" << turn_uturn.tangent_point_on_incoming_segment.x << ", " << turn_uturn.tangent_point_on_incoming_segment.y << ")" << std::endl;
        std::cout << "  T2: (" << turn_uturn.tangent_point_on_outgoing_segment.x << ", " << turn_uturn.tangent_point_on_outgoing_segment.y << ")" << std::endl;
        std::cout << "  路径点数: " << turn_uturn.path_points.size() << std::endl;
    } else { std::cout << "  生成失败: " << turn_uturn.info << std::endl; }
    std::cout << "------------------------------------" << std::endl;


    Point2D wa = {0,0}, wb={10,0}, wc={5, 8.66}; // 左转60度 beta_corner=120度
    std::cout << "## 测试左转120度 (beta_corner=120)" << std::endl;
    TurnGeometry turn_left = generator.generateTurnTrajectory(wa, wb, wc, v_hat, c_rate, 0.3, points_per_arc);
    if (turn_left.is_valid) {
        std::cout << turn_left.info << std::endl;
        std::cout << "  T1: (" << turn_left.tangent_point_on_incoming_segment.x << ", " << turn_left.tangent_point_on_incoming_segment.y << ")" << std::endl;
        std::cout << "  T2: (" << turn_left.tangent_point_on_outgoing_segment.x << ", " << turn_left.tangent_point_on_outgoing_segment.y << ")" << std::endl;
        std::cout << "  路径点数: " << turn_left.path_points.size() << std::endl;
    } else { std::cout << "  生成失败: " << turn_left.info << std::endl;}


    return 0;
}