#include <iostream>
#include <vector>
#include <cmath>    // 数学函数
#include <numeric>  // 未来可能用于向量操作，暂时未直接使用
#include <stdexcept> // 用于运行时错误
#include <algorithm> // 用于 std::min 和 std::max
#include <string>    // 用于 std::to_string

// 定义一个简单的二维点结构 (可以替换为您项目中已有的 Point/Vector 类)
struct Point2D {
    double x, y;

    Point2D operator-(const Point2D& other) const {
        return {x - other.x, y - other.y};
    }

    Point2D operator+(const Point2D& other) const {
        return {x + other.x, y + other.y};
    }

    Point2D operator*(double scalar) const {
        return {x * scalar, y * scalar};
    }

    Point2D operator/(double scalar) const {
        if (scalar == 0) throw std::runtime_error("二维点标量除法时发生除零错误。");
        return {x / scalar, y / scalar};
    }

    // 计算向量的模（长度）
    double norm() const {
        return std::sqrt(x*x + y*y);
    }

    // 单位化向量
    Point2D normalize() const {
        double n = norm();
        if (n == 0) throw std::runtime_error("无法单位化零向量。");
        return {x / n, y / n};
    }

    // 计算两个向量的点积
    static double dot(const Point2D& v1, const Point2D& v2) {
        return v1.x * v2.x + v1.y * v2.y;
    }
};

class KappaTrajectorySolver {
public:
    // 构造函数
    KappaTrajectorySolver() = default;

    /**
     * @brief 计算最优的 kappa* 值。
     *
     * @param w_prev 上一个航点。
     * @param w_i 当前航点 (转弯点)。
     * @param w_next 下一个航点。
     * @param v_hat 预设速度。
     * @param c 转向速率 (rad/s)。
     * @param eps kappa* 的目标精度。
     * @return 最优的 kappa* 值。
     */
    double findKappaStar(const Point2D& w_prev, const Point2D& w_i, const Point2D& w_next,
                           double v_hat, double c, double eps = 1e-6) {
        if (c == 0) {
            throw std::invalid_argument("转向速率 'c' 不能为零。");
        }
        R_ = v_hat / c; // 转弯半径

        Point2D v1 = w_i - w_prev; // 从 w_prev 指向 w_i 的向量
        Point2D v2 = w_next - w_i; // 从 w_i 指向 w_next 的向量

        if (v1.norm() == 0 || v2.norm() == 0) {
            // 如果任意一段路径长度为零，则无法定义转弯
            // 这意味着航点重合，beta_ 角未定义
            throw std::invalid_argument("航点导致零长度路径段，无法计算转弯角。");
        }

        Point2D t_hat_im1_i = v1.normalize(); // 前一段路径的单位方向向量
        Point2D t_hat_i_ip1 = v2.normalize(); // 后一段路径的单位方向向量

        double dot_product = Point2D::dot(t_hat_i_ip1, t_hat_im1_i);
        // 将点积限制在 [-1, 1] 范围内，以避免 acos 的定义域错误 (由于浮点数不精确)
        dot_product = std::max(-1.0, std::min(1.0, dot_product));
        beta_ = std::acos(dot_product); // 两段路径间的转角 (0 到 PI)
                                       // beta_ = 0 表示直线, beta_ = PI 表示180度掉头

        // 处理 beta_ 接近 0 (直线) 或 PI (U型转弯) 的情况
        // 文档中的 Lambda(kappa) 公式在 beta=0 时包含 cot(beta/2)，会导致未定义
        // 算法前提是 Lambda(0)>0 和 Lambda(1)<0，这适用于一个明确的转弯
        const double angle_tolerance = 1e-7; // 用于比较角度的容差
        if (beta_ < angle_tolerance) { // 路径几乎是直线
            // std::cerr << "提示: beta 非常接近 0 (路径几乎是直线)。 Lambda 函数的特定形式可能不适用。"
            //          << " 对于直线，L(kappa) = L_orig，Lambda(kappa) 理论上恒为 0。"
            //          << " 返回 kappa = 0.0 (或 1.0) 作为约定。" << std::endl;
            // 对于直线，长度本身就匹配，不需要kappa调整，或者说任何kappa都行。
            // 但二分法依赖 Lambda(0)>0, Lambda(1)<0。
            // 如果是直线，Lambda 应该总是0。
            // 为了DTS系统有一个值，可以返回0或1，或根据系统要求处理。
            // 这里返回0，因为此时不需要“缩短”转弯来匹配长度。
            return 0.0;
        }
        if (std::abs(beta_ - M_PI) < angle_tolerance) { // 路径几乎是U型转弯
            // std::cerr << "提示: beta 非常接近 PI (U型转弯)。" << std::endl;
            // 对于理想的U型转弯，公式也可能达到某些极限情况，但通常是明确定义的。
            // Xi 对于 kappa=1 且 sin(beta/2)=1 (即 beta=PI) 时为 1。
            // 此时 2*acos(Xi) = 0, Xi*sqrt(1/Xi^2-1) = 0
            // (PI-beta)/2 = 0
            // -(1-kappa)*cos(beta/2) = -(1-kappa)*0 = 0
            // -kappa*cot(beta/2) = -kappa*0 = 0
            // 所以 Lambda(kappa) = 0 for beta=PI. 任何 kappa 都可以。
            // 类似于直线情况，返回一个约定值。
            return 0.0; // 或者 1.0
        }


        double a = 0.0;
        double b = 1.0;

        double lambda_a = calculateLambda(a);
        double lambda_b = calculateLambda(b);

        // 根据文档 Lemma 7: Lambda(0) > 0 且 Lambda(1) < 0
        // 为浮点数比较增加一个小容差
        double check_eps = 1e-9; // 用于检查符号条件的容差
        if (lambda_a < -check_eps || lambda_b > check_eps) {
            // 如果端点条件不满足，二分法无法按预期进行。
            // 这可能表示输入几何形状的问题，或者算法假设不适用于此特定情况。
            std::string error_msg = "Lambda 函数不满足端点条件 Lambda(0)>0 和 Lambda(1)<0。\n";
            error_msg += "  请检查航点几何形状或此情况下的算法假设。\n";
            error_msg += "  beta = " + std::to_string(beta_ * 180.0 / M_PI) + " 度\n";
            error_msg += "  Lambda(0) = " + std::to_string(lambda_a) + "\n";
            error_msg += "  Lambda(1) = " + std::to_string(lambda_b);
            throw std::runtime_error(error_msg);
        }
        // 如果 lambda_a 或 lambda_b 已经非常接近0，那么对应的 a 或 b 就是解。
        if (std::abs(lambda_a) < eps) return 0.0;
        if (std::abs(lambda_b) < eps) return 1.0;


        // 二分查找流程
        while ((b - a) > eps) {
            double m = (a + b) / 2.0;
            double lambda_m = calculateLambda(m);

            if (lambda_m > 0) { // 如果 Lambda(m) > 0, 解在 [m, b] 区间
                a = m;
            } else { // 否则解在 [a, m] 区间
                b = m;
            }
        }
        return (a + b) / 2.0; // 返回区间的中间值作为 kappa*
    }

private:
    double R_;    // 转弯半径: v_hat / c
    double beta_; // 转弯角 (弧度)

    /**
     * @brief 计算路径长度差函数 Lambda(kappa)。
     * Lambda(kappa) = L(kappa) - L_orig
     *
     * @param kappa 参数 kappa 取值范围 [0,1]。
     * @return Lambda(kappa) 的值。
     */
    double calculateLambda(double kappa) const {
        // Xi(kappa, beta) = ((1+kappa) + (1-kappa)*sin(beta/2)) / 2
        double sin_beta_half = std::sin(beta_ / 2.0);
        double Xi = ((1.0 + kappa) + (1.0 - kappa) * sin_beta_half) / 2.0;

        // 数值稳定性：对 Xi 进行限幅，避免 acos 和 sqrt 的参数超出定义域
        // delta 值用于防止 Xi 精确等于 0 或 1 导致的问题
        const double Xi_delta = 1e-12; // 调整后的 delta 以增强鲁棒性
        Xi = std::max(Xi_delta, std::min(Xi, 1.0 - Xi_delta));
        // 当 Xi 趋近于 1 时，acos(Xi) 趋近于 0，Xi*sqrt(1/Xi^2 - 1) 趋近于 0。
        // 当 Xi 趋近于 0 时 (在 kappa=[0,1] 和 beta_=[0,PI] 条件下，Xi 的最小值为 0.5*(1+0)=0.5，所以 Xi 不会接近0)。

        double term_acos_Xi = std::acos(Xi);
        double term_sqrt_Xi_expr; // 即 Xi*sqrt(1/Xi^2 - 1)

        // 对 Xi*sqrt(1/Xi^2 - 1) 项进行稳定计算
        // 这个表达式等价于 sqrt(1 - Xi^2) 当 Xi > 0
        // 或者 Xi * sqrt((1-Xi)(1+Xi))/Xi^2 = sqrt(1-Xi^2) for Xi in (0,1)
        if (Xi >= 1.0 - Xi_delta) { // 如果 Xi 非常接近 1 (或等于1)
            term_sqrt_Xi_expr = 0.0; // Xi * sqrt(1/Xi^2 - 1) 在 Xi->1 时的极限是 0
                                     // 同时 acos(Xi) 此时也为 0
        } else if (Xi <= Xi_delta) { // 如果 Xi 非常接近 0 (理论上不应发生)
             term_sqrt_Xi_expr = 0.0; // 避免潜在的除零或 sqrt 负数
        }
        else {
            // term_sqrt_Xi_expr = Xi * std::sqrt(1.0 / (Xi * Xi) - 1.0);
            // 更稳定的形式: sqrt(1 - Xi*Xi)
            // 因为 1/Xi^2 - 1 = (1-Xi^2)/Xi^2, 所以 Xi * sqrt((1-Xi^2)/Xi^2) = Xi * sqrt(1-Xi^2)/|Xi|
            // 由于 Xi > 0, 所以是 sqrt(1-Xi^2)
            if (1.0 - Xi * Xi < 0) { // 进一步的保护，防止浮点误差导致根号内为负
                term_sqrt_Xi_expr = 0.0;
            } else {
                term_sqrt_Xi_expr = std::sqrt(1.0 - Xi * Xi);
            }
        }

        double cos_beta_half = std::cos(beta_ / 2.0);
        double cot_beta_half;
        // beta_ 已经在 findKappaStar 中检查过，不应为0。
        // sin(beta_/2.0) 也因此不为0。
        if (std::abs(sin_beta_half) < 1e-9) { // 额外的保护层，理论上不应触发
            // 如果 beta_ 非常接近 0 或 2*PI, cot_beta_half 会变得非常大或未定义
            // 此情况应已由 findKappaStar 中的 beta_ 检查处理
            throw std::runtime_error("计算 cot(beta/2) 时出错：beta/2 的正弦值过小。");
        } else {
            cot_beta_half = cos_beta_half / sin_beta_half;
        }

        /* 根据文档中的公式 (42)
           Λ(κ) = 2R * ( (π-β)/2
                         + 2*acos(Ξ(κ,β))
                         - Ξ(κ,β)*sqrt(1/Ξ(κ,β)^2 - 1)
                         - (1-κ)*cos(β/2)
                         - κ*cot(β/2) )
        */
        double lambda_val = (M_PI - beta_) / 2.0
                          + 2.0 * term_acos_Xi
                          - term_sqrt_Xi_expr // 使用了 sqrt(1-Xi^2) 的形式
                          - (1.0 - kappa) * cos_beta_half
                          - kappa * cot_beta_half;

        return 2.0 * R_ * lambda_val;
    }

// 定义 PI (如果cmath中没有定义 M_PI)
#ifndef M_PI
    static constexpr double M_PI = 3.14159265358979323846;
#endif
};

// --- 示例用法 ---
// 编译 (例如, 使用 g++):
// g++ -std=c++17 your_file_name.cpp -o kappa_solver_2d
// ./kappa_solver_2d
int main() {
    KappaTrajectorySolver solver;

    // 测试用例 1: 90度转弯
    Point2D w_prev1 = {0, 0};
    Point2D w_i1    = {5, 0};
    Point2D w_next1 = {5, 5}; // 沿x轴然后沿y轴，beta = PI/2 (90度)

    double v_hat = 10.0; // 预设速度 (m/s)
    double c     = 1.0;  // 转向速率 (rad/s)
    double eps   = 1e-7; // 精度

    std::cout << "测试用例 1: 90度转弯" << std::endl;
    try {
        double kappa_star1 = solver.findKappaStar(w_prev1, w_i1, w_next1, v_hat, c, eps);
        std::cout << "  最优 kappa*: " << kappa_star1 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  错误: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // 测试用例 2: 45度转弯
    // w_i - w_prev = (5,0) -> t_hat_im1_i = (1,0)
    // w_next - w_i = (cos(45), sin(45)) = (1/sqrt(2), 1/sqrt(2))
    // beta = acos( (1,0) dot (1/sqrt(2), 1/sqrt(2)) ) = acos(1/sqrt(2)) = PI/4 (45度)
    Point2D w_prev2 = {0,0};
    Point2D w_i2    = {5,0};
    Point2D w_next2 = {5.0 + std::cos(M_PI/4.0), 0.0 + std::sin(M_PI/4.0)};

    std::cout << "测试用例 2: 45度转弯 (beta = 45 度)" << std::endl;
    try {
        double kappa_star2 = solver.findKappaStar(w_prev2, w_i2, w_next2, v_hat, c, eps);
        std::cout << "  最优 kappa*: " << kappa_star2 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  错误: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // 测试用例 3: 较缓和的转弯 (beta = 135度)
    // t_hat_im1_i = (1,0)
    // t_hat_i_ip1 = (-1/sqrt(2), 1/sqrt(2)) (即从x正向转到左上135度方向)
    // beta = acos( (1,0) dot (-1/sqrt(2), 1/sqrt(2)) ) = acos(-1/sqrt(2)) = 3*PI/4 (135度)
    Point2D w_prev3 = {0,0};
    Point2D w_i3    = {5,0};
    Point2D w_next3 = {5.0 + std::cos(3.0*M_PI/4.0), 0.0 + std::sin(3.0*M_PI/4.0)};
    std::cout << "测试用例 3: 缓和转弯 (beta = 135 度)" << std::endl;
    try {
        double kappa_star3 = solver.findKappaStar(w_prev3, w_i3, w_next3, v_hat, c, eps);
        std::cout << "  最优 kappa*: " << kappa_star3 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  错误: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // 测试用例 4: 几乎直线 (beta 接近 0)
    Point2D w_prev4 = {0,0};
    Point2D w_i4    = {5,0};
    Point2D w_next4 = {10, 0.0001}; // beta 会非常小但非零
    std::cout << "测试用例 4: 几乎直线 (beta 接近 0)" << std::endl;
    try {
        // 预期 findKappaStar 中的 beta_ < angle_tolerance 条件会捕获此情况并返回 0.0
        double kappa_star4 = solver.findKappaStar(w_prev4, w_i4, w_next4, v_hat, c, eps);
        std::cout << "  最优 kappa*: " << kappa_star4 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  '几乎直线'时发生错误: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // 测试用例 5: 几乎U型转弯 (beta 接近 PI)
    Point2D w_prev5 = {0,0};
    Point2D w_i5    = {5,0};
    Point2D w_next5 = {0, 0.0001}; // (w_next - w_i) 约等于 (-5, 0.0001)
                                  // t_hat_im1_i = (1,0)
                                  // t_hat_i_ip1 约等于 (-1, small_val)
                                  // 点积接近 -1, beta 接近 PI
    std::cout << "测试用例 5: 几乎U型转弯 (beta 接近 PI)" << std::endl;
    try {
        // 预期 findKappaStar 中的 std::abs(beta_ - M_PI) < angle_tolerance 条件会捕获
        double kappa_star5 = solver.findKappaStar(w_prev5, w_i5, w_next5, v_hat, c, eps);
        std::cout << "  最优 kappa*: " << kappa_star5 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  '几乎U型转弯'时发生错误: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // 测试用例 6: 重合航点 (应该抛出异常)
    Point2D w_prev6 = {0,0};
    Point2D w_i6    = {0,0};
    Point2D w_next6 = {5,5};
    std::cout << "测试用例 6: 重合航点 (w_prev = w_i)" << std::endl;
    try {
        double kappa_star6 = solver.findKappaStar(w_prev6, w_i6, w_next6, v_hat, c, eps);
        std::cout << "  最优 kappa*: " << kappa_star6 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  捕获到预期错误: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;


    return 0;
}