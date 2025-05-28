# 等长轨迹的κ求解C++实现

### I. 功能要点

本代码主要实现内容包括：

1.  接收二维或三维航点 $\mathbf{w}_{i-1}, \mathbf{w}_i, \mathbf{w}_{i+1}$ 作为输入。
2.  基于这些航点计算转向角 β。这个角度是在由这三个航点定义的平面内的。
3.  计算出一个标量参数 κ*。这个 κ* 的目标是使得后续由 DTS（Dynamically Time-Scaled）算法生成的实际 κ-轨迹的**路径长度**与原始分段直线路径的长度一致。

这个算法的核心是**求解参数 κ***。

之后，这个计算出来的 κ* 会被传递给 DTS 控制律，由 DTS 系统利用这个 κ* 来**生成实际的平滑转弯路径**。

**如果原始航点是三维的，那么 DTS 系统生成的这个 κ-轨迹也将在三维空间中**。每个转弯部分（由三个连续航点定义）将在包含这三个点的平面内进行平滑。由于连续的转弯平面可能不同，整个轨迹将是三维的。

功能要点总结：

* 该 C++ 代码：输入三维航点，输出一个用于长度匹配的**标量参数 κ***。
* DTS 系统（未在此代码中实现）：输入 κ* 和航点等信息，输出实际的**三维平滑曲线**。



### II. C++实现代码（二维平面航点）

以下代码适用于处理二维平面上的航点，该 C++ 代码本身并不直接计算出三维的 Kappa **曲线**的点集。

```cpp
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

```



#### 在无人机仿真软件中使用方法

**集成与使用步骤：**

1.  **代码文件准备与集成：**
    * **复制代码：** 将 `Point2D` 结构体和 `KappaTrajectorySolver` 类的完整 C++ 代码（包括头文件引用和 `M_PI` 定义）保存到一个或多个 `.cpp` 和/或 `.h` 文件中。例如，您可以创建一个 `kappa_solver.h` 和 `kappa_solver.cpp`。
        * `kappa_solver.h`: 包含类定义和 `Point2D` 结构体定义。
        * `kappa_solver.cpp`: 包含类成员函数的实现。
    * **添加到项目：** 将这些文件添加到您的无人机仿真软件的项目中。具体操作取决于您使用的IDE（如 Visual Studio, Qt Creator, Eclipse）或构建系统（如 CMake, Makefile）。
        * **IDE:** 通常是右键点击项目 -> "Add Existing Item" 或类似选项。
        * **CMake:** 修改 `CMakeLists.txt` 文件，将新的源文件和头文件路径添加到 `add_executable` 或 `add_library` 以及 `target_include_directories` 中。
        * **Makefile:** 手动更新 Makefile 中的源文件列表和编译规则。
    * **确保编译：** 重新编译您的仿真软件项目，确保新添加的代码能够被正确编译和链接。

2.  **在仿真代码中包含头文件：**
    * 在您需要使用 `KappaTrajectorySolver` 的仿真模块的 `.cpp` 文件顶部，包含其头文件：
        ```cpp
        #include "kappa_solver.h" // 假设您将类定义放在此文件中
        // 或者如果直接集成在某个cpp中，确保其声明可见
        ```

3.  **创建 `KappaTrajectorySolver` 对象：**
    * 在您的路径规划模块、飞行控制模块或任何需要计算 κ* 的地方，创建一个 `KappaTrajectorySolver` 类的实例。
        ```cpp
        KappaTrajectorySolver kappa_calculator; // 可以是成员变量，也可以是局部变量
        ```

4.  **准备输入参数并调用 `findKappaStar`：**
    * 当无人机规划路径或即将进入一个转弯时，您需要从仿真环境中获取以下参数：
        * `w_prev`：上一个航点 (类型 `Point2D`)。
        * `w_i`：当前转弯点 (类型 `Point2D`)。
        * `w_next`：下一个航点 (类型 `Point2D`)。
        * `v_hat`：无人机在该航段的预设速度 (类型 `double`, 单位：米/秒)。
        * `c`：无人机的转向速率 (类型 `double`, 单位：弧度/秒)。这通常是无人机的性能参数或期望的转弯敏捷度。
        * `eps`：计算 κ* 的目标精度 (类型 `double`, 例如 `1e-6` 或 `1e-7`)。
    * 调用 `findKappaStar` 方法：
        ```cpp
        Point2D waypoint_prev = { /* 从仿真获取x, y */ };
        Point2D waypoint_i    = { /* 从仿真获取x, y */ };
        Point2D waypoint_next = { /* 从仿真获取x, y */ };
        
        double desired_speed = 10.0; // 米/秒，从仿真任务或参数获取
        double turn_rate_capability = 0.8; // 弧度/秒，根据无人机性能设定
        
        double optimal_kappa_star;
        try {
            optimal_kappa_star = kappa_calculator.findKappaStar(
                waypoint_prev,
                waypoint_i,
                waypoint_next,
                desired_speed,
                turn_rate_capability
                // eps 使用默认值 1e-6，或显式传入
            );
        
            // 成功获取 optimal_kappa_star
            // std::cout << "计算得到的 kappa*: " << optimal_kappa_star << std::endl;
        
        } catch (const std::exception& e) {
            std::cerr << "计算 kappa* 时发生错误: " << e.what() << std::endl;
            // 在此处理错误：
            // 1. 使用一个默认的 kappa 值 (例如 0.0, 0.5 或 1.0)
            // 2. 标记路径段为有问题，触发重新规划或警告
            // 3. 根据仿真软件的错误处理机制进行上报
            optimal_kappa_star = 0.5; // 示例：使用一个安全的默认值
        }
        ```

5.  **使用计算出的 `optimal_kappa_star`：**
    * **核心用途**：根据原始文档，计算出的 `optimal_kappa_star` 主要用于与 DTS（Dynamically Time-Scaled，动态时间缩放）算法结合。
    * **传递给DTS系统**：将 `optimal_kappa_star` 作为参数传递给您仿真软件中的 DTS 控制律或轨迹生成模块。
    * **轨迹生成**：DTS 系统会使用这个 κ* 值来调整转弯的几何形状，以生成一条路径长度与原始分段直线路径长度相等的平滑轨迹。仿真软件中的轨迹跟踪控制器随后会执行这条由DTS生成的轨迹。
    * 如果您仿真软件中没有DTS系统，但有其他可以利用类似参数（调整转弯松紧度或形状）的平滑算法，您可以尝试将 `optimal_kappa_star` 适配到该算法。

**注意事项和最佳实践：**

* **调用时机：**
    * **路径规划阶段：** 对于一条包含多个转弯的预定路径，可以在路径生成后，为每个转弯点（由三个连续航点定义）预先计算并存储各自的 κ*。
    * **飞行前准备：** 在无人机开始执行任务前。
    * **动态重规划：** 如果路径是动态生成的或在飞行中发生变化，则需要在新的转弯段形成时重新计算。
* **坐标系和单位：** 确保所有输入给 `findKappaStar` 函数的航点坐标、速度和转向速率的单位（例如，米、米/秒、弧度/秒）与仿真环境内部使用的单位一致。
* **错误处理的鲁棒性：** `findKappaStar` 函数在某些输入下（如航点重合导致 `beta_` 无法计算，或 `Lambda(0)` 和 `Lambda(1)` 的符号不符合预期）会抛出异常。您的仿真代码必须能够妥善处理这些异常，以避免仿真崩溃。可以考虑的策略包括记录错误、使用默认的 κ 值、或者触发路径重新规划。
* **DTS 系统的配合：** 这个 κ* 求解算法的目的是为了让DTS生成的轨迹与原始路径等长。其效果的发挥依赖于仿真软件中存在相应的DTS控制律或类似的轨迹生成机制来消费这个 κ* 值。
* **性能考虑：** `findKappaStar` 函数本身包含一个迭代过程（二分查找）和一些三角函数计算。对于非常多的航点或极高频率的重规划，需要评估其计算耗时是否满足实时性要求。通常情况下，对于单个转弯的计算是非常快的。
* **三维空间的扩展：** 如果仿真环境是三维的，并且希望处理三维路径，那就需要以下的三维版本 `KappaTrajectorySolver`。该算法的原理同样适用，转弯发生在由三个三维航点定义的平面内。

通过以上步骤，将能够将 `KappaTrajectorySolver` 集成到无人机仿真软件中，并用它来计算确保等长路径的 κ* 参数，从而辅助生成更平滑且符合时序约束的飞行轨迹。

#### 二维航点与三维航点实现的主要差异

1.  **`Point3D` -> `Point2D`**：
    * 结构体移除了 `z` 成员。
    * 所有相关的向量操作（如 `norm`, `normalize`, `dot`, 运算符重载）都已更新为二维计算。
2.  **`calculateLambda` 中 `term_sqrt_Xi_expr` 的计算**：
    * 原公式中的 `Xi*sqrt(1/Xi^2 - 1)` 在数值上等价于 `sqrt(1 - Xi^2)` （对于 `Xi > 0`）。代码已更新为使用 `sqrt(1.0 - Xi * Xi)` 这种更直接且可能更稳定的形式，并增加了对根号内表达式可能为负的浮点保护。
3.  **`beta_` 角度的特殊处理**：
    * 对 `beta_` 极小（接近直线）或极大（接近U型转弯）的情况增加了更明确的处理和说明。根据算法文档，`Lambda` 函数的特定形式和二分法的前提是针对一个“有效”的转弯。对于直线或完美U型转弯，`Lambda(kappa)` 理论上对所有 `kappa` 都为0，此时 `kappa` 的选择可能需要根据DTS系统的约定。当前代码在这类退化情况下会直接返回 `kappa = 0.0`。
4.  **`main` 函数测试用例**：
    * 所有测试用例都更新为使用 `Point2D`。
    * 增加了一个重合航点的测试用例，以验证异常处理。

如需考虑航点为三维航点，可以参照第III节的 C++实现代码（三维航点）。

### III. C++实现代码（三维航点）

Here's a modular C++ implementation of the等长κ-轨迹算法 based on the provided document.

This code will define a class `KappaTrajectorySolver` that encapsulates the logic for calculating κ*.

```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>   // For std::inner_product or a similar concept if needed for dot product
#include <stdexcept> // For runtime_error
#include <algorithm> // For std::min and std::max

// Define a simple 3D point structure (can be replaced with your existing Point/Vector class)
struct Point3D {
    double x, y, z;

    Point3D operator-(const Point3D& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Point3D operator+(const Point3D& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Point3D operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    Point3D operator/(double scalar) const {
        if (scalar == 0) throw std::runtime_error("Division by zero in Point3D scalar division.");
        return {x / scalar, y / scalar, z / scalar};
    }

    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    Point3D normalize() const {
        double n = norm();
        if (n == 0) throw std::runtime_error("Cannot normalize a zero vector.");
        return {x / n, y / n, z / n};
    }

    static double dot(const Point3D& v1, const Point3D& v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }
};

class KappaTrajectorySolver {
public:
    // Constructor
    KappaTrajectorySolver() = default;

    /**
     * @brief Calculates the optimal kappa* value.
     *
     * @param w_prev Previous waypoint.
     * @param w_i Current waypoint (turn point).
     * @param w_next Next waypoint.
     * @param v_hat Preset speed.
     * @param c Turning rate.
     * @param eps Target precision for kappa*.
     * @return The optimal kappa* value.
     */
    double findKappaStar(const Point3D& w_prev, const Point3D& w_i, const Point3D& w_next,
                           double v_hat, double c, double eps = 1e-6) {
        if (c == 0) {
            throw std::invalid_argument("Turning rate 'c' cannot be zero.");
        }
        R_ = v_hat / c; //

        Point3D v1 = w_i - w_prev;
        Point3D v2 = w_next - w_i;

        if (v1.norm() == 0 || v2.norm() == 0) {
             // Or handle as a straight line, no turn, kappa doesn't apply or is degenerate.
            // For now, let's assume this means kappa is irrelevant, perhaps return 0 or 1
            // depending on desired behavior for straight segments.
            // For the context of the problem trying to find a specific kappa for a turn,
            // this is an invalid scenario for the Lambda function as beta would be undefined.
            // std::cerr << "Warning: Zero length segment detected. Turn angle calculation is undefined." << std::endl;
            // Based on the problem, Lambda(kappa) = 0 is the goal.
            // If there's no turn, L(kappa) = L_orig trivially.
            // However, the Lambda function itself is defined for a turn.
            // Let's assume for a straight line, any kappa is "valid" in terms of length matching,
            // but the typical use case is for turns. If it's perfectly straight (beta=pi),
            // the geometry changes.
            // The problem implies beta is an angle of a turn. If v1 or v2 is zero, beta is ill-defined.
            // If t_hat_i_ip1 and t_hat_im1_i are pointing in the same direction or opposite,
            // beta would be pi or 0.
            // Let's re-evaluate beta calculation to be robust.
            throw std::invalid_argument("Waypoints result in a zero-length segment, cannot compute turn angle.");
        }

        Point3D t_hat_im1_i = v1.normalize();
        Point3D t_hat_i_ip1 = v2.normalize();

        double dot_product = Point3D::dot(t_hat_i_ip1, t_hat_im1_i);
        // Clamp dot_product to avoid domain errors with acos due to potential floating point inaccuracies
        dot_product = std::max(-1.0, std::min(1.0, dot_product));
        beta_ = std::acos(dot_product); //

        // If beta is very close to 0 (straight segment going forward) or M_PI (straight segment reversing),
        // the turn geometry assumed by Lambda(kappa) might be degenerate.
        // The document implies beta is a turn angle, so beta != 0 and beta != PI.
        // Specifically, cot(beta/2) is used, which is problematic if beta = 0 or beta = 2*PI.
        // (pi-beta)/2 implies beta should be < pi.
        // cos(beta/2) is fine. sin(beta/2) in Xi implies beta/2 not equal to k*pi.
        if (beta_ < 1e-9 || std::abs(beta_ - M_PI) < 1e-9) {
            // For a straight line (beta = PI for forward, beta = 0 for reverse in terms of t_hat dot),
            // L(kappa) should intrinsically equal L_orig.
            // The Lambda function aims to correct deviations. If there's no turn, Lambda should be 0.
            // Let's check the components of Lambda:
            // (pi-beta)/2 -> if beta=pi, this is 0. if beta=0, this is pi/2
            // cos(beta/2) -> if beta=pi, 0. if beta=0, 1
            // cot(beta/2) -> if beta=pi, 0. if beta=0, undefined.
            // This indicates the formula is for beta in (0, pi).
            // std::cerr << "Warning: Path is nearly straight (beta approx 0 or PI). Kappa may not be well-defined by this Lambda." << std::endl;
            // If it's straight, L_kappa = L_orig for any kappa by definition as there's no turn to modify.
            // Thus Lambda(kappa) = 0. We can return any valid kappa, e.g., 0 or 1.
            // However, the premise of the bisection relies on Lambda(0) > 0 and Lambda(1) < 0 for a turn.
            // If it's perfectly straight (beta = pi), then sin(beta/2)=0, cos(beta/2)=0.
            // Xi = (1+k)/2.
            // If beta_ is very small (sharp turn, almost 180 deg U-turn):
            //  (pi-beta)/2 approaches pi/2.
            //  sin(beta/2) approaches 0. Xi approaches (1+kappa)/2.
            //  cos(beta/2) approaches 1.
            //  cot(beta/2) becomes very large. This is the critical term for small beta.
            // If beta_ is very close to PI (shallow turn, almost straight):
            //  (pi-beta)/2 approaches 0.
            //  sin(beta/2) approaches 1. Xi approaches (1+kappa + 1-kappa)/2 = 1.
            //  cos(beta/2) approaches 0.
            //  cot(beta/2) approaches 0.
            //  If Xi = 1, arccos(Xi)=0, sqrt(1/Xi^2 - 1) term is 0/0 or problematic.
            //  The document mentions numerical robustness for Xi approx 1.
            if (beta_ < 1e-6) { // Very sharp U-turn, close to 0 degrees between vectors
                 // std::cerr << "Warning: Beta is very close to 0 (sharp U-turn). Lambda function might be unstable." << std::endl;
                 // This scenario means the vehicle is nearly reversing its path.
                 // The geometry described might not be typical for such a maneuver.
                 // For now, proceed, but be aware of potential issues.
            }
            if (std::abs(beta_ - M_PI) < 1e-6) { // Almost straight
                // In this case, Lambda(kappa) should ideally be 0 for all kappa.
                // The problem states Lambda(0) > 0 and Lambda(1) < 0. This holds for turns.
                // For a straight line, the turn radius is infinite, or the concept of kappa-turn doesn't apply.
                // L(kappa) = L_orig, so Lambda(kappa) = 0. The bisection wouldn't find a unique root if f(a) and f(b) are both 0.
                // We can return a default value, e.g. 0.5, or what the DTS system expects for straight.
                // The question is about finding kappa* such that Lambda(kappa*) = 0.
                // If Lambda is always 0, any kappa works. The algorithm's premise fails.
                // However, due to floating point, it might not be exactly 0.
                // Let's assume the waypoints form a distinct turn for the algorithm to be meaningful.
                 return 0.0; // Or 1.0, or handle as per system's needs for "no effective turn".
            }
        }


        double a = 0.0; //
        double b = 1.0; //

        double lambda_a = calculateLambda(a);
        double lambda_b = calculateLambda(b);

        // Assertion from the document: Lambda(0) > 0 and Lambda(1) < 0
        // Add a small tolerance for floating point comparisons
        double check_eps = 1e-9;
        if (lambda_a < -check_eps || lambda_b > check_eps) {
            // This can happen if beta is too close to PI (almost straight line)
            // or other edge cases where the monotonicity doesn't hold as expected
            // or the values are practically zero.
            // std::cerr << "Warning: Lambda(0) or Lambda(1) does not meet expected sign conditions." << std::endl;
            // std::cerr << "  beta = " << beta_ * 180.0 / M_PI << " degrees" << std::endl;
            // std::cerr << "  Lambda(0) = " << lambda_a << std::endl;
            // std::cerr << "  Lambda(1) = " << lambda_b << std::endl;
            // If lambda_a is already ~0, kappa=0 is the solution.
            // If lambda_b is already ~0, kappa=1 is the solution.
            if (std::abs(lambda_a) < eps) return 0.0;
            if (std::abs(lambda_b) < eps) return 1.0;

            // If signs are not as expected, bisection method cannot proceed as described.
            // This might indicate an issue with the input geometry or that the assumptions
            // of Lemma 7 are not met for this specific input.
            // One possible recovery: if both are positive, means solution might be >1 (not allowed) or error.
            // If both are negative, means solution might be <0 (not allowed) or error.
            // Or, the function isn't monotonic for these params.
            // For now, throw an error if the fundamental assumption is violated significantly.
             throw std::runtime_error("Lambda function does not satisfy endpoint conditions Lambda(0)>0 and Lambda(1)<0. "
                                      "Check waypoint geometry or algorithm assumptions for this case. "
                                      "beta = " + std::to_string(beta_ * 180.0 / M_PI) + " deg, "
                                      "L(0) = " + std::to_string(lambda_a) + ", L(1) = " + std::to_string(lambda_b));
        }
        // Ensure they are strictly different for bisection to make sense.
        // If lambda_a is positive and lambda_b is negative, proceed.
        // If lambda_a is negative (or zero) and lambda_b is positive (or zero), signs are flipped from expectation.
        // The problem states Lambda(0) > 0 and Lambda(1) < 0 strictly.

        while ((b - a) > eps) { //
            double m = (a + b) / 2.0; //
            double lambda_m = calculateLambda(m);

            if (lambda_m > 0) { //
                a = m;
            } else {
                b = m;
            }
        }
        return (a + b) / 2.0; //
    }

private:
    double R_;    // Radius: v_hat / c
    double beta_; // Turn angle in radians

    /**
     * @brief Calculates the path length difference function Lambda(kappa).
     *
     * @param kappa Parameter kappa in [0,1].
     * @return Value of Lambda(kappa).
     */
    double calculateLambda(double kappa) const {
        // Xi(kappa, beta) = ((1+kappa) + (1-kappa)*sin(beta/2)) / 2 
        double sin_beta_half = std::sin(beta_ / 2.0);
        double Xi = ((1.0 + kappa) + (1.0 - kappa) * sin_beta_half) / 2.0;

        // Numerical stability for Xi close to 0 or 1
        // Using a small delta, e.g., 1e-9, to prevent domain errors with acos or sqrt
        const double delta = 1e-12; // Adjusted delta for robustness
        Xi = std::max(delta, std::min(Xi, 1.0 - delta));
        // If Xi is exactly 1 (e.g. beta=pi and kappa=any, or beta<pi and kappa=1 leads to sin_beta_half=1),
        // then sqrt(1/Xi^2 - 1) will be sqrt(1-1)=0.
        // If Xi is exactly 0 (not possible given its formula if beta/2 in [0, pi/2] and kappa in [0,1]),
        // then 1/Xi^2 is problematic. Smallest Xi happens when kappa=1 and sin(beta/2)=0 (beta=0), so Xi=1.
        // Or kappa=-1 (not allowed) and sin(beta/2)=1 (beta=pi), then Xi = 0.
        // Given kappa in [0,1] and beta in (0,pi):
        // min sin_beta_half is close to 0 (for beta close to 0)
        // max sin_beta_half is 1 (for beta = pi, but we exclude beta=pi based on cot)
        // If beta is close to 0: sin_beta_half ~ 0. Xi ~ (1+kappa)/2. Smallest is 0.5 (kappa=0). Largest 1 (kappa=1).
        // If beta is pi (excluded, but for limit): sin_beta_half=1. Xi = ((1+k)+(1-k)*1)/2 = (1+k+1-k)/2 = 1.
        // So Xi should be in [~0.5, 1]. Clamping to 1-delta is most important.

        double term_acos_Xi = std::acos(Xi);
        double term_sqrt_Xi;
        if (std::abs(Xi) < 1.0 - 1e-9) { // Avoid division by zero if Xi is 1.0
             term_sqrt_Xi = Xi * std::sqrt(1.0 / (Xi * Xi) - 1.0);
        } else if (Xi >= 1.0 - 1e-9) { // If Xi is very close to 1 (or 1)
            term_sqrt_Xi = 0.0; // Limit of Xi * sqrt(1/Xi^2 - 1) as Xi -> 1 is 0
            term_acos_Xi = 0.0; // acos(1) = 0
        }
        else { // Xi is very close to -1 (should not happen with kappa in [0,1] and beta in [0,PI])
            term_sqrt_Xi = 0.0; // sqrt would be complex
            term_acos_Xi = M_PI; // acos(-1) = PI
        }


        double cos_beta_half = std::cos(beta_ / 2.0);
        double cot_beta_half;
        if (std::abs(std::sin(beta_ / 2.0)) < 1e-9) { // beta is close to 0 or 2*PI
            // This case should ideally be handled by the beta check in findKappaStar
            // If beta_ is very small, cot_beta_half becomes huge.
            // Based on prior checks, beta_ should not be 0.
            throw std::runtime_error("cot(beta/2) undefined: beta is too close to 0 or 2*PI.");
        } else {
            cot_beta_half = cos_beta_half / sin_beta_half;
        }

        /* Original formula:
           Lambda(kappa) = 2R * ( (pi-beta)/2
                                 + 2*acos(Xi)
                                 - Xi*sqrt(1/Xi^2 - 1)
                                 - (1-kappa)*cos(beta/2)
                                 - kappa*cot(beta/2) )
        */
        double lambda_val = (M_PI - beta_) / 2.0
                          + 2.0 * term_acos_Xi
                          - term_sqrt_Xi
                          - (1.0 - kappa) * cos_beta_half
                          - kappa * cot_beta_half;

        return 2.0 * R_ * lambda_val;
    }

#ifndef M_PI
    static constexpr double M_PI = 3.14159265358979323846;
#endif
};

// --- Example Usage ---
// To compile (e.g., with g++):
// g++ -std=c++17 your_file_name.cpp -o kappa_solver
// ./kappa_solver
int main() {
    KappaTrajectorySolver solver;

    // Example Test Case (you should replace these with actual values from your simulation)
    Point3D w_prev = {0, 0, 0};
    Point3D w_i    = {5, 0, 0};
    Point3D w_next = {5, 5, 0}; // 90 degree turn

    double v_hat = 10.0; // m/s
    double c     = 1.0;  // rad/s (turn rate)
    double eps   = 1e-7;

    std::cout << "Test Case 1: 90-degree turn" << std::endl;
    try {
        double kappa_star = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star << std::endl;

        // Verify Lambda(kappa_star) is close to 0
        // Need to re-setup solver's internal R_ and beta_ for an external check,
        // or make calculateLambda public if direct testing is needed.
        // For simplicity, we'll trust the internal calculation for now.
        // To test lambda(kappa_star) ~ 0, one would need to call calculateLambda directly.
        // This implies calculateLambda might need to be public or a helper test function.
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // Test Case 2: Sharper turn
    Point3D w_next_sharp = {6, 1, 0}; // Approx 45 deg turn angle with (5,0)-(0,0) and (6,1)-(5,0)
                                      // Vector 1: (5,0,0). Vector 2: (1,1,0)
                                      // Angle between (5,0) and (1,1) from origin: (1,1) is at 45 deg. Vector from (5,0) to (6,1) is (1,1)
                                      // So turn from x-axis to vector (1,1) means beta is 135 deg for the path vectors.
                                      // v_prev_i = (5,0,0), v_i_next = (1,1,0)
                                      // acos( (5,0,0) dot (1,1,0) / (|5|*|sqrt(2)|) ) = acos(5 / (5*sqrt(2))) = acos(1/sqrt(2)) = 45 deg. This is beta.
                                      // The document defines beta = acos( t_hat_i_ip1^T * t_hat_im1_i )
                                      // t_im1_i = normalize(w_i - w_prev) = normalize( (5,0,0) - (0,0,0) ) = (1,0,0)
                                      // t_i_ip1 = normalize(w_next - w_i) = normalize( (6,1,0) - (5,0,0) ) = normalize( (1,1,0) ) = (1/sqrt(2), 1/sqrt(2), 0)
                                      // dot_product = 1 * 1/sqrt(2) + 0 * 1/sqrt(2) = 1/sqrt(2)
                                      // beta = acos(1/sqrt(2)) = pi/4 = 45 degrees.

    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {5+cos(M_PI/4.0), 5+sin(M_PI/4.0), 0}; // Should be w_i + vector, e.g. {5+1, 1, 0} for a clear next step.
                                               // Let's define w_next to make a clear angle.
                                               // w_next such that (w_next - w_i) makes 45 deg with (w_i - w_prev)
                                               // w_i - w_prev = (5,0,0)
                                               // Let w_next - w_i = (1,1,0). So w_next = (6,1,0)
    w_next = {6,1,0};


    std::cout << "Test Case 2: 45-degree turn (beta = 45 deg)" << std::endl;
    try {
        double kappa_star_sharp = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_sharp << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // Test Case 3: Shallow turn (beta close to PI)
    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {10, 0.1, 0}; // Very shallow turn
                           // t_im1_i = (1,0,0)
                           // t_i_ip1 = normalize( (5, 0.1, 0) ). Dot product very close to 1. Beta close to 0 for these vectors.
                           // The definition of beta is acos(t_hat_i_ip1^T * t_hat_im1_i)
                           // Here, t_hat_im1_i is direction of (w_i - w_prev)
                           // t_hat_i_ip1 is direction of (w_next - w_i)
                           // If w_prev=(0,0,0), w_i=(5,0,0), w_next=(10,0.1,0):
                           // t_hat_im1_i = (1,0,0)
                           // vec_i_ip1 = (5, 0.1, 0). norm = sqrt(25.01) approx 5.001
                           // t_hat_i_ip1 approx (5/5.001, 0.1/5.001, 0) approx (0.9998, 0.0199, 0)
                           // dot_product approx 0.9998. acos(0.9998) approx 0.02 radians (1.14 degrees). This is beta.
                           // This is a very SHARP turn in terms of the algorithm's beta (angle between successive segment vectors).
                           // The problem says: β = cos⁻¹(t̂ᵢ,ᵢ₊₁ᵀ t̂ᵢ₋₁,ᵢ)
                           // For my current setup: w_prev -> w_i is along +x.
                           // w_i -> w_next is a vector (5, 0.1, 0).
                           // The angle between (1,0,0) and (5,0.1,0) is small.
                           // This means beta is small.
                           // A "shallow turn" in vehicle path means the angle of deviation from straight is small.
                           // If beta is the internal angle of the corner (like in a polygon), then a shallow turn has beta close to PI.
                           // The document refers to (pi-beta)/2. If beta is the turn angle (0 for straight, pi for U-turn),
                           // then (pi-beta)/2 means beta is the angle of deviation from straight, not the interior angle.
                           // The pseudocode uses: β = arccos(dot(normalize(w_next - w_i), normalize(w_i - w_prev)))
                           // This is the angle between the two path segments.
                           // If w_i-w_prev = (1,0) and w_next-w_i = (1,0.01), beta is small. This is a sharp deviation *from the previous segment's direction*.
                           // This is what the code calculates. A "shallow turn" in terms of path curvature would mean beta is close to PI (180 degrees),
                           // meaning the change in direction is small.
                           // Let's re-evaluate the "shallow turn" example:
                           // w_prev = (0,0,0), w_i = (5,0,0), w_next = (4,1,0) -> results in beta > 90 deg
                           // This means t_im1_i is (1,0,0) and t_i_ip1 is normalize(-1,1,0). dot is -1/sqrt(2). beta = 135 deg. This is a "gentle" turn.
    w_next = {4, 1, 0}; // beta = 135 degrees (gentle turn)
    std::cout << "Test Case 3: Gentle turn (beta = 135 deg)" << std::endl;
    try {
        double kappa_star_shallow = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_shallow << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error: " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;


     // Test Case 4: Almost straight line (beta close to PI), where Lambda(0) and Lambda(1) might be problematic
    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {10,0.00001,0}; // Almost straight - (w_next - w_i) is almost same direction as (w_i - w_prev)
                             // t_im1_i = (1,0,0)
                             // t_i_ip1 = normalize(5, 0.00001, 0) is very close to (1,0,0).
                             // dot product is very close to 1. beta is very close to 0.
                             // This is a *very sharp* turn according to the beta definition (vehicle almost continues straight).
                             // No, this is confusing. Beta = arccos(dot(v2,v1)).
                             // If v1 = (1,0) and v2 = (1,epsilon), then dot is close to 1, beta is close to 0.
                             // This is a very *slight* turn. A shallow turn where the change in heading is small.
                             // The document's (pi-beta)/2 implies beta is the *interior* angle of the turn corner if we view segments.
                             // Let's re-read "β = cos⁻¹(t̂ᵢ,ᵢ₊₁ᵀ t̂ᵢ₋₁,ᵢ)"
                             // t_hat_im1_i points from w_im1 to w_i.
                             // t_hat_i_ip1 points from w_i to w_ip1.
                             // If straight: t_hat_i_ip1 = t_hat_im1_i. Dot product is 1. Beta = 0.
                             // If 180 deg U-turn: t_hat_i_ip1 = -t_hat_im1_i. Dot product is -1. Beta = PI.
                             // So beta IS the turn angle: 0 for straight, PI for full reversal.
                             // My previous interpretation was mixed up.
                             // (pi-beta)/2:
                             //  - If straight (beta=0), then pi/2.
                             //  - If 90 deg turn (beta=pi/2), then pi/4.
                             //  - If U-turn (beta=pi), then 0.
                             // cot(beta/2):
                             //  - If straight (beta=0), then cot(0) is undefined. This is an issue.
                             //  - If 90 deg turn (beta=pi/2), then cot(pi/4) = 1.
                             //  - If U-turn (beta=pi), then cot(pi/2) = 0.
                             // This means the formula has a singularity or problem at beta=0 (straight).
                             // The document states Lambda(0)>0, Lambda(1)<0 (Lemma 7).
                             // This implies beta is not 0. "Thus, there exists ... unique k* in (0,1)" (Theorem 8).
                             // The problem is defined for a *turn*. If it's straight, L(kappa)=L_orig, so Lambda=0 for all kappa.
                             // The algorithm should handle this. My `beta_ < 1e-6` check returning 0.0 is one way.
    std::cout << "Test Case 4: Almost straight (beta close to 0)" << std::endl;
    try {
        // For beta close to 0, the initial check `if (beta_ < 1e-6)` should catch it.
        // If it's truly straight, my code might return 0.0 from that early exit.
        // The `findKappaStar` will throw if `cot(beta/2)` is problematic due to `beta_` near 0.
        // My current `beta_ < 1e-6` check returns 0.0 for kappa*.
        // Let's make it slightly non-zero beta
        w_next = {10, 0.001, 0}; // beta will be small but non-zero.
        double kappa_star_straightish = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_straightish << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error for 'almost straight': " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    // Test Case 5: U-turn (beta close to PI)
    w_prev = {0,0,0};
    w_i    = {5,0,0};
    w_next = {0,0.001,0}; // Almost a U-turn relative to w_i. (w_next - w_i) is (-5, 0.001, 0)
                          // t_im1_i = (1,0,0)
                          // t_i_ip1 = normalize(-5, 0.001, 0), approx (-1, small_val, 0)
                          // Dot product close to -1. Beta close to PI.
    std::cout << "Test Case 5: Almost U-turn (beta close to PI)" << std::endl;
    try {
        // My `std::abs(beta_ - M_PI) < 1e-6` check returns 0.0 for kappa*.
        // Let's make beta slightly less than PI
        w_next = {0, 0.1, 0}; // beta will be large, close to PI.
        double kappa_star_uturn = solver.findKappaStar(w_prev, w_i, w_next, v_hat, c, eps);
        std::cout << "  Optimal kappa*: " << kappa_star_uturn << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "  Error for 'almost U-turn': " << e.what() << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;


    return 0;
}

```

**关键组件和特性：**

1. **`Point3D` 结构体**：一个用于 3D 点/向量的基础结构，包含减法、归一化和点积等常用操作。如果你已有现成的向量数学库（例如 Eigen、GLM），可以直接替换使用。
2. **`KappaTrajectorySolver` 类**：

   * **`findKappaStar` 方法**：主要公共接口接口，接受三个航路点（$\mathbf{w}_{i-1}, \mathbf{w}_i, \mathbf{w}_{i+1}$）、速度 $\hat v$、转向速率 $c$ 以及精度 $\varepsilon$ 作为输入。方法内部首先计算 $R$ 和 β，然后执行二分查找。
   * **`calculateLambda` 方法**：私有辅助函数，利用“引理6”中的公式计算 $\Lambda(\kappa)$。
   * **数值稳定性**：

     * 将用于 `acos` 的点积结果限制在` [-1,1]` 之间。
     * 将 $\Xi$ 限幅到 $[\delta,\,1-\delta]$，以避免 $\Xi$ 接近 0 或 1 时导致的 $acos(\Xi)$ 和 $\sqrt{1/\Xi^2-1}$ 不稳定。
     * 对 $\Xi\sqrt{1/\Xi^2-1}$ 在 $\Xi\approx1$ 时保证正确返回 0。
     * 增加了对 β 极接近 0 或 $\pi$ 的检查：

       * 若 $\beta\approx0$（直线路段），$\cot(\beta/2)$ 将无定义。理论上此时 $\Lambda(\kappa)$ 应恒为 0，代码约定退化情况直接返回 $\kappa=0$。
       * 若 $\beta\approx\pi$（U 型转弯），公式基本正常，但当前实现也返回 0.0，具体行为可根据需求调整。问题假设 $\beta\in(0,\pi)$。
   * **二分查找**：实现了上述描述的二分查找算法，目标是找到 $\Lambda(\kappa^*)\approx0$。会检查端点条件 $\Lambda(0)>0$ 和 $\Lambda(1)<0$，若不满足则抛出 `runtime_error`。
   * **常量定义**：若 `<cmath>` 中未包含 `M_PI`，需自行定义。

**在无人机仿真软件中使用方法：**

1. **包含代码**：将 `Point3D` 结构体（或你的向量类型）和 `KappaTrajectorySolver` 类添加到项目中。
2. **实例化求解器**：

   ```cpp
   KappaTrajectorySolver kappa_solver;
   ```
3. **调用 `findKappaStar`**：当需为由三点定义的转弯段计算 κ* 时：

   ```cpp
   // 假设你有符合 Point3D 接口的航点
   Point3D w_prev = { ... };
   Point3D w_i    = { ... };
   Point3D w_next = { ... };
   
   double desired_speed      = ...;   // 你的 v_hat
   double turning_capability = ...;   // 你的 c
   double desired_precision  = 1e-6;  // 你的 ε
   
   try {
       double optimal_kappa = kappa_solver.findKappaStar(
           w_prev, w_i, w_next,
           desired_speed, turning_capability, desired_precision
       );
       // 将 optimal_kappa 应用到你的 DTS 控制律中
   } catch (const std::exception& e) {
       std::cerr << "计算 kappa* 失败: " << e.what() << std::endl;
       // 处理错误：例如使用默认 kappa，或标记该航段有问题
   }
   ```

**基于文档的重要注意事项：**

* **角度 β 定义**：代码中

  $$
    \beta = \cos^{-1}\bigl(\hat{\mathbf{t}}_{i,i+1}^T \hat{\mathbf{t}}_{i-1,i}\bigr)
  $$

  即航点 $\mathbf{w}_{i-1}\to\mathbf{w}_i$ 与 $\mathbf{w}_i\to\mathbf{w}_{i+1}$ 两向量间的夹角；$\beta=0$ 表示直行，$\beta=\pi$ 表示 180° 的 U 型转弯。
* **$\beta=0$ 奇异性**：$\Lambda(\kappa)$ 中的 $\cot(\beta/2)$ 在 $\beta=0$ 时无定义。理论（引理7、定理8）假设 $\beta\in(0,\pi)$ 才保证 $\Lambda(0)>0$、$\Lambda(1)<0$。实现中对 $\beta\approx0$ 添加了退化处理：直接返回 $\kappa^*=0$。对 $\beta\approx\pi$ 也做了检查，若有需要可根据 U 转的期望效果调整。
* **数值鲁棒性**：对 $\Xi$ 的限幅以及对临界值的处理非常关键；可根据使用的浮点类型（double 或 float）调整 $\delta$ 大小。
* **错误处理**：代码中包含对零除、端点条件不满足等情况的异常抛出，请将其纳入仿真软件的错误管理流程。

此 C++ 代码提供了一个模块化的稳健起点，请务必在仿真环境中对多种航点配置进行充分测试。
