

# Kappa 轨迹生成器 (Kappa Trajectory Generator) Smooth Turn Path Generation for UAVs


这是一个基于 C++ 实现的无人机（UAV）平滑转弯轨迹生成模块。它根据 E. Anderson, R. Beard, T. McLain 在论文 "Real-time dynamic trajectory smoothing for unmanned aerial vehicles" 中提出的 $\\kappa$-轨迹几何定义，在给定的三个航路点之间生成平滑的转弯路径。该模块特别考虑了参数 $\\kappa=0$（路径通过中间航点）和 $\\kappa=1$（最短时间/单圆弧转弯）的特殊情况。

## 核心功能 ✨

  * **平滑转弯生成**：输入前一个航点 ($w\_{i-1}$)、当前转弯点 ($w\_i$) 和下一个航点 ($w\_{i+1}$)，以及飞行速度、最大转向速率和 $\\kappa$ 参数，生成连接入航线段和出航线段的平滑曲线。
  * **$\\kappa$ 参数控制**：通过调整 $\\kappa \\in [0,1]$ 参数，可以控制转弯的紧凑程度和路径特性。
      * $\\kappa=0$：生成的平滑路径将通过（或非常接近）中间航点 $w\_i$。
      * $\\kappa=1$：生成一个单一圆弧的、最短过渡时间的转弯路径。
      * $0 \< \\kappa \< 1$：生成由三个连续圆弧构成的平滑路径。
  * **几何精确**：严格按照论文中描述的几何关系（如图4和图9）进行计算。
  * **离散化输出**：输出构成平滑转弯路径的一系列离散二维点，方便集成到仿真系统或实际飞行控制器中。
  * **模块化设计**：代码结构清晰，易于理解和集成到现有的无人机仿真或控制系统中。

## 文件结构 

  * `kappa_trajectory_generator.h` (或 `.cpp`): 包含 `Point2D`, `TurnGeometry` 结构体和 `KappaTrajectoryGenerator` 类的定义与实现。
  * `main.cpp` (示例): 展示了如何使用 `KappaTrajectoryGenerator` 类来生成轨迹。

## 如何使用 

1.  **包含代码**: 将 `kappa_trajectory_generator.h` 和 `kappa_trajectory_generator.cpp` (如果分离的话) 添加到您的 C++ 项目中。
2.  **创建对象**:
    ```cpp
    #include "kappa_trajectory_generator.h" // (或者您的头文件名)

    KappaTrajectoryGenerator generator;
    ```
3.  **准备输入参数**:
      * `w_im1`, `w_i`, `w_ip1`: 三个连续的 `Point2D` 航点。
      * `v_hat`: 期望飞行速度 (m/s)。
      * `turn_rate_c`: 最大转向速率 (rad/s)。
      * `kappa`: $\\kappa$ 参数 (0.0 到 1.0)。
      * `num_points_per_arc`: 每个圆弧段生成的离散点数量（可选，默认为20）。
4.  **生成轨迹**:
    ```cpp
    Point2D w0 = {0, 0};
    Point2D w1 = {10, 0};
    Point2D w2 = {10, 10};

    double v_hat = 5.0;
    double c_rate = 0.5; // R = v_hat / c_rate = 10m
    double kappa_value = 0.5;

    TurnGeometry turn_geom = generator.generateTurnTrajectory(w0, w1, w2, v_hat, c_rate, kappa_value);

    if (turn_geom.is_valid) {
        // 使用 turn_geom.tangent_point_on_incoming_segment (T1)
        // 使用 turn_geom.tangent_point_on_outgoing_segment (T2)
        // 使用 turn_geom.path_points (连接T1和T2的圆弧点序列)
        std::cout << "转弯信息: " << turn_geom.info << std::endl;
        std::cout << "T1: (" << turn_geom.tangent_point_on_incoming_segment.x << ", "
                  << turn_geom.tangent_point_on_incoming_segment.y << ")" << std::endl;
        // ...
    } else {
        std::cerr << "错误: 无法生成轨迹 - " << turn_geom.info << std::endl;
    }
    ```
5.  **路径拼接**:
    生成的 `TurnGeometry` 对象包含了转弯的入切点 `T1`、出切点 `T2` 以及连接它们的平滑曲线点 `path_points`。在完整的路径规划中，您需要：
      * 将前一段直线路径的终点连接到当前转弯的 `T1` 点。
      * 使用 `path_points` 作为转弯部分的路径。
      * 将当前转弯的 `T2` 点连接到下一段直线路径的起点。

## 理论基础 

本代码的实现基于以下论文中描述的 $\\kappa$-轨迹几何构造方法：

  * Anderson, E. P., Beard, R. W., & McLain, T. W. (2003). *Real-time dynamic trajectory smoothing for unmanned aerial vehicles*. (Submitted to IEEE Transactions on Control Systems Technology, February, 2003).

具体参考论文的第三章 "EXTREMAL TRAJECTORIES" (特别是 Definition 3, Figure 4) 和第六章 "SIMULATION" (Figure 9 的几何细节在 Lemma 6 的证明中)。

## 编译示例 (g++) 

```bash
g++ -std=c++17 kappa_trajectory_generator.cpp main.cpp -o kappa_path_test -lm
./kappa_path_test
```

(假设所有代码都在一个 `kappa_trajectory_generator.cpp` 文件中，或者相应地链接其他 `.cpp` 文件。)

## 注意事项 

  * 代码目前仅支持二维平面内的航点和轨迹生成。
  * 对于非常接近直线或180度U型转弯的极端情况，$\\kappa$-轨迹模型可能会退化。代码包含对这些情况的初步处理，但用户应根据具体应用场景进行验证和调整。
  * 生成的离散点密度由 `num_points_per_arc` 参数控制。

