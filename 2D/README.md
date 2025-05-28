# 二维平面上的航点等长轨迹的κ求解



本代码适用于处理二维平面上的航点，该 C++ 代码本身并不直接计算出三维的 Kappa **曲线**的点集。

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

如需考虑航点为三维航点，可以参照第III节的 C++实现代码（三维航点）。

此 C++ 代码提供了一个模块化的稳健起点，请务必在仿真环境中对多种航点配置进行充分测试。
