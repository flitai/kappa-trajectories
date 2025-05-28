# 三维航点等长轨迹的$\kappa$求解C++实现

本代码主要实现内容包括：

1.  接收三维航点 $\mathbf{w}_{i-1}, \mathbf{w}_i, \mathbf{w}_{i+1}$ 作为输入。
2.  基于这些三维航点计算转向角 $\beta$。这个角度是在由这三个点定义的平面内的。
3.  计算出一个标量参数 $\kappa^*$。这个 $\kappa^*$ 的目标是使得后续由 DTS（Dynamically Time-Scaled）算法生成的实际 $\kappa$-轨迹的**路径长度**与原始分段直线路径的长度一致。

这个算法的核心是**求解参数 $\kappa^*$**。

之后，这个计算出来的 $\kappa^*$ 会被传递给 DTS 控制律，由 DTS 系统利用这个 $\kappa^*$ 来**生成实际的平滑转弯路径**。

**如果原始航点是三维的，那么 DTS 系统生成的这个 $\kappa$-轨迹也将在三维空间中**。每个转弯部分（由三个连续航点定义）将在包含这三个点的平面内进行平滑。由于连续的转弯平面可能不同，整个轨迹将是三维的。

功能要点总结：

* 该 C++ 代码：输入三维航点，输出一个用于长度匹配的**标量参数 $\kappa^*$**。
* DTS 系统（未在此代码中实现）：输入 $\kappa^*$ 和航点等信息，输出实际的**三维平滑曲线**。



This code will define a class `KappaTrajectorySolver` that encapsulates the logic for calculating $\kappa^*$


**关键组件和特性：**

1. **`Point3D` 结构体**：一个用于 3D 点/向量的基础结构，包含减法、归一化和点积等常用操作。如果你已有现成的向量数学库（例如 Eigen、GLM），可以直接替换使用。
2. **`KappaTrajectorySolver` 类**：

   * **`findKappaStar` 方法**：主要公共接口接口，接受三个航路点（$\mathbf{w}_{i-1}, \mathbf{w}_i, \mathbf{w}_{i+1}$）、速度 $\hat v$、转向速率 $c$ 以及精度 $\varepsilon$ 作为输入。方法内部首先计算 $R$ 和 $\beta$，然后执行二分查找。
   * **`calculateLambda` 方法**：私有辅助函数，利用“引理6”中的公式计算 $\Lambda(\kappa)$。
   * **数值稳定性**：

     * 将用于 `acos` 的点积结果限制在` [-1,1]` 之间。
     * 将 $\Xi$ 限幅到 $[\delta,\,1-\delta]$，以避免 $\Xi$ 接近 0 或 1 时导致的 $acos(\Xi)$ 和 $\sqrt{1/\Xi^2-1}$ 不稳定。
     * 对 $\Xi\sqrt{1/\Xi^2-1}$ 在 $\Xi\approx1$ 时保证正确返回 0。
     * 增加了对 $\beta$ 极接近 0 或 $\pi$ 的检查：

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
3. **调用 `findKappaStar`**：当需为由三点定义的转弯段计算 $\kappa^*$ 时：

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

* **角度 $\beta$ 定义**：代码中

  $$
    \beta = \cos^{-1}\bigl(\hat{\mathbf{t}}_{i,i+1}^T \hat{\mathbf{t}}_{i-1,i}\bigr)
  $$

  即航点 $\mathbf{w}_{i-1}\to\mathbf{w}_i$ 与 $\mathbf{w}_i\to\mathbf{w}_{i+1}$ 两向量间的夹角；$\beta=0$ 表示直行，$\beta=\pi$ 表示 180° 的 U 型转弯。
* **$\beta=0$ 奇异性**：$\Lambda(\kappa)$ 中的 $\cot(\beta/2)$ 在 $\beta=0$ 时无定义。理论（引理7、定理8）假设 $\beta\in(0,\pi)$ 才保证 $\Lambda(0)>0$、$\Lambda(1)<0$。实现中对 $\beta\approx0$ 添加了退化处理：直接返回 $\kappa^*=0$。对 $\beta\approx\pi$ 也做了检查，若有需要可根据 U 转的期望效果调整。
* **数值鲁棒性**：对 $\Xi$ 的限幅以及对临界值的处理非常关键；可根据使用的浮点类型（double 或 float）调整 $\delta$ 大小。
* **错误处理**：代码中包含对零除、端点条件不满足等情况的异常抛出，请将其纳入仿真软件的错误管理流程。

此 C++ 代码提供了一个模块化的稳健起点，请务必在仿真环境中对多种航点配置进行充分测试。
