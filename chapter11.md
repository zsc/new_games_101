# 第11章：动画与物理模拟

在计算机图形学中，静态的美丽固然重要，但赋予场景以生命的动画才是真正让虚拟世界活起来的关键。本章将深入探讨动画的数学基础，从简单的关键帧插值到复杂的物理模拟，涵盖质点系统、刚体动力学以及流体模拟的核心概念。我们将特别关注数值方法的稳定性和准确性，这对于创建可信的动画效果至关重要。

## 11.1 动画基本概念

### 11.1.1 时间与运动的数学表示

动画的本质是物体属性随时间的连续变化。在数学上，我们可以将任何可动画的属性表示为时间的函数：

$$\mathbf{p}(t) = f(t)$$

其中 $\mathbf{p}$ 可以是位置、旋转、颜色或任何其他属性。对于位置动画，我们有：

$$\mathbf{x}(t) = [x(t), y(t), z(t)]^T$$

速度和加速度通过微分获得：

$$\mathbf{v}(t) = \frac{d\mathbf{x}}{dt}, \quad \mathbf{a}(t) = \frac{d^2\mathbf{x}}{dt^2}$$

### 11.1.2 关键帧动画与插值

关键帧动画是最基础也是最常用的动画技术。给定一系列关键帧 $\{(\mathbf{p}_i, t_i)\}_{i=0}^n$，我们需要在这些关键帧之间进行插值。

**线性插值**：
$$\mathbf{p}(t) = \mathbf{p}_i + \frac{t - t_i}{t_{i+1} - t_i}(\mathbf{p}_{i+1} - \mathbf{p}_i), \quad t \in [t_i, t_{i+1}]$$

**球面线性插值（用于旋转）**：
对于四元数 $\mathbf{q}_i$ 和 $\mathbf{q}_{i+1}$：

$$\text{slerp}(\mathbf{q}_i, \mathbf{q}_{i+1}, s) = \frac{\sin((1-s)\theta)}{\sin\theta}\mathbf{q}_i + \frac{\sin(s\theta)}{\sin\theta}\mathbf{q}_{i+1}$$

其中 $\cos\theta = \mathbf{q}_i \cdot \mathbf{q}_{i+1}$，$s = \frac{t - t_i}{t_{i+1} - t_i}$。

### 11.1.3 样条曲线在动画中的应用

为了获得更平滑的动画，我们常使用样条曲线。Catmull-Rom样条是动画中的常用选择，因为它通过所有控制点：

$$\mathbf{p}(t) = \mathbf{T}(t) \mathbf{M}_{CR} \mathbf{G}$$

其中：
$$\mathbf{T}(t) = [1, t, t^2, t^3]$$

$$\mathbf{M}_{CR} = \frac{1}{2}\begin{bmatrix}
0 & 2 & 0 & 0 \\
-1 & 0 & 1 & 0 \\
2 & -5 & 4 & -1 \\
-1 & 3 & -3 & 1
\end{bmatrix}$$

$$\mathbf{G} = [\mathbf{p}_{i-1}, \mathbf{p}_i, \mathbf{p}_{i+1}, \mathbf{p}_{i+2}]^T$$

### 11.1.4 动画压缩与优化

实际应用中，动画数据可能非常庞大。常用的压缩技术包括：

1. **关键帧约简**：使用Douglas-Peucker算法移除冗余关键帧
2. **曲线拟合**：用低阶多项式或样条拟合原始动画
3. **主成分分析（PCA）**：对于多个相关的动画通道

误差度量通常使用：
$$E = \sum_{t} \|\mathbf{p}_{original}(t) - \mathbf{p}_{compressed}(t)\|^2$$

## 11.2 质点弹簧系统与运动学

### 11.2.1 质点系统的物理模型

质点系统是物理模拟的基础。每个质点 $i$ 具有：
- 质量 $m_i$
- 位置 $\mathbf{x}_i$
- 速度 $\mathbf{v}_i$

系统的运动由牛顿第二定律控制：
$$m_i \mathbf{a}_i = \mathbf{f}_i = \sum_j \mathbf{f}_{ij} + \mathbf{f}_{ext,i}$$

### 11.2.2 弹簧-阻尼器模型

弹簧力的计算遵循胡克定律：
$$\mathbf{f}_{spring} = -k_s (|\mathbf{x}_i - \mathbf{x}_j| - l_0) \frac{\mathbf{x}_i - \mathbf{x}_j}{|\mathbf{x}_i - \mathbf{x}_j|}$$

阻尼力用于能量耗散：
$$\mathbf{f}_{damping} = -k_d (\mathbf{v}_i - \mathbf{v}_j) \cdot \frac{\mathbf{x}_i - \mathbf{x}_j}{|\mathbf{x}_i - \mathbf{x}_j|} \frac{\mathbf{x}_i - \mathbf{x}_j}{|\mathbf{x}_i - \mathbf{x}_j|}$$

组合后的总力：
$$\mathbf{f}_{total} = \mathbf{f}_{spring} + \mathbf{f}_{damping}$$

### 11.2.3 约束与求解器

约束可以表示为：
$$C(\mathbf{x}) = 0$$

对于距离约束：
$$C_{ij} = |\mathbf{x}_i - \mathbf{x}_j| - d_{ij} = 0$$

使用拉格朗日乘数法，运动方程变为：
$$\mathbf{M}\ddot{\mathbf{x}} = \mathbf{f} + \mathbf{J}^T \boldsymbol{\lambda}$$
$$\mathbf{J}\ddot{\mathbf{x}} = -\dot{\mathbf{J}}\dot{\mathbf{x}}$$

其中 $\mathbf{J} = \frac{\partial C}{\partial \mathbf{x}}$ 是约束的雅可比矩阵。

### 11.2.4 逆运动学(IK)与正运动学(FK)

**正运动学**：给定关节角度 $\boldsymbol{\theta}$，计算末端执行器位置：
$$\mathbf{x}_{end} = f_{FK}(\boldsymbol{\theta})$$

**逆运动学**：给定目标位置 $\mathbf{x}_{target}$，求解关节角度：
$$\boldsymbol{\theta} = f_{IK}^{-1}(\mathbf{x}_{target})$$

IK通常没有解析解，需要迭代求解。使用雅可比转置法：
$$\Delta\boldsymbol{\theta} = \alpha \mathbf{J}^T (\mathbf{x}_{target} - \mathbf{x}_{current})$$

或伪逆法：
$$\Delta\boldsymbol{\theta} = \mathbf{J}^+ (\mathbf{x}_{target} - \mathbf{x}_{current})$$

其中 $\mathbf{J}^+ = \mathbf{J}^T(\mathbf{J}\mathbf{J}^T)^{-1}$ 是雅可比矩阵的伪逆。

### 11.2.5 雅可比矩阵与奇异性

雅可比矩阵定义为：
$$\mathbf{J} = \frac{\partial \mathbf{x}}{\partial \boldsymbol{\theta}} = \begin{bmatrix}
\frac{\partial x}{\partial \theta_1} & \cdots & \frac{\partial x}{\partial \theta_n} \\
\frac{\partial y}{\partial \theta_1} & \cdots & \frac{\partial y}{\partial \theta_n} \\
\frac{\partial z}{\partial \theta_1} & \cdots & \frac{\partial z}{\partial \theta_n}
\end{bmatrix}$$

当 $\det(\mathbf{J}\mathbf{J}^T) \approx 0$ 时，系统接近奇异配置。处理奇异性的方法：

1. **阻尼最小二乘（DLS）**：
$$\Delta\boldsymbol{\theta} = \mathbf{J}^T(\mathbf{J}\mathbf{J}^T + \lambda^2\mathbf{I})^{-1}\Delta\mathbf{x}$$

2. **奇异值分解（SVD）**：
$$\mathbf{J} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T$$
$$\mathbf{J}^+ = \mathbf{V}\boldsymbol{\Sigma}^+\mathbf{U}^T$$

## 11.3 常微分方程求解

### 11.3.1 物理模拟中的ODE问题

物理模拟通常归结为求解以下形式的ODE系统：
$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y})$$

其中状态向量 $\mathbf{y} = [\mathbf{x}, \mathbf{v}]^T$，导数为：
$$\frac{d\mathbf{y}}{dt} = \begin{bmatrix} \mathbf{v} \\ \mathbf{M}^{-1}\mathbf{f} \end{bmatrix}$$

### 11.3.2 显式与隐式积分方法

**显式欧拉法**：
$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\mathbf{f}(t_n, \mathbf{y}_n)$$

优点：简单快速
缺点：条件稳定，需要小时间步长

**隐式欧拉法**：
$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\mathbf{f}(t_{n+1}, \mathbf{y}_{n+1})$$

需要求解非线性方程，通常使用牛顿迭代：
$$(\mathbf{I} - h\frac{\partial \mathbf{f}}{\partial \mathbf{y}})\Delta\mathbf{y} = h\mathbf{f}(t_{n+1}, \mathbf{y}_n)$$

**中点法（RK2）**：
$$\mathbf{k}_1 = h\mathbf{f}(t_n, \mathbf{y}_n)$$
$$\mathbf{k}_2 = h\mathbf{f}(t_n + \frac{h}{2}, \mathbf{y}_n + \frac{\mathbf{k}_1}{2})$$
$$\mathbf{y}_{n+1} = \mathbf{y}_n + \mathbf{k}_2$$

**RK4（四阶Runge-Kutta）**：
$$\mathbf{k}_1 = h\mathbf{f}(t_n, \mathbf{y}_n)$$
$$\mathbf{k}_2 = h\mathbf{f}(t_n + \frac{h}{2}, \mathbf{y}_n + \frac{\mathbf{k}_1}{2})$$
$$\mathbf{k}_3 = h\mathbf{f}(t_n + \frac{h}{2}, \mathbf{y}_n + \frac{\mathbf{k}_2}{2})$$
$$\mathbf{k}_4 = h\mathbf{f}(t_n + h, \mathbf{y}_n + \mathbf{k}_3)$$
$$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{1}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)$$

### 11.3.3 稳定性分析

考虑测试方程 $\dot{y} = \lambda y$，其中 $\lambda < 0$。

显式欧拉的稳定条件：
$$|1 + h\lambda| < 1$$
$$h < \frac{2}{|\lambda|}$$

对于弹簧系统，$\lambda \approx -\sqrt{k/m}$，因此：
$$h < 2\sqrt{\frac{m}{k}}$$

隐式欧拉总是稳定的（A-稳定），但可能过度阻尼。

### 11.3.4 自适应时间步长

使用嵌入式RK方法（如RK45）估计误差：
$$\epsilon = \|\mathbf{y}_{n+1}^{(5)} - \mathbf{y}_{n+1}^{(4)}\|$$

调整时间步长：
$$h_{new} = h \cdot \min\left(2, \max\left(0.5, 0.9\left(\frac{\epsilon_{tol}}{\epsilon}\right)^{1/5}\right)\right)$$

### 11.3.5 能量守恒与数值耗散

理想的哈密顿系统应该保持总能量：
$$H = T + V = \frac{1}{2}\mathbf{v}^T\mathbf{M}\mathbf{v} + V(\mathbf{x})$$

辛积分器（如Verlet方法）更好地保持能量：
$$\mathbf{x}_{n+1} = 2\mathbf{x}_n - \mathbf{x}_{n-1} + h^2\mathbf{M}^{-1}\mathbf{f}_n$$

速度通过有限差分计算：
$$\mathbf{v}_n = \frac{\mathbf{x}_{n+1} - \mathbf{x}_{n-1}}{2h}$$

## 11.4 刚体与流体模拟

### 11.4.1 刚体动力学基础

刚体的状态由以下量描述：
- 质心位置：$\mathbf{x}$
- 方向（四元数）：$\mathbf{q}$
- 线速度：$\mathbf{v}$
- 角速度：$\boldsymbol{\omega}$

运动方程：
$$m\dot{\mathbf{v}} = \mathbf{f}$$
$$\mathbf{I}\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega}) = \boldsymbol{\tau}$$

其中 $\mathbf{I}$ 是惯性张量，$\boldsymbol{\tau}$ 是力矩。

世界坐标系下的惯性张量：
$$\mathbf{I}_{world} = \mathbf{R}\mathbf{I}_{body}\mathbf{R}^T$$

### 11.4.2 碰撞检测与响应

**碰撞检测**分为两个阶段：

1. **粗检测（Broad Phase）**：使用AABB、OBB或球形包围盒
2. **细检测（Narrow Phase）**：GJK算法或SAT（分离轴定理）

**碰撞响应**使用冲量方法：

接触点的相对速度：
$$\mathbf{v}_{rel} = (\mathbf{v}_A + \boldsymbol{\omega}_A \times \mathbf{r}_A) - (\mathbf{v}_B + \boldsymbol{\omega}_B \times \mathbf{r}_B)$$

冲量大小：
$$j = \frac{-(1 + e)\mathbf{v}_{rel} \cdot \mathbf{n}}{\frac{1}{m_A} + \frac{1}{m_B} + (\mathbf{I}_A^{-1}(\mathbf{r}_A \times \mathbf{n})) \times \mathbf{r}_A) \cdot \mathbf{n} + ((\mathbf{I}_B^{-1}(\mathbf{r}_B \times \mathbf{n})) \times \mathbf{r}_B) \cdot \mathbf{n}}$$

其中 $e$ 是恢复系数，$\mathbf{n}$ 是接触法线。

### 11.4.3 流体的欧拉与拉格朗日描述

**拉格朗日描述**：跟踪流体粒子
$$\frac{D\mathbf{u}}{Dt} = \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u}$$

**欧拉描述**：固定网格上的场
$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\nabla^2\mathbf{u} + \mathbf{g}$$

### 11.4.4 Navier-Stokes方程

不可压缩流体的NS方程：
$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\nabla^2\mathbf{u} + \mathbf{f}$$
$$\nabla \cdot \mathbf{u} = 0$$

使用投影法求解：
1. **对流步**：$\mathbf{u}^* = \mathbf{u}^n - \Delta t(\mathbf{u} \cdot \nabla)\mathbf{u}$
2. **扩散步**：$\mathbf{u}^{**} = \mathbf{u}^* + \Delta t\nu\nabla^2\mathbf{u}^*$
3. **投影步**：求解 $\nabla^2 p = \frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}^{**}$
4. **校正步**：$\mathbf{u}^{n+1} = \mathbf{u}^{**} - \frac{\Delta t}{\rho}\nabla p$

### 11.4.5 SPH与格子方法

**光滑粒子流体动力学（SPH）**：

密度估计：
$$\rho_i = \sum_j m_j W(|\mathbf{x}_i - \mathbf{x}_j|, h)$$

压力：
$$p_i = k(\rho_i - \rho_0)$$

加速度：
$$\mathbf{a}_i = -\sum_j m_j\left(\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}\right)\nabla W_{ij} + \nu\sum_j m_j\frac{\mathbf{u}_j - \mathbf{u}_i}{\rho_j}\nabla^2 W_{ij}$$

**格子Boltzmann方法（LBM）**：

分布函数演化：
$$f_i(\mathbf{x} + \mathbf{e}_i\Delta t, t + \Delta t) = f_i(\mathbf{x}, t) + \Omega_i(\mathbf{x}, t)$$

其中 $\Omega_i$ 是碰撞算子。

宏观量通过矩获得：
$$\rho = \sum_i f_i, \quad \rho\mathbf{u} = \sum_i f_i\mathbf{e}_i$$

## 本章小结

本章深入探讨了计算机图形学中的动画与物理模拟技术：

1. **动画基础**：
   - 时间参数化：$\mathbf{p}(t) = f(t)$
   - 关键帧插值：线性插值、球面线性插值（slerp）
   - 样条曲线：Catmull-Rom样条 $\mathbf{p}(t) = \mathbf{T}(t)\mathbf{M}_{CR}\mathbf{G}$

2. **质点弹簧系统**：
   - 运动方程：$m_i\mathbf{a}_i = \sum_j \mathbf{f}_{ij} + \mathbf{f}_{ext,i}$
   - 弹簧力：$\mathbf{f} = -k_s(|\mathbf{x}| - l_0)\hat{\mathbf{x}}$
   - IK雅可比：$\Delta\boldsymbol{\theta} = \mathbf{J}^+\Delta\mathbf{x}$

3. **数值积分**：
   - 显式欧拉：$\mathbf{y}_{n+1} = \mathbf{y}_n + h\mathbf{f}(t_n, \mathbf{y}_n)$
   - 隐式欧拉：$\mathbf{y}_{n+1} = \mathbf{y}_n + h\mathbf{f}(t_{n+1}, \mathbf{y}_{n+1})$
   - 稳定性条件：$h < 2\sqrt{m/k}$（显式）
   - Verlet积分：$\mathbf{x}_{n+1} = 2\mathbf{x}_n - \mathbf{x}_{n-1} + h^2\mathbf{M}^{-1}\mathbf{f}_n$

4. **刚体动力学**：
   - 欧拉方程：$\mathbf{I}\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega}) = \boldsymbol{\tau}$
   - 碰撞冲量：$j = \frac{-(1+e)\mathbf{v}_{rel} \cdot \mathbf{n}}{1/m_A + 1/m_B + ...}$

5. **流体模拟**：
   - Navier-Stokes：$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\nabla^2\mathbf{u}$
   - 不可压缩条件：$\nabla \cdot \mathbf{u} = 0$
   - SPH密度：$\rho_i = \sum_j m_j W(|\mathbf{x}_i - \mathbf{x}_j|, h)$

## 练习题

### 基础题

**练习11.1** 给定两个关键帧位置 $\mathbf{p}_0 = [0, 0, 0]^T$ 在 $t_0 = 0$，$\mathbf{p}_1 = [1, 2, 0]^T$ 在 $t_1 = 1$。使用线性插值计算 $t = 0.3$ 时的位置。

<details>
<summary>答案</summary>

使用线性插值公式：
$$\mathbf{p}(0.3) = \mathbf{p}_0 + \frac{0.3 - 0}{1 - 0}(\mathbf{p}_1 - \mathbf{p}_0) = [0, 0, 0]^T + 0.3([1, 2, 0]^T - [0, 0, 0]^T)$$
$$= [0, 0, 0]^T + 0.3[1, 2, 0]^T = [0.3, 0.6, 0]^T$$

</details>

**练习11.2** 一个质量为 $m = 1$ kg 的质点通过刚度 $k = 100$ N/m 的弹簧连接到原点。使用显式欧拉法，计算保证稳定性的最大时间步长。

<details>
<summary>答案</summary>

稳定性条件：$h < 2\sqrt{m/k}$

代入数值：
$$h < 2\sqrt{1/100} = 2 \times 0.1 = 0.2 \text{ s}$$

因此最大时间步长约为 0.2 秒。

</details>

**练习11.3** 两个刚体发生完全弹性碰撞（$e = 1$），质量分别为 $m_A = 2$ kg 和 $m_B = 3$ kg，碰撞前速度为 $v_A = 5$ m/s 和 $v_B = -2$ m/s（一维情况）。计算碰撞后的速度。

<details>
<summary>答案</summary>

动量守恒：$m_A v_A + m_B v_B = m_A v_A' + m_B v_B'$
能量守恒（完全弹性）：$\frac{1}{2}m_A v_A^2 + \frac{1}{2}m_B v_B^2 = \frac{1}{2}m_A v_A'^2 + \frac{1}{2}m_B v_B'^2$

代入数值：
$2 \times 5 + 3 \times (-2) = 2v_A' + 3v_B'$
$10 - 6 = 2v_A' + 3v_B'$
$4 = 2v_A' + 3v_B'$ ... (1)

相对速度关系（$e = 1$）：
$v_B' - v_A' = -(v_B - v_A) = -(-2 - 5) = 7$ ... (2)

解方程组得：
$v_A' = -2.6$ m/s
$v_B' = 4.4$ m/s

</details>

### 挑战题

**练习11.4** 推导二维平面上三连杆机械臂的雅可比矩阵。设三个关节角度为 $\theta_1, \theta_2, \theta_3$，每段长度为 $l_1, l_2, l_3$。

*提示：使用链式法则和旋转矩阵。*

<details>
<summary>答案</summary>

末端执行器位置：
$$\mathbf{x} = \begin{bmatrix} x \\ y \end{bmatrix} = \begin{bmatrix}
l_1\cos\theta_1 + l_2\cos(\theta_1 + \theta_2) + l_3\cos(\theta_1 + \theta_2 + \theta_3) \\
l_1\sin\theta_1 + l_2\sin(\theta_1 + \theta_2) + l_3\sin(\theta_1 + \theta_2 + \theta_3)
\end{bmatrix}$$

雅可比矩阵：
$$\mathbf{J} = \begin{bmatrix}
-l_1\sin\theta_1 - l_2\sin(\theta_1 + \theta_2) - l_3\sin(\theta_1 + \theta_2 + \theta_3) & -l_2\sin(\theta_1 + \theta_2) - l_3\sin(\theta_1 + \theta_2 + \theta_3) & -l_3\sin(\theta_1 + \theta_2 + \theta_3) \\
l_1\cos\theta_1 + l_2\cos(\theta_1 + \theta_2) + l_3\cos(\theta_1 + \theta_2 + \theta_3) & l_2\cos(\theta_1 + \theta_2) + l_3\cos(\theta_1 + \theta_2 + \theta_3) & l_3\cos(\theta_1 + \theta_2 + \theta_3)
\end{bmatrix}$$

</details>

**练习11.5** 证明Verlet积分保持二阶精度，并说明为什么它比显式欧拉法更好地保持能量。

*提示：使用泰勒展开分析局部截断误差。*

<details>
<summary>答案</summary>

泰勒展开：
$$\mathbf{x}(t + h) = \mathbf{x}(t) + h\dot{\mathbf{x}}(t) + \frac{h^2}{2}\ddot{\mathbf{x}}(t) + \frac{h^3}{6}\dddot{\mathbf{x}}(t) + O(h^4)$$
$$\mathbf{x}(t - h) = \mathbf{x}(t) - h\dot{\mathbf{x}}(t) + \frac{h^2}{2}\ddot{\mathbf{x}}(t) - \frac{h^3}{6}\dddot{\mathbf{x}}(t) + O(h^4)$$

相加得：
$$\mathbf{x}(t + h) = 2\mathbf{x}(t) - \mathbf{x}(t - h) + h^2\ddot{\mathbf{x}}(t) + O(h^4)$$

局部截断误差为 $O(h^4)$，因此是二阶精度。

Verlet积分是时间可逆的（symplectic），这保证了相空间体积守恒，从而更好地保持长期能量行为。

</details>

**练习11.6** 设计一个自适应时间步长算法，用于模拟刚度变化很大的弹簧系统。要求在保持稳定性的同时最大化效率。

*提示：考虑使用局部误差估计和刚度检测。*

<details>
<summary>答案</summary>

算法框架：

1. 估计局部刚度：$k_{max} = \max_i k_i$
2. 计算稳定性限制：$h_{stable} = C\sqrt{m_{min}/k_{max}}$，其中 $C < 2$
3. 使用嵌入式RK45估计误差：$\epsilon = \|\mathbf{x}_5 - \mathbf{x}_4\|$
4. 调整步长：
   $$h_{new} = \min\left(h_{stable}, h \cdot \left(\frac{\epsilon_{tol}}{\epsilon}\right)^{1/5}\right)$$
5. 对于接近奇异的配置，切换到隐式方法

关键点是同时考虑稳定性和精度要求。

</details>

**练习11.7** 推导SPH方法中的压力梯度项，并解释为什么使用对称形式可以保证动量守恒。

*提示：从连续性方程开始，考虑核函数的性质。*

<details>
<summary>答案</summary>

从压力梯度开始：
$$-\frac{1}{\rho}\nabla p$$

SPH近似：
$$-\frac{1}{\rho_i}\nabla p_i = -\frac{1}{\rho_i}\sum_j m_j \frac{p_j}{\rho_j}\nabla W_{ij}$$

但这不保证动量守恒。对称形式：
$$\mathbf{f}_i^{pressure} = -m_i\sum_j m_j\left(\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}\right)\nabla W_{ij}$$

由于 $\nabla W_{ij} = -\nabla W_{ji}$，有：
$$\mathbf{f}_i + \mathbf{f}_j = 0$$

这保证了动量守恒。

</details>

**练习11.8** 分析投影法求解不可压缩流体时的数值误差来源，并提出改进方案。

*提示：考虑分裂误差和边界条件处理。*

<details>
<summary>答案</summary>

误差来源：
1. **分裂误差**：将NS方程分步求解引入 $O(\Delta t)$ 误差
2. **投影误差**：压力Poisson方程的离散化误差
3. **边界条件**：速度和压力边界条件的不一致性

改进方案：
1. 使用高阶分裂格式（如Strang splitting）
2. 交错网格（MAC grid）减少压力震荡
3. 一致的边界条件处理
4. 使用增量压力格式：
   $$\nabla^2 \phi = \frac{1}{\Delta t}\nabla \cdot \mathbf{u}^*$$
   $$p^{n+1} = p^n + \phi - \frac{\mu}{2}\nabla \cdot \mathbf{u}^*$$

</details>

## 常见陷阱与错误

### 1. 数值积分的稳定性问题

**陷阱**：使用显式欧拉法模拟刚性系统（如高刚度弹簧）时出现爆炸。

**解决方案**：
- 检查时间步长：确保 $h < 2\sqrt{m/k}$
- 对刚性系统使用隐式方法
- 实现自适应时间步长

### 2. 四元数插值的错误

**陷阱**：直接对四元数分量进行线性插值导致非单位四元数。

**解决方案**：
- 始终使用球面线性插值（slerp）
- 插值后归一化四元数
- 注意处理 $\theta \approx 0$ 的特殊情况

### 3. 约束求解的数值漂移

**陷阱**：位置约束在迭代过程中逐渐违反。

**解决方案**：
- 使用Baumgarte稳定化：$\mathbf{J}\ddot{\mathbf{x}} = -\alpha C - \beta \dot{C}$
- 实现位置校正步骤
- 使用更稳定的约束公式

### 4. 碰撞检测的遗漏

**陷阱**：快速移动的物体穿透（tunneling）。

**解决方案**：
- 实现连续碰撞检测（CCD）
- 使用保守前进算法
- 限制最大速度或细分时间步

### 5. 流体模拟的质量损失

**陷阱**：不可压缩流体出现体积变化。

**解决方案**：
- 确保散度自由：$\nabla \cdot \mathbf{u} = 0$
- 使用高精度压力求解器
- 实现FLIP/PIC混合方法减少数值耗散

### 6. SPH的粒子聚集

**陷阱**：SPH粒子倾向于聚集成团。

**解决方案**：
- 添加人工压力项
- 使用XSPH速度平滑
- 调整核函数和搜索半径

## 最佳实践检查清单

### 动画系统设计

- [ ] 选择合适的时间表示（固定时间步 vs 可变时间步）
- [ ] 实现插值方法的缓存机制
- [ ] 支持动画混合和过渡
- [ ] 考虑动画压缩需求
- [ ] 实现动画事件系统

### 物理模拟架构

- [ ] 分离碰撞检测和动力学求解
- [ ] 实现多种积分器供选择
- [ ] 支持约束的增量添加/删除
- [ ] 设计清晰的力/约束接口
- [ ] 实现调试可视化

### 数值方法选择

- [ ] 根据刚度选择积分方法
- [ ] 实现误差估计和自适应
- [ ] 考虑能量守恒需求
- [ ] 平衡精度和性能
- [ ] 提供稳定性监控

### 性能优化

- [ ] 使用空间分割加速碰撞检测
- [ ] 实现多级时间步进
- [ ] 利用SIMD指令
- [ ] 考虑GPU加速
- [ ] 实现LOD系统

### 鲁棒性保证

- [ ] 处理数值退化情况
- [ ] 实现故障恢复机制
- [ ] 添加参数有效性检查
- [ ] 提供稳定性警告
- [ ] 记录诊断信息

### 用户接口

- [ ] 提供直观的参数调节
- [ ] 实现实时预览
- [ ] 支持动画导入/导出
- [ ] 提供性能分析工具
- [ ] 文档化所有参数含义