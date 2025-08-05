# 第12章：物质点法

物质点法（Material Point Method, MPM）是一种混合欧拉-拉格朗日数值方法，在计算机图形学中被广泛用于模拟各种材料的大变形行为，包括雪、沙、泥浆、弹性体等。MPM结合了粒子方法的灵活性和网格方法的计算效率，特别适合处理拓扑变化、碎裂、混合等复杂现象。本章将深入探讨MPM的理论基础、移动最小二乘MPM（MLS-MPM）的改进，以及高级技术与应用。

## 章节大纲

### 12.1 物质点法基础理论
- 连续介质力学基础
- MPM的基本框架
- 粒子-网格交互
- 时间积分方案

### 12.2 MLS-MPM（移动最小二乘物质点法）
- 传统MPM的局限性
- MLS-MPM的理论推导
- APIC与MLS-MPM的关系
- 实现细节与优化

### 12.3 高级MPM技术与应用
- 多材料MPM
- 自适应MPM
- GPU加速技术
- 工程应用案例

## 12.1 物质点法基础理论

### 12.1.1 连续介质力学基础

在介绍MPM之前，我们需要回顾连续介质力学的基本概念。考虑一个连续体 $\Omega$，其运动可以用映射 $\phi: \Omega_0 \times [0,T] \rightarrow \mathbb{R}^d$ 描述，其中 $\Omega_0$ 是初始构型，$d$ 是空间维度。

**变形梯度**定义为：
$$\mathbf{F} = \frac{\partial \phi}{\partial \mathbf{X}}$$

其中 $\mathbf{X}$ 是材料坐标。变形梯度 $\mathbf{F}$ 包含了所有的变形信息，包括旋转和拉伸。

**应变张量**常用的有：
- Green应变张量：$\mathbf{E} = \frac{1}{2}(\mathbf{F}^T\mathbf{F} - \mathbf{I})$
- 小应变张量：$\boldsymbol{\epsilon} = \frac{1}{2}(\nabla\mathbf{u} + \nabla\mathbf{u}^T)$

**柯西应力张量** $\boldsymbol{\sigma}$ 描述了当前构型下的应力状态。

**守恒方程**：
1. 质量守恒：$\frac{D\rho}{Dt} + \rho\nabla \cdot \mathbf{v} = 0$
2. 动量守恒：$\rho\frac{D\mathbf{v}}{Dt} = \nabla \cdot \boldsymbol{\sigma} + \rho\mathbf{b}$
3. 能量守恒：$\rho\frac{De}{Dt} = \boldsymbol{\sigma} : \mathbf{D} - \nabla \cdot \mathbf{q}$

其中 $\mathbf{D} = \frac{1}{2}(\nabla\mathbf{v} + \nabla\mathbf{v}^T)$ 是变形率张量。

### 12.1.2 MPM的基本框架

MPM将连续体离散为一组物质点（粒子），每个粒子携带质量、动量、变形梯度等物理量。同时使用背景欧拉网格进行动量方程的求解。

**基本变量**：
- 粒子变量：位置 $\mathbf{x}_p$，速度 $\mathbf{v}_p$，质量 $m_p$，体积 $V_p$，变形梯度 $\mathbf{F}_p$
- 网格变量：速度 $\mathbf{v}_i$，质量 $m_i$

**MPM的计算循环**包含以下步骤：

1. **粒子到网格的传输（P2G）**：
   $$m_i = \sum_p w_{ip} m_p$$
   $$m_i \mathbf{v}_i = \sum_p w_{ip} m_p \mathbf{v}_p$$
   
   其中 $w_{ip} = N_i(\mathbf{x}_p)$ 是插值权重，$N_i$ 是形函数。

2. **网格上的动量更新**：
   $$m_i \mathbf{v}_i^{n+1} = m_i \mathbf{v}_i^n + \Delta t \mathbf{f}_i$$
   
   其中力 $\mathbf{f}_i$ 包含内力和外力：
   $$\mathbf{f}_i = -\sum_p V_p \boldsymbol{\sigma}_p \nabla w_{ip} + \mathbf{f}_i^{\text{ext}}$$

3. **网格到粒子的传输（G2P）**：
   $$\mathbf{v}_p^{n+1} = \sum_i w_{ip} \mathbf{v}_i^{n+1}$$
   $$\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \Delta t \mathbf{v}_p^{n+1}$$

4. **变形梯度更新**：
   $$\mathbf{F}_p^{n+1} = (\mathbf{I} + \Delta t \sum_i \mathbf{v}_i^{n+1} \otimes \nabla w_{ip}) \mathbf{F}_p^n$$

### 12.1.3 粒子-网格交互

**形函数选择**是MPM的关键。常用的形函数包括：

1. **线性形函数**（三线性插值）：
   $$N(\xi) = \prod_{d=1}^D N^1(\xi_d)$$
   其中 $N^1(\xi) = 1 - |\xi|$ 当 $|\xi| < 1$，否则为0。

2. **二次B样条**：
   $$N^2(\xi) = \begin{cases}
   \frac{3}{4} - \xi^2 & |\xi| < \frac{1}{2} \\
   \frac{1}{2}(\frac{3}{2} - |\xi|)^2 & \frac{1}{2} \leq |\xi| < \frac{3}{2} \\
   0 & \text{otherwise}
   \end{cases}$$

3. **三次B样条**提供更高的连续性，但计算成本更高。

**权重计算**涉及粒子位置到网格节点的映射：
$$w_{ip} = N\left(\frac{\mathbf{x}_p - \mathbf{x}_i}{h}\right)$$
其中 $h$ 是网格间距。

### 12.1.4 时间积分方案

MPM中常用的时间积分方案包括：

1. **显式积分**：
   - 前向欧拉：简单但稳定性受限
   - 辛欧拉：保持能量守恒性质
   - Runge-Kutta方法：高阶精度

2. **半隐式积分**：
   - 对某些项（如弹性力）进行隐式处理
   - 提高稳定性，允许更大的时间步长

**CFL条件**：
$$\Delta t \leq C \frac{h}{c}$$
其中 $c = \sqrt{E/\rho}$ 是材料中的声速，$C$ 是CFL数（通常取0.1-0.5）。

### 12.1.5 本构模型

MPM可以处理各种材料模型：

1. **线弹性**：
   $$\boldsymbol{\sigma} = \lambda \text{tr}(\boldsymbol{\epsilon})\mathbf{I} + 2\mu\boldsymbol{\epsilon}$$

2. **超弹性**（Neo-Hookean）：
   $$\boldsymbol{\sigma} = \mu(\mathbf{B} - \mathbf{I}) + \lambda\ln(J)\mathbf{I}$$
   其中 $\mathbf{B} = \mathbf{F}\mathbf{F}^T$，$J = \det(\mathbf{F})$。

3. **塑性模型**：
   - von Mises屈服准则
   - Drucker-Prager模型（适用于颗粒材料）
   - 雪的本构模型

**应力计算**通常涉及：
1. 计算应变或变形梯度
2. 应用本构关系获得应力
3. 考虑塑性修正（如果适用）

## 12.2 MLS-MPM（移动最小二乘物质点法）

### 12.2.1 传统MPM的局限性

传统MPM虽然在处理大变形问题上有优势，但存在一些固有问题：

1. **数值耗散**：粒子-网格-粒子的传输过程引入人工粘性
2. **角动量守恒**：FLIP（Fluid Implicit Particle）方法保持角动量但引入噪声
3. **计算效率**：需要存储和更新大量粒子信息
4. **精度问题**：低阶形函数限制了应力计算的精度

### 12.2.2 MLS-MPM的理论推导

MLS-MPM的核心思想是使用移动最小二乘（Moving Least Squares）方法在粒子周围重构速度场，从而获得更精确的变形梯度更新。

**速度场重构**：
假设在粒子 $p$ 附近的速度场可以表示为：
$$\mathbf{v}(\mathbf{x}) = \mathbf{c}_0 + \mathbf{C}_1 (\mathbf{x} - \mathbf{x}_p)$$

其中 $\mathbf{c}_0$ 是常数项，$\mathbf{C}_1$ 是一阶导数矩阵。

**最小二乘拟合**：
通过最小化加权误差：
$$E = \sum_i w_{ip} m_i ||\mathbf{v}_i - \mathbf{c}_0 - \mathbf{C}_1(\mathbf{x}_i - \mathbf{x}_p)||^2$$

求解得到：
$$\mathbf{c}_0 = \sum_i w_{ip} \mathbf{v}_i$$
$$\mathbf{C}_1 = \left(\sum_i w_{ip} \mathbf{v}_i \otimes (\mathbf{x}_i - \mathbf{x}_p)\right) \mathbf{D}_p^{-1}$$

其中 $\mathbf{D}_p = \sum_i w_{ip} (\mathbf{x}_i - \mathbf{x}_p) \otimes (\mathbf{x}_i - \mathbf{x}_p)$ 是惯性张量。

**APIC矩阵**：
MLS-MPM引入仿射粒子-网格（APIC）矩阵 $\mathbf{C}_p$ 来存储速度梯度信息：
$$\mathbf{C}_p = \mathbf{B}_p \mathbf{D}_p^{-1}$$
其中 $\mathbf{B}_p = \sum_i w_{ip} \mathbf{v}_i \otimes (\mathbf{x}_i - \mathbf{x}_p)$。

### 12.2.3 MLS-MPM算法流程

1. **初始化**：
   - 设置粒子位置 $\mathbf{x}_p$、速度 $\mathbf{v}_p$、质量 $m_p$
   - 初始化变形梯度 $\mathbf{F}_p = \mathbf{I}$
   - 计算初始 APIC 矩阵 $\mathbf{C}_p = \mathbf{0}$

2. **P2G传输（改进版）**：
   $$m_i = \sum_p w_{ip} m_p$$
   $$m_i \mathbf{v}_i = \sum_p w_{ip} m_p [\mathbf{v}_p + \mathbf{C}_p(\mathbf{x}_i - \mathbf{x}_p)]$$

3. **网格动量更新**（与传统MPM相同）：
   $$\mathbf{v}_i^{n+1} = \mathbf{v}_i^n + \frac{\Delta t}{m_i} \mathbf{f}_i$$

4. **G2P传输（改进版）**：
   $$\mathbf{v}_p^{n+1} = \sum_i w_{ip} \mathbf{v}_i^{n+1}$$
   $$\mathbf{C}_p^{n+1} = \sum_i w_{ip} \mathbf{v}_i^{n+1} \otimes (\mathbf{x}_i - \mathbf{x}_p) \cdot \mathbf{D}_p^{-1}$$

5. **粒子更新**：
   $$\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \Delta t \mathbf{v}_p^{n+1}$$
   $$\mathbf{F}_p^{n+1} = (\mathbf{I} + \Delta t \mathbf{C}_p^{n+1}) \mathbf{F}_p^n$$

### 12.2.4 优化技巧

1. **惯性张量的高效计算**：
   对于规则网格和对称核函数：
   $$\mathbf{D}_p = \frac{4h^2}{3}\mathbf{I}$$（二维情况）
   $$\mathbf{D}_p = \frac{h^2}{2}\mathbf{I}$$（三维情况）

2. **体积守恒**：
   通过修正变形梯度保持体积：
   $$\mathbf{F}_p^{\text{corrected}} = J^{1/d} \mathbf{F}_p / \det(\mathbf{F}_p)^{1/d}$$

3. **数值稳定性增强**：
   - 限制变形梯度的奇异值
   - 使用隐式积分处理刚性材料
   - 自适应时间步长

### 12.2.5 与APIC的关系

MLS-MPM可以看作是APIC（Affine Particle-In-Cell）方法的推广：

1. **APIC**：保存仿射速度场 $\mathbf{v}(\mathbf{x}) = \mathbf{v}_p + \mathbf{C}_p(\mathbf{x} - \mathbf{x}_p)$
2. **MLS-MPM**：通过MLS重构获得最优的 $\mathbf{C}_p$
3. **理论联系**：MLS-MPM在数学上等价于使用特定权重的APIC

**优势对比**：
- APIC：实现简单，但精度受限
- MLS-MPM：更高精度，更好的动量守恒
- 计算成本相近

## 12.3 高级MPM技术与应用

### 12.3.1 多材料MPM

实际应用中经常需要模拟多种材料的相互作用，如流体-固体耦合、多相流等。

**多材料表示方法**：

1. **分离粒子法**：
   - 不同材料使用不同的粒子集
   - 在网格上进行耦合
   - 适用于界面清晰的情况

2. **混合粒子法**：
   - 单个粒子可以包含多种材料
   - 使用体积分数 $\phi_k$ 表示材料 $k$ 的占比
   - 适用于混合材料

**界面处理**：

1. **接触算法**：
   $$\mathbf{v}_i^{k,\text{final}} = \mathbf{v}_i^k + \frac{\Delta t}{m_i^k} \mathbf{f}_i^{\text{contact}}$$
   
   接触力通过罚函数或拉格朗日乘子计算。

2. **混合理论**：
   应力通过体积加权平均：
   $$\boldsymbol{\sigma} = \sum_k \phi_k \boldsymbol{\sigma}_k$$

**表面张力**：
对于流体界面，需要考虑表面张力：
$$\mathbf{f}_{\text{surface}} = \gamma \kappa \mathbf{n}$$
其中 $\gamma$ 是表面张力系数，$\kappa$ 是曲率，$\mathbf{n}$ 是法向量。

### 12.3.2 自适应MPM

为了提高计算效率和精度，可以使用自适应技术：

1. **空间自适应**：
   - **自适应网格细化**（AMR）：在需要高精度的区域使用更细的网格
   - **粒子分裂与合并**：根据变形程度调整粒子密度
   
   分裂准则：
   $$\text{split if } \det(\mathbf{F}_p) > \theta_{\text{split}}$$
   
   合并准则：
   $$\text{merge if } ||\mathbf{x}_p - \mathbf{x}_q|| < r_{\text{merge}} \text{ and } \det(\mathbf{F}_p) < \theta_{\text{merge}}$$

2. **时间自适应**：
   根据局部CFL条件调整时间步长：
   $$\Delta t_p = \min\left(C_{\text{CFL}} \frac{h}{||\mathbf{v}_p|| + c_p}, \Delta t_{\text{max}}\right)$$

3. **模型自适应**：
   - 在小变形区域使用线性模型
   - 在大变形区域使用非线性模型
   - 动态切换本构模型

### 12.3.3 GPU加速技术

MPM天然适合GPU并行化，关键优化策略包括：

1. **数据结构优化**：
   - **粒子排序**：按空间位置排序减少内存访问冲突
   - **网格哈希**：使用空间哈希表快速查找粒子-网格映射
   - **稀疏网格**：只分配活跃网格节点

2. **并行化策略**：
   - **P2G阶段**：每个粒子独立计算，使用原子操作累加到网格
   - **网格更新**：每个网格节点独立更新
   - **G2P阶段**：每个粒子独立更新

3. **性能优化**：
   - **共享内存使用**：缓存频繁访问的网格数据
   - **Warp级别优化**：确保线程束内的一致性
   - **混合精度计算**：在不影响精度的地方使用float16

**典型加速比**：
- CPU单线程 → GPU：100-1000倍
- 实时模拟成为可能（>30 FPS）

### 12.3.4 工程应用案例

1. **计算机动画**：
   - 雪崩模拟（迪士尼《冰雪奇缘》）
   - 沙尘效果
   - 流体-固体交互

2. **地质工程**：
   - 滑坡预测
   - 土壤液化分析
   - 爆破效果模拟

3. **生物力学**：
   - 软组织变形
   - 血流动力学
   - 细胞力学

4. **工业设计**：
   - 3D打印过程模拟
   - 材料成型分析
   - 碰撞安全测试

### 12.3.5 前沿研究方向

1. **机器学习增强**：
   - 使用神经网络学习本构模型
   - 加速结构预测
   - 参数自动调优

2. **多尺度MPM**：
   - 宏观-介观-微观耦合
   - 分子动力学与连续介质耦合
   - 自适应尺度转换

3. **拓扑优化**：
   - 基于MPM的形状优化
   - 材料分布优化
   - 多目标优化

4. **实时交互系统**：
   - VR/AR应用
   - 触觉反馈
   - 实时物理编辑
