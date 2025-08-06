# 第8章：全局光照

全局光照（Global Illumination）是计算机图形学中最具挑战性的问题之一。与局部光照不同，全局光照考虑了光线在场景中的多次反射、折射和散射，能够生成逼真的视觉效果，如软阴影、颜色渗透（color bleeding）、焦散（caustics）等。本章将从辐射度量学的基础出发，推导渲染方程，并介绍蒙特卡洛方法在求解全局光照中的应用。

## 8.1 辐射度量学

辐射度量学（Radiometry）是描述光能传播的物理框架。理解这些概念对于正确实现基于物理的渲染至关重要。

### 8.1.1 基本概念

辐射度量学中的物理量形成了一个层次结构，从能量开始，逐步引入时间、空间和方向的微分。理解这些量之间的关系对于正确实现基于物理的渲染至关重要。

**辐射能量（Radiant Energy）**
- 符号：$Q$
- 单位：焦耳（J）
- 物理意义：光子携带的总能量
- 与光子的关系：$Q = \sum_i h\nu_i$，其中 $h$ 是普朗克常数，$\nu_i$ 是第 $i$ 个光子的频率
- 量子力学视角：每个光子能量 $E = h\nu = \frac{hc}{\lambda}$
- 波长与颜色：可见光波长范围 380-780nm，对应不同颜色感知

**辐射功率/辐射通量（Radiant Power/Flux）**
- 符号：$\Phi = \frac{dQ}{dt}$
- 单位：瓦特（W）或流明（lm）
- 物理意义：单位时间内通过某个表面或空间区域的能量
- 实例：100W灯泡每秒辐射100焦耳的能量
- 光度学对应：光通量（Luminous Flux），考虑人眼响应曲线
- 光效转换：典型白炽灯约15 lm/W，LED可达100+ lm/W
- 积分形式：$\Phi = \int_A \int_{\Omega} L(\mathbf{p}, \omega) \cos\theta \, d\omega \, dA$

**辐射强度（Radiant Intensity）**
- 符号：$I = \frac{d\Phi}{d\omega}$
- 单位：W/sr（瓦特/球面度）
- 物理意义：点光源在某个方向上单位立体角的功率
- 适用条件：光源尺寸远小于到观察点的距离（$r \gg d_{source}$）
- 与点光源的关系：对于各向同性点光源，$I = \frac{\Phi}{4\pi}$
- 实际应用：描述LED、远距离恒星等近似点光源
- 方向分布：实际光源常有方向性，如聚光灯 $I(\theta) = I_0 \cos^n\theta$
- 距离无关性：强度不随距离变化（但照度会以 $1/r^2$ 衰减）

**辐照度（Irradiance）**
- 符号：$E = \frac{d\Phi}{dA}$
- 单位：W/m²
- 物理意义：单位面积接收到的功率（所有方向的总和）
- 与入射方向的关系：$E = \int_{\Omega} L_i(\omega) \cos\theta \, d\omega$
- 实际测量：照度计测量的就是辐照度的光度学对应量
- 太阳常数：地球大气层外的太阳辐照度约为 1361 W/m²
- 点光源辐照度：$E = \frac{I}{r^2} \cos\theta$（平方反比定律）
- 环境光近似：常用常数辐照度 $E_{ambient}$ 模拟天空光
- 半球辐照度：$E = \int_{\Omega^+} L_i \cos\theta \, d\omega$（只考虑上半球）

**辐射出射度（Radiant Exitance）**
- 符号：$M = \frac{d\Phi}{dA}$
- 单位：W/m²
- 物理意义：单位面积发出的功率
- 与辐照度的区别：方向相反，一个描述入射，一个描述出射
- 与辐射度的关系：$M = \int_{\Omega} L_o(\omega) \cos\theta \, d\omega$
- Lambertian发射体：$M = \pi L$（各向同性发射）
- Stefan-Boltzmann定律：黑体 $M = \sigma T^4$
- 面光源建模：用出射度描述发光表面

**辐射度（Radiance）**
- 符号：$L = \frac{d^2\Phi}{dA \cos\theta \, d\omega}$
- 单位：W/(m²·sr)
- 物理意义：单位面积、单位立体角、单位投影面积上的功率
- 核心地位：辐射度是渲染方程的基本量，因为：
  1. 沿光线传播不变（在真空中）
  2. 相机传感器响应与辐射度成正比
  3. 可以通过积分得到其他所有辐射度量
  4. 人眼感知亮度直接对应辐射度
- 与其他量的关系：
  - $d\Phi = L \, dA \cos\theta \, d\omega$
  - $dE = L \cos\theta \, d\omega$
  - $dI = L \, dA \cos\theta$
- 光线束（Ray Bundle）：$d\Phi = L \, dA_{\perp} \, d\omega$
- 典型值：太阳表面 $\approx 2 \times 10^7$ W/(m²·sr)

### 8.1.2 辐射度的性质

辐射度是渲染中最重要的量，具有以下关键性质：

1. **方向性与五维函数**
   - 辐射度是位置和方向的函数：$L(\mathbf{p}, \omega)$ 或 $L(x, y, z, \theta, \phi)$
   - 这是一个五维函数（3D位置 + 2D方向）
   - 光场（Light Field）：固定 $z$，得到4D函数 $L(x, y, \theta, \phi)$
   - 环境贴图：固定位置，得到2D函数 $L(\theta, \phi)$
   - 全光函数（Plenoptic Function）：$L(x, y, z, \theta, \phi, \lambda, t)$ 包含波长和时间
   - 光场相机：通过微透镜阵列捕捉 4D 光场信息

2. **传播不变性（Radiance Invariance）**
   在真空或均匀介质中传播时，辐射度沿光线保持不变：
   $$L(\mathbf{p}, \omega) = L(\mathbf{p} + t\omega, \omega), \quad \forall t > 0$$
   
   **证明思路**：考虑两个微分面元之间的能量传输
   - 从面元1发出：$d\Phi_1 = L_1 \, dA_1 \cos\theta_1 \, d\omega_1$
   - 到达面元2：$d\Phi_2 = L_2 \, dA_2 \cos\theta_2 \, d\omega_2$
   - 由于 $d\omega_1 = \frac{dA_2 \cos\theta_2}{r^2}$ 和 $d\omega_2 = \frac{dA_1 \cos\theta_1}{r^2}$
   - 能量守恒要求：$d\Phi_1 = d\Phi_2$，因此 $L_1 = L_2$
   
   **几何光学含义**：
   - 光线是辐射度的特征线
   - 这个性质是射线跟踪算法的理论基础
   - 在介质界面处会发生变化（折射定律）

3. **线性叠加性**
   - 多个光源的贡献可以线性叠加：$L_{total} = \sum_i L_i$
   - 这是因为光子之间不相互作用（在线性光学范围内）
   - 使得我们可以独立计算每个光源的贡献
   - 非线性效应：高强度激光、双光子吸收等特殊情况
   - 波动光学效应：干涉、衍射需要考虑相位

4. **与相机响应的关系**
   - 相机传感器的响应正比于入射辐射度
   - 像素值 $\propto \int_{\text{pixel}} \int_{\text{aperture}} L \cos\theta \, dA \, d\omega \, dt$
   - 这解释了为什么辐射度是渲染的核心量
   - 曝光时间积分：$E_{sensor} = \int_0^T L(t) \, dt$
   - 景深效应：孔径积分产生模糊
   - HDR成像：扩展辐射度动态范围

5. **能量维度分析**
   - 检查单位的一致性：$[L] = \frac{[\Phi]}{[A][\omega]} = \frac{W}{m^2 \cdot sr}$
   - 积分恢复功率：$\Phi = \int_A \int_{\Omega} L \cos\theta \, dA \, d\omega$
   - 单位立体角的物理含义：1 sr ≈ 65.5°的圆锥角
   - 全球立体角：$4\pi$ sr，半球 $2\pi$ sr
   - 小立体角近似：$\omega \approx \frac{A}{r^2}$ 当 $A \ll r^2$

6. **参与介质中的衰减**
   在参与介质（如雾、烟）中，辐射度会衰减：
   $$\frac{dL}{ds} = -\sigma_t L + \sigma_s L_{scattered} + L_{emitted}$$
   其中 $\sigma_t$ 是消光系数，$\sigma_s$ 是散射系数
   
   **传输方程的解**：
   $$L(s) = L(0) e^{-\tau(s)} + \int_0^s \sigma_s L_{in}(s') e^{-\tau(s,s')} ds'$$
   其中光学厚度 $\tau(s) = \int_0^s \sigma_t(s') ds'$
   
   **相函数**：$p(\omega_i \to \omega_o)$ 描述散射方向分布
   - Rayleigh散射：小粒子，蓝天效应
   - Mie散射：大粒子，云雾效应

### 8.1.3 立体角与投影立体角

立体角是三维空间中角度的推广，理解它对于正确处理方向积分至关重要。

**立体角的几何意义**
- 定义：从原点看某个表面张成的"角度大小"
- 单位：球面度（steradian, sr）
- 整个球面的立体角：$4\pi$ sr
- 半球的立体角：$2\pi$ sr
- 计算公式：$\omega = \frac{A}{r^2}$（球面上面积$A$除以半径平方）
- 平面角类比：2D中弧度 $\theta = \frac{s}{r}$，3D中立体角 $\omega = \frac{A}{r^2}$
- 直观理解：圆锥角 $\omega = 2\pi(1 - \cos\theta)$，其中 $\theta$ 是半角
- 单位立体角的物理尺度：1 sr ≈ 65.5°圆锥角，约占全球面的 $\frac{1}{4\pi} \approx 8\%$

**立体角的拓扑结构**
- 球面 $S^2$ 是二维流形，局部同胚于 $\mathbb{R}^2$
- 参数化奇点：球坐标在极点处退化（$\sin\theta = 0$）
- 替代参数化：立方体映射、八面体映射避免奇点
- Gauss-Bonnet定理：$\int_{S^2} K \, dA = 4\pi$（$K=1/r^2$ 是高斯曲率）

**立体角的微分形式**
在球坐标系 $(\theta, \phi)$ 中：
$$d\omega = \sin\theta \, d\theta \, d\phi$$

推导：球面上的微分面元
- 沿 $\theta$ 方向：$r \, d\theta$
- 沿 $\phi$ 方向：$r \sin\theta \, d\phi$（纬度圈半径为 $r\sin\theta$）
- 面积：$dA = r^2 \sin\theta \, d\theta \, d\phi$
- 立体角：$d\omega = \frac{dA}{r^2} = \sin\theta \, d\theta \, d\phi$

**度量张量视角**
球面上的度量张量：
$$g = \begin{pmatrix} 1 & 0 \\ 0 & \sin^2\theta \end{pmatrix}$$

面积元素：$dA = \sqrt{\det(g)} \, d\theta \, d\phi = \sin\theta \, d\theta \, d\phi$

**笛卡尔坐标表示**：
对于方向向量 $\omega = (\omega_x, \omega_y, \omega_z)$：
- $\omega_x = \sin\theta\cos\phi$
- $\omega_y = \sin\theta\sin\phi$  
- $\omega_z = \cos\theta$
- 约束：$\omega_x^2 + \omega_y^2 + \omega_z^2 = 1$

**坐标转换的数值稳定性**：
- 避免 $\arccos$ 的数值误差：使用 $\theta = \text{atan2}(\sqrt{\omega_x^2 + \omega_y^2}, \omega_z)$
- 处理 $\phi$ 的奇点：当 $\sin\theta \approx 0$ 时，$\phi$ 不确定
- 四元数表示：避免万向锁问题，用于旋转插值
- Rodrigues公式：旋转向量 $\mathbf{v}$ 绕轴 $\mathbf{k}$ 角度 $\theta$：
  $$\mathbf{v}_{rot} = \mathbf{v}\cos\theta + (\mathbf{k} \times \mathbf{v})\sin\theta + \mathbf{k}(\mathbf{k} \cdot \mathbf{v})(1-\cos\theta)$$

**半球积分**
$$\int_{\Omega} d\omega = \int_0^{2\pi} \int_0^{\pi/2} \sin\theta \, d\theta \, d\phi = 2\pi \int_0^{\pi/2} \sin\theta \, d\theta = 2\pi[-\cos\theta]_0^{\pi/2} = 2\pi$$

**常见积分公式**：
- 全球：$\int_{4\pi} d\omega = 4\pi$
- 圆锥：$\int_{\theta_0} d\omega = 2\pi(1 - \cos\theta_0)$
- 小立体角近似：$\Delta\omega \approx \sin\theta \Delta\theta \Delta\phi$
- 球冠（spherical cap）：高度 $h$ 的球冠立体角 $\omega = 2\pi h/r$
- 椭圆锥：需要椭圆积分，无闭式解
- 任意多边形：使用 Girard's theorem，$\omega = \sum_i \alpha_i - (n-2)\pi$

**特殊函数与立体角**：
- Legendre多项式展开：$f(\theta) = \sum_{l=0}^{\infty} a_l P_l(\cos\theta)$
- 球谐函数：$Y_l^m(\theta, \phi)$ 构成 $L^2(S^2)$ 的正交基
- Zonal harmonics：旋转对称函数 $f(\theta) = \sum_l a_l P_l(\cos\theta)$
- 快速球谐变换：$O(N^2\log N)$ 算法

**投影立体角（Projected Solid Angle）**
$$d\omega^\perp = \cos\theta \, d\omega = \cos\theta \sin\theta \, d\theta \, d\phi$$

物理意义：
- 考虑了表面法线与入射方向的夹角
- Lambert余弦定律的体现
- 垂直入射的贡献最大，掠射角贡献趋于零
- 几何解释：在法线方向投影的立体角

**投影立体角的不变量**：
- 形状因子（Form Factor）的互易性：$F_{ij} A_i = F_{ji} A_j$
- Nusselt类比：投影立体角 = 半球上的"有效面积"
- 守恒定律：封闭系统中 $\sum_j F_{ij} = 1$

**微分几何视角**：
- 投影算子：$\mathcal{P}_\mathbf{n} = I - \mathbf{n}\mathbf{n}^T$
- 立体角的外微分：$d\omega = \frac{\mathbf{r} \cdot d\mathbf{A}}{r^3}$
- Stokes定理应用：$\oint_{\partial S} \mathbf{v} \cdot d\mathbf{l} = \int_S (\nabla \times \mathbf{v}) \cdot d\mathbf{A}$

**应用场景**：
- 辐照度计算：$E = \int_{\Omega} L \cos\theta \, d\omega = \int_{\Omega} L \, d\omega^\perp$
- 形状因子（Form Factor）计算
- 半球采样权重

**关键积分结果**
1. 半球上的投影立体角积分：
   $$\int_{\Omega} \cos\theta \, d\omega = \int_0^{2\pi} \int_0^{\pi/2} \cos\theta \sin\theta \, d\theta \, d\phi = \pi$$
   
   **推导过程**：
   $$\int_0^{\pi/2} \cos\theta \sin\theta \, d\theta = \int_0^{\pi/2} \frac{1}{2}\sin(2\theta) \, d\theta = \frac{1}{2}$$
   
   **几何意义**：半球在其底面的投影面积为 $\pi r^2$，除以 $r^2$ 得 $\pi$
   
   **变分原理**：最大化 $\int f(\theta) \cos\theta \, d\omega$ 的分布是 $\delta(\theta)$（法线方向）

2. 余弦的$n$次方积分（用于BRDF分析）：
   $$\int_{\Omega} \cos^n\theta \, d\omega = \frac{2\pi}{n+1}$$
   
   **推导**：使用Beta函数
   $$\int_0^{\pi/2} \cos^n\theta \sin\theta \, d\theta = \int_0^1 u^n \, du = \frac{1}{n+1}$$
   
   **特例**：
   - $n=0$: $2\pi$ （半球立体角）
   - $n=1$: $\pi$ （投影立体角）
   - $n=2$: $\frac{2\pi}{3}$ （Phong模型常用）
   - $n \to \infty$: $0$ （完美镜面）
   
   **一般化**：
   $$\int_{\Omega} \cos^m\theta \sin^n\theta \, d\omega = \frac{\pi B(\frac{m+1}{2}, \frac{n+1}{2})}{B(\frac{1}{2}, \frac{1}{2})}$$
   其中 $B$ 是Beta函数

3. 立体角与面积的转换：
   $$d\omega = \frac{\cos\theta' \, dA'}{||\mathbf{p} - \mathbf{p}'||^2}$$
   其中 $\theta'$ 是 $\mathbf{p}'$ 处表面法线与连线的夹角
   
   **双向形式**：
   $$dA \, d\omega = \frac{\cos\theta \cos\theta' \, dA \, dA'}{||\mathbf{p} - \mathbf{p}'||^2}$$
   
   **几何因子的完整形式**：
   $$G(\mathbf{p} \leftrightarrow \mathbf{p}') = \frac{\cos\theta \cos\theta'}{||\mathbf{p} - \mathbf{p}'||^2} V(\mathbf{p} \leftrightarrow \mathbf{p}')$$
   
   **退化情况处理**：
   - $\theta' > \pi/2$：背向，$G = 0$
   - $||\mathbf{p} - \mathbf{p}'|| \to 0$：需要正则化
   - 共面情况：$\cos\theta = 0$ 或 $\cos\theta' = 0$

**实际应用**
- **辐照度计算**：$E = \int_{\Omega} L_i \cos\theta \, d\omega$
- **环境贴图采样**：需要考虑 $\sin\theta$ 项进行重要性采样
  - 纬度-经度映射：$(u,v) \to (\theta,\phi) = (\pi v, 2\pi u)$
  - PDF转换：$p(\theta,\phi) = \frac{p(u,v)}{2\pi^2 \sin\theta}$
  - 立方体贴图：避免极点奇异性，6个面均匀采样
  - 八面体映射：更均匀的采样分布
- **面光源采样**：从立体角域转换到面积域进行采样
  - 均匀采样面积：$p_A = 1/A_{light}$
  - 转换到立体角：$p_\omega = p_A \cdot \frac{r^2}{\cos\theta'}$
  - 球面光源：直接在立体角域采样更高效
  - 线性变换光源（Linearly Transformed Cosines）：解析积分
- **半球均匀采样**：
  - 随机变量：$(\xi_1, \xi_2) \in [0,1]^2$
  - 转换：$\theta = \arccos(\xi_1)$, $\phi = 2\pi\xi_2$
  - PDF: $p(\omega) = 1/(2\pi)$
- **余弦加权采样**：
  - 转换：$\theta = \arcsin(\sqrt{\xi_1})$, $\phi = 2\pi\xi_2$
  - PDF: $p(\omega) = \cos\theta/\pi$
  - Malley's method：在单位圆盘采样后投影到半球

**高级采样技术**：
- **球面分层采样**：将球面划分为等面积区域
- **螺旋采样**：Fibonacci螺旋提供准均匀分布
- **球面Voronoi图**：最优传输理论应用
- **蓝噪声采样**：改善视觉质量，减少聚集

### 8.1.4 BRDF与反射方程

双向反射分布函数（Bidirectional Reflectance Distribution Function, BRDF）是描述表面反射特性的核心概念。

**BRDF的定义**
$$f_r(\omega_i \to \omega_o) = \frac{dL_o(\omega_o)}{dE(\omega_i)} = \frac{dL_o(\omega_o)}{L_i(\omega_i) \cos\theta_i \, d\omega_i}$$

物理解释：
- 入射辐照度引起的出射辐射度变化率
- 单位：$1/sr$（因为 $[L]/[E] = [W/(m^2 \cdot sr)]/[W/m^2] = 1/sr$）
- 描述了材质的反射特性，与几何形状无关
- 局部线性假设：输出正比于输入（在线性光学范围内）

**微观解释**：
- BRDF描述大量微观交互的统计平均
- 微表面模型提供物理基础
- 波长依赖性：$f_r(\omega_i \to \omega_o, \lambda)$
- 空间变化：Spatially Varying BRDF (SVBRDF): $f_r(\mathbf{x}, \omega_i \to \omega_o)$

**散射理论视角**：
- 光子与表面的散射截面：$\sigma(\omega_i \to \omega_o)$
- BRDF与散射截面的关系：$f_r = \frac{\sigma}{\cos\theta_i}$
- 量子力学解释：电子跃迁和光子重新发射

**一般化到BSDF**：
- BSDF = BRDF + BTDF（双向散射分布函数）
- BSSRDF：次表面散射，$S(\mathbf{x}_i, \omega_i \to \mathbf{x}_o, \omega_o)$
- 各向异性：$f_r(\theta_i, \phi_i \to \theta_o, \phi_o)$ 需要完整的4D函数

**BRDF的基本性质**

1. **非负性**：$f_r(\omega_i \to \omega_o) \geq 0$
   - 物理意义：能量不能为负
   - 数学结果：正定算子

2. **Helmholtz互易性（Reciprocity）**：
   $$f_r(\omega_i \to \omega_o) = f_r(\omega_o \to \omega_i)$$
   - 物理基础：时间反演对称性
   - 实际意义：光路可逆
   - 应用：双向路径追踪中路径权重计算
   - 注意：磁光效应下不成立（非互易材料）
   - 实验验证：使用测角光度计交换入射和观察方向

3. **能量守恒（Energy Conservation）**：
   $$\rho(\omega_i) = \int_{\Omega^+} f_r(\omega_i \to \omega_o) \cos\theta_o \, d\omega_o \leq 1$$
   - 反射的能量不能超过入射能量
   - 等号成立：完美反射（如理想镜面）
   - 小于1：部分能量被吸收
   - 方向半球反射率：$\rho(\omega_i)$ 依赖于入射方向
   - 半球-半球反射率：$\bar{\rho} = \frac{1}{\pi}\int_{\Omega^+} \rho(\omega_i) \cos\theta_i \, d\omega_i$
   
   **白炉测试**：
   在均匀入射照明下，出射辐射度应等于入射辐射度：
   $$L_o = \int_{\Omega^+} f_r \cdot L_i \cos\theta_i \, d\omega_i = L_i \int_{\Omega^+} f_r \cos\theta_i \, d\omega_i \leq L_i$$
   
   **广义互易性**：
   $$\rho(\omega_o) = \int_{\Omega^+} f_r(\omega_i \to \omega_o) \cos\theta_i \, d\omega_i$$
   由于Helmholtz互易性，$\rho(\omega_o) = \rho(\omega_i)$ 当交换角色

4. **线性性**：
   - BRDF对入射光强度是线性的
   - 允许叠加原理的应用
   - 多光源：$L_o = \sum_k \int f_r L_{i,k} \cos\theta_i \, d\omega_i$
   - 傅立叶分解：可以在频域分析BRDF

5. **各向同性 vs 各向异性**：
   - 各向同性：BRDF不依赖于旋转角 $\phi$
     - 简化：$f_r(\theta_i, \theta_o, |\phi_i - \phi_o|)$
     - 参数化：3D函数而非4D
   - 各向异性：如拉丝金属、织物、木材纹理
     - Ward模型：椭圆高斯分布
     - 切线空间参数化：使用局部坐标系

6. **正定性**：
   - 作为算子，BRDF必须是正定的
   - 保证积分方程的解存在且唯一

**常见BRDF模型**

1. **Lambertian（漫反射）**：
   $$f_r = \frac{\rho}{\pi}$$
   - $\rho$：反射率（albedo），$0 \leq \rho \leq 1$
   - 各向同性，与方向无关
   - 能量守恒验证：$\int_{\Omega} \frac{\rho}{\pi} \cos\theta \, d\omega = \rho$
   - 典型材质：石膏、粗糙纸张、未抛光的石头
   - 为什么除以$\pi$：保证 $\int \cos\theta \, d\omega = \pi$ 后总反射率为 $\rho$
   - 物理基础：完全粗糙表面，微表面法线均匀分布
   - 局限性：真实世界中没有完美的Lambertian表面

2. **理想镜面反射**：
   $$f_r(\omega_i \to \omega_o) = \frac{\delta(\omega_i - R(\omega_o))}{\cos\theta_i}$$
   - $R(\omega_o)$：反射方向 $R(\omega_o) = 2(\mathbf{n} \cdot \omega_o)\mathbf{n} - \omega_o$
   - 使用Dirac delta函数表示
   - 积分后：$L_o(\omega_o) = L_i(R(\omega_o))$
   - 实际实现：显式跟踪镜面反射，而非采样
   - Fresnel项：$F(\omega_i, \mathbf{n})$ 描述反射比例
   - 完整形式：$f_r = \frac{F(\omega_i, \mathbf{n})\delta(\omega_i - R(\omega_o))}{\cos\theta_i}$

3. **Phong模型**（经验模型）：
   $$f_r = \frac{k_d}{\pi} + k_s \frac{n+2}{2\pi} \cos^n\alpha$$
   - $\alpha$：反射方向与视线方向的夹角
   - $n$：光泽度指数（典型值：5-1000）
   - 注意：不满足能量守恒
   - 修正的Phong：使用半向量 $\mathbf{h}$ 代替反射向量
   - 归一化因子 $\frac{n+2}{2\pi}$ 保证积分为1
   
   **Blinn-Phong模型**：
   $$f_r = \frac{k_d}{\pi} + k_s \frac{n+8}{8\pi} (\mathbf{n} \cdot \mathbf{h})^n$$
   - 半向量：$\mathbf{h} = \frac{\omega_i + \omega_o}{||\omega_i + \omega_o||}$
   - 更物理正确，计算更高效
   - 与微表面模型的联系：$D(\mathbf{h}) \propto (\mathbf{n} \cdot \mathbf{h})^n$

4. **微表面模型（Microfacet）**：
   $$f_r = \frac{D(\mathbf{h}) G(\omega_i, \omega_o) F(\omega_i, \mathbf{h})}{4 \cos\theta_i \cos\theta_o}$$
   - $D$：法线分布函数（如GGX）
   - $G$：几何遮蔽函数
   - $F$：Fresnel项
   - $\mathbf{h}$：半向量 $\mathbf{h} = \frac{\omega_i + \omega_o}{||\omega_i + \omega_o||}$
   
   **法线分布函数$D$**：
   - GGX/Trowbridge-Reitz: $D(\mathbf{h}) = \frac{\alpha^2}{\pi((\mathbf{n} \cdot \mathbf{h})^2(\alpha^2 - 1) + 1)^2}$
     - 长尾分布，更真实的高光
     - 参数转换：$\alpha = \text{roughness}^2$
   - Beckmann: $D(\mathbf{h}) = \frac{1}{\pi\alpha^2\cos^4\theta_h} e^{-\frac{\tan^2\theta_h}{\alpha^2}}$
     - 基于高斯统计
     - 计算更复杂，但物理基础清晰
   - $\alpha$：粗糙度参数
   
   **几何遮蔽函数$G$**：
   - Smith G函数：$G(\omega_i, \omega_o) = G_1(\omega_i) G_1(\omega_o)$
   - $G_1$的GGX形式：$G_1(\omega) = \frac{2(\mathbf{n} \cdot \omega)}{(\mathbf{n} \cdot \omega) + \sqrt{\alpha^2 + (1-\alpha^2)(\mathbf{n} \cdot \omega)^2}}$
   - 高度相关形式：考虑遮蔽-阴影相关性
   
   **Fresnel项$F$**：
   - Schlick近似：$F(\omega, \mathbf{h}) = F_0 + (1-F_0)(1-\omega \cdot \mathbf{h})^5$
   - $F_0$：垂直入射时的反射率
   - 完整Fresnel方程：需要复折射率
   
   **为什么除以4**：
   - 来自雅可比行列式的转换
   - $\frac{\partial \omega_h}{\partial \omega_o} = \frac{1}{4(\omega_i \cdot \mathbf{h})}$
   - 保证从半向量分布到BRDF的正确映射

**反射方程（Reflection Equation）**
$$L_o(\mathbf{p}, \omega_o) = L_e(\mathbf{p}, \omega_o) + \int_{\Omega^+} f_r(\mathbf{p}, \omega_i \to \omega_o) L_i(\mathbf{p}, \omega_i) \cos\theta_i \, d\omega_i$$

各项含义：
- $L_e$：自发光项（emission）
- $L_i$：入射辐射度
- $\cos\theta_i$：Lambert余弦项，几何因子
- 积分域 $\Omega^+$：上半球

**分层材质模型**：
- 漫反射层 + 镜面层：$f_r = f_{diffuse} + f_{specular}$
- 能量分配：使用Fresnel项调节
- Disney BRDF：基于艺术家友好的参数
- 分层模型的互易性：每层单独满足

**BRDF的测量与存储**
- 测量设备：测角光度计（gonioreflectometer）
  - 机械臂系统：精确控制入射和观察角度
  - 图像基的测量：使用相机和球面光源
  - 时间成本：完整测量一个材质需要数小时
- 存储挑战：4D函数（2D入射 + 2D出射）
  - 离散化：$90 \times 90 \times 180 \times 360$ 网格
  - 内存需求：每个材质约100MB
- 压缩方法：
  - 球谐函数：$f_r \approx \sum_{l,m} a_{lm} Y_l^m$
  - 小波分解：多分辨率表示
  - 因子分解：$f_r \approx \sum_k u_k(\omega_i) v_k(\omega_o)$
  - 神经网络：学习压缩表示
- 数据库：
  - MERL BRDF：100种材质，各向同性
  - UTIA BTF：包含空间变化
  - RGL BRDF：包含各向异性材料

**BRDF的未来研究方向**：
- 波动光学BRDF：考虑干涉和衍射
- 时变BRDF：模拟老化和风化
- 参与介质BSSRDF：次表面散射
- 量子光学效应：超表面和纳米材料

## 8.2 渲染方程

### 8.2.1 渲染方程的推导

Kajiya在1986年提出的渲染方程统一了计算机图形学中的光照计算，是全局光照的数学基础。这个方程的美在于它用简洁的数学形式捕捉了光传输的复杂性。

**从局部到全局的演进**

1. **局部光照模型**（只考虑直接光照）：
   $$L_o = L_e + \sum_{\text{lights}} f_r \cdot L_{light} \cdot \cos\theta_i$$
   
   历史背景：
   - 1967年 Phong模型
   - 1975年 Blinn-Phong改进
   - 足以处理实时渲染，但缺乏真实感

2. **问题**：忽略了间接光照
   - 没有软阴影
   - 没有颜色渗透（color bleeding）
   - 没有镜面间接反射
   - 没有焦散（caustics）
   - 没有环境光遮蔽（ambient occlusion）

3. **关键洞察**：入射辐射度来自其他表面的出射辐射度
   $$L_i(\mathbf{p}, \omega_i) = L_o(\mathbf{p}', -\omega_i)$$
   
   其中：
   - $\mathbf{p}'$ 是射线 $(\mathbf{p}, \omega_i)$ 的第一个交点
   - 射线追踪函数：$\mathbf{p}' = \text{raycast}(\mathbf{p}, \omega_i)$
   - 负号表示方向反转（从 $\mathbf{p}'$ 指向 $\mathbf{p}$）
   
   **光线传播的拓扑性质**：
   - 光线是相空间中的流线
   - 在非均匀介质中弯曲（海市蜃楼效应）
   - 在引力场中弯曲（广义相对论）

**完整的渲染方程**

将入射辐射度的表达式代入反射方程：
$$L_o(\mathbf{p}, \omega_o) = L_e(\mathbf{p}, \omega_o) + \int_{\Omega^+} f_r(\mathbf{p}, \omega_i \to \omega_o) L_o(\mathbf{p}', -\omega_i) \cos\theta_i \, d\omega_i$$

这是一个**递归积分方程**：
- 未知量 $L_o$ 同时出现在方程两边
- 在不同位置和方向上相互耦合
- 解析求解几乎不可能

**数学性质**：
- 这是Fredholm第二类积分方程
- 解的存在性和唯一性由Banach不动点定理保证
- 收敛条件：平均反射率 < 1

**物理意义**：
- 表达了光的多次反射
- 统一了各种光照现象
- 是能量守恒的数学体现

**边界条件**
1. **真空/环境**：$L_o = L_{env}(\omega_o)$
   - 无限远处的环境光
   - HDR环境贴图
   - 天空模型（如Hosek-Wilkie）

2. **光源表面**：$L_e > 0$
   - 面光源：均匀或空间变化的发射
   - 温度辐射：Planck定律
   - IES光源数据

3. **非发光表面**：$L_e = 0$
   - 绝大多数物体属于此类
   - 但可能有荧光效应

4. **参与介质边界**：
   - 需要考虑体积散射
   - 辐射传输方程（RTE）

**简化形式**

定义反射算子 $\mathcal{T}$：
$$(\mathcal{T}L)(\mathbf{p}, \omega_o) = \int_{\Omega^+} f_r(\mathbf{p}, \omega_i \to \omega_o) L(\mathbf{p}', -\omega_i) \cos\theta_i \, d\omega_i$$

渲染方程变为：
$$L = L_e + \mathcal{T}L$$

这是Fredholm第二类积分方程的一个实例。

**历史上的其他形式**：
- **辐射度方程**（Radiosity）：漫反射表面的特例
- **体积渲染方程**：包含参与介质
- **波动方程**：考虑干涉和衍射

### 8.2.2 积分算子形式

渲染方程可以使用算子理论进行分析，这提供了深入的数学洞察和求解方法。这种形式化不仅优雅，而且为数值方法提供了理论基础。

**光传输算子的定义**

定义算子 $\mathcal{K}$：
$$(\mathcal{K}L)(\mathbf{p}, \omega) = \int_{\Omega^+} f_r(\mathbf{p}, \omega_i \to \omega) L(\mathbf{p}', -\omega_i) \cos\theta_i \, d\omega_i$$

这是一个线性算子：
- $\mathcal{K}(aL_1 + bL_2) = a\mathcal{K}L_1 + b\mathcal{K}L_2$
- 作用于函数空间 $L^2(\mathcal{M} \times S^2)$
- 是紧算子（compact operator）在适当条件下

**函数空间的选择**：
- $L^2$：平方可积，适合能量分析
- $L^\infty$：有界函数，适合最大值分析
- Sobolev空间：考虑光滑性

**算子方程**

渲染方程可写为：
$$L = L_e + \mathcal{K}L$$

或等价地：
$$(I - \mathcal{K})L = L_e$$

其中 $I$ 是恒等算子。

**Neumann级数解**

形式解为：
$$L = (I - \mathcal{K})^{-1}L_e = \sum_{n=0}^{\infty} \mathcal{K}^n L_e$$

收敛条件：谱半径 $\rho(\mathcal{K}) < 1$
- 物理意义：平均反射率小于1
- 保证能量最终衰减为零

**收敛速度分析**：
- 几何收敛：$||\mathcal{K}^n|| \leq \rho(\mathcal{K})^n$
- 实际场景：$\rho(\mathcal{K}) \approx 0.5-0.8$
- 3-5次迭代后误差 < 1%

**谱分解**：
- 特征函数：$\mathcal{K}\phi_i = \lambda_i \phi_i$
- 展开：$L = \sum_i \frac{\langle L_e, \phi_i \rangle}{1-\lambda_i} \phi_i$
- 主特征值决定收敛性

**物理解释**

每一项的含义：
- $\mathcal{K}^0 L_e = L_e$：直接看到光源
- $\mathcal{K}^1 L_e$：一次反射（直接光照）
- $\mathcal{K}^2 L_e$：二次反射（一次间接光照）
- $\mathcal{K}^n L_e$：$n$次反射

这也解释了为什么：
- 直接光照通常最亮
- 随着反射次数增加，贡献逐渐减小
- 大多数场景3-5次反射后可以忽略

**不同场景的反射次数**：
- 室外场景：2-3次通常足够
- 室内场景：5-10次可能需要
- 焦散效果：可能需要更多

**能量分布统计**：
- 直接光照：60-80%
- 一次间接：15-30%
- 高次间接：< 10%

**算子的核函数表示**

$\mathcal{K}$ 可以用核函数表示：
$$K(\mathbf{p} \to \mathbf{p}', \omega \to \omega') = f_r(\mathbf{p}, \omega' \to \omega) \cos\theta' G(\mathbf{p} \leftrightarrow \mathbf{p}') \delta(\omega' - \omega_{\mathbf{p} \to \mathbf{p}'})$$

其中：
- $G(\mathbf{p} \leftrightarrow \mathbf{p}')$：几何因子
- $\delta$：Dirac delta函数，确保方向一致

**核函数的性质**：
- 对称性（由于互易性）：$K(\mathbf{p} \to \mathbf{p}', \omega \to \omega') = K(\mathbf{p}' \to \mathbf{p}, \omega' \to \omega)$
- 正定性：$K \geq 0$
- 可积性：$\int \int K < \infty$ 在适当条件下

**离散化的核矩阵**：
- 有限元方法：$K_{ij} = \langle \phi_i, \mathcal{K}\phi_j \rangle$
- 矩阵大小：$N^2$，其中$N$是基函数个数
- 稀疏性：利用可见性和距离衰减

**迭代求解方法**

1. **固定点迭代**：
   $$L^{(k+1)} = L_e + \mathcal{K}L^{(k)}$$
   
2. **Jacobi迭代**：
   $$L^{(k+1)} = L_e + \mathcal{K}L^{(k)}$$
   
3. **Gauss-Seidel迭代**：
   更新时立即使用新值

**与有限元方法的联系**

离散化后，渲染方程变为线性系统：
$$(I - K)L = L_e$$

其中 $K$ 是离散化的传输矩阵。

### 8.2.3 路径积分形式

路径积分提供了另一种理解渲染方程的视角，将光线传输看作路径的集合。

**路径空间的定义**

一条长度为 $k$ 的路径：
$$\bar{x} = \mathbf{x}_0 \mathbf{x}_1 ... \mathbf{x}_k$$

其中：
- $\mathbf{x}_0$：相机位置
- $\mathbf{x}_k$：光源上的点
- $\mathbf{x}_1, ..., \mathbf{x}_{k-1}$：中间的反射点

路径空间：
$$\Omega = \bigcup_{k=2}^{\infty} \mathcal{M}^{k+1}$$

**路径贡献函数**

一条路径的贡献：
$$f(\bar{x}) = L_e(\mathbf{x}_k \to \mathbf{x}_{k-1}) \prod_{i=1}^{k-1} f_r(\mathbf{x}_{i+1} \to \mathbf{x}_i \to \mathbf{x}_{i-1}) G(\mathbf{x}_i \leftrightarrow \mathbf{x}_{i+1})$$

分解：
1. **发射项**：$L_e(\mathbf{x}_k \to \mathbf{x}_{k-1})$
2. **BRDF链**：$\prod f_r$，每个反射点的BRDF
3. **几何链**：$\prod G$，点之间的几何关系

**几何因子的详细分析**

$$G(\mathbf{x} \leftrightarrow \mathbf{y}) = \frac{\cos\theta_x \cos\theta_y}{||\mathbf{x} - \mathbf{y}||^2} V(\mathbf{x} \leftrightarrow \mathbf{y})$$

组成部分：
1. **余弦项**：$\cos\theta_x \cos\theta_y$
   - $\theta_x$：$\mathbf{x}$ 处法线与连线的夹角
   - 体现Lambert定律

2. **距离衰减**：$1/r^2$
   - 球面波的几何扩散
   - 能量守恒的体现

3. **可见性**：$V(\mathbf{x} \leftrightarrow \mathbf{y}) \in \{0, 1\}$
   - 1：两点互相可见
   - 0：存在遮挡

**渲染方程的路径积分形式**

$$L(\mathbf{x}_0 \to \mathbf{x}_1) = \sum_{k=2}^{\infty} \int_{\mathcal{M}^{k-1}} f(\mathbf{x}_0 \mathbf{x}_1 ... \mathbf{x}_k) \, dA(\mathbf{x}_2) ... dA(\mathbf{x}_k)$$

这是对所有可能路径的积分：
- 外层求和：对所有路径长度
- 内层积分：对所有路径配置

**路径分类（Heckbert记法）**

用正则表达式描述路径类型：
- L：光源
- E：相机（眼睛）
- S：镜面反射/折射
- D：漫反射

例子：
- L(D|S)*E：所有路径
- LDE：直接光照
- L(D|S)SE：焦散
- LS*DS*E：镜面间接光照

**测度论视角**

路径空间上的测度：
$$d\mu(\bar{x}) = dA(\mathbf{x}_0) ... dA(\mathbf{x}_k)$$

渲染的目标是计算：
$$I = \int_{\Omega} f(\bar{x}) \, d\mu(\bar{x})$$

这为蒙特卡洛方法提供了理论基础。

### 8.2.4 测度变换

在渲染中，经常需要在不同的积分域之间转换。理解测度变换对于正确实现重要性采样至关重要。

**立体角到面积的变换**

基本关系：
$$d\omega = \frac{\cos\theta' \, dA'}{||\mathbf{p} - \mathbf{p}'||^2}$$

推导：
1. 从 $\mathbf{p}$ 看 $dA'$ 张成的立体角
2. 投影面积：$dA' \cos\theta'$（$\theta'$ 是 $\mathbf{p}'$ 处的角度）
3. 立体角 = 投影面积 / 距离平方

**面积形式的渲染方程**

将测度变换应用于渲染方程：
$$L_o(\mathbf{p}, \omega_o) = L_e(\mathbf{p}, \omega_o) + \int_{\mathcal{M}} f_r(\mathbf{p}, \mathbf{p}' \to \omega_o) L_o(\mathbf{p}', \mathbf{p} - \mathbf{p}') G(\mathbf{p} \leftrightarrow \mathbf{p}') \, dA'$$

这里几何因子 $G$ 包含了测度变换的所有项。

**其他常见的测度变换**

1. **球面坐标到笛卡尔坐标**：
   $$d\omega = \sin\theta \, d\theta \, d\phi$$
   对应的向量：
   $$\omega = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$$

2. **半球均匀采样到余弦加权采样**：
   - 均匀采样：$p(\omega) = 1/(2\pi)$
   - 余弦加权：$p(\omega) = \cos\theta/\pi$
   - Jacobian：$\frac{\cos\theta/\pi}{1/(2\pi)} = 2\cos\theta$

3. **光源采样的测度变换**：
   - 在光源表面采样：$p_A(\mathbf{x}) = 1/A_{light}$
   - 转换到立体角：$p_\omega(\omega) = \frac{r^2}{A_{light}\cos\theta_{light}}$

**变换的雅可比行列式**

一般的变量替换公式：
$$p_Y(y) = p_X(x) \left|\frac{\partial x}{\partial y}\right|$$

多维情况：
$$p_Y(\mathbf{y}) = p_X(\mathbf{x}) \left|\det\left(\frac{\partial \mathbf{x}}{\partial \mathbf{y}}\right)\right|$$

**实际应用中的陷阱**

1. **忘记余弦项**：
   - 错误：$p_\omega = \frac{r^2}{A}$
   - 正确：$p_\omega = \frac{r^2}{A\cos\theta}$

2. **双重计算几何项**：
   - 错误：在BRDF和测度变换中都包含 $\cos\theta$
   - 正确：明确区分哪些项属于BRDF，哪些属于测度变换

3. **正负号错误**：
   - 注意法线方向和光线方向的关系
   - 使用 $|\cos\theta|$ 或确保角度在 $[0, \pi/2]$

## 8.3 蒙特卡洛积分与路径追踪

渲染方程是一个高维积分方程，解析求解几乎不可能。蒙特卡洛方法提供了数值求解的框架。

### 8.3.1 蒙特卡洛积分基础

对于积分 $I = \int_{\Omega} f(x) \, dx$，蒙特卡洛估计量为：

$$\langle I \rangle = \frac{1}{N} \sum_{i=1}^{N} \frac{f(X_i)}{p(X_i)}$$

其中 $X_i \sim p(x)$ 是从概率密度函数 $p(x)$ 采样的随机变量。

**期望与方差**：
- 期望：$E[\langle I \rangle] = I$（无偏估计）
- 方差：$V[\langle I \rangle] = \frac{1}{N} V\left[\frac{f(X)}{p(X)}\right]$

**收敛速度**：误差以 $O(1/\sqrt{N})$ 的速度收敛，与维度无关。

### 8.3.2 重要性采样

选择合适的概率密度函数 $p(x)$ 可以显著降低方差。理想情况下：
$$p(x) \propto |f(x)|$$

对于渲染，常用的采样策略包括：

1. **余弦加权采样**：$p(\omega) \propto \cos\theta$
2. **BRDF采样**：$p(\omega) \propto f_r(\omega_i \to \omega_o)$
3. **光源采样**：直接在光源上采样
4. **多重重要性采样（MIS）**：组合多种采样策略

### 8.3.3 路径追踪算法

基本路径追踪通过递归求解渲染方程：

```
radiance(ray):
    hit = intersect(ray, scene)
    if not hit:
        return background
    
    L_e = hit.emission
    
    // 直接光照
    L_d = sample_direct_lighting(hit)
    
    // 间接光照（俄罗斯轮盘赌）
    if random() < P_rr:
        wi, pdf = sample_brdf(hit)
        L_i = radiance(spawn_ray(hit, wi))
        L_indirect = f_r * L_i * cos(theta) / (pdf * P_rr)
    else:
        L_indirect = 0
    
    return L_e + L_d + L_indirect
```

### 8.3.4 俄罗斯轮盘赌

为了得到无偏估计同时避免无限递归，使用俄罗斯轮盘赌：

$$L = \begin{cases}
\frac{L_{continue}}{P_{rr}} & \text{以概率 } P_{rr} \\
0 & \text{以概率 } 1-P_{rr}
\end{cases}$$

期望值保持不变：$E[L] = P_{rr} \cdot \frac{L_{continue}}{P_{rr}} = L_{continue}$

### 8.3.5 直接光照采样

对于面光源，可以直接在光源表面采样：

$$L_d = \int_{A_{light}} f_r(\omega_i \to \omega_o) L_e \cos\theta_i \frac{\cos\theta_{light}}{||\mathbf{x} - \mathbf{x}_{light}||^2} V(\mathbf{x} \leftrightarrow \mathbf{x}_{light}) \, dA_{light}$$

采样策略：
1. 在光源表面均匀采样：$p(x_{light}) = 1/A_{light}$
2. 转换为立体角域的pdf：$p(\omega) = \frac{||\mathbf{x} - \mathbf{x}_{light}||^2}{A_{light} \cos\theta_{light}}$

### 8.3.6 多重重要性采样（MIS）

当有多种采样策略时，MIS提供了组合它们的最优方式。平衡启发式权重：

$$w_i(x) = \frac{p_i(x)}{\sum_j p_j(x)}$$

功率启发式（通常效果更好）：
$$w_i(x) = \frac{p_i^2(x)}{\sum_j p_j^2(x)}$$

组合估计量：
$$\langle I \rangle = \sum_i \frac{1}{n_i} \sum_{j=1}^{n_i} w_i(X_{ij}) \frac{f(X_{ij})}{p_i(X_{ij})}$$

### 8.3.7 路径空间的测度

在路径空间中，一条长度为 $k$ 的路径的概率密度为：

$$p(\bar{x}) = p(\mathbf{x}_0) \prod_{i=1}^{k} p(\mathbf{x}_i | \mathbf{x}_{i-1})$$

面积测度下的转换：
$$p(\mathbf{x}_i | \mathbf{x}_{i-1}) = p(\omega_i) \frac{\cos\theta_i}{||\mathbf{x}_i - \mathbf{x}_{i-1}||^2}$$

### 8.3.8 双向路径追踪

双向路径追踪（BDPT）从相机和光源同时追踪路径，然后连接它们：

1. 从相机追踪子路径：$\mathbf{x}_0, \mathbf{x}_1, ..., \mathbf{x}_s$
2. 从光源追踪子路径：$\mathbf{y}_0, \mathbf{y}_1, ..., \mathbf{y}_t$
3. 连接策略：连接 $\mathbf{x}_s$ 和 $\mathbf{y}_t$ 形成完整路径

对于每种连接策略 $(s,t)$，使用MIS组合所有可能的路径。

### 8.3.9 渐进式光子映射

光子映射是另一种全局光照算法，特别适合处理焦散等效果：

1. **光子发射阶段**：从光源发射光子，在场景中追踪并存储
2. **密度估计**：使用k近邻搜索估计辐射度
   $$L(\mathbf{x}, \omega) \approx \frac{1}{\pi r^2} \sum_{p \in N_k(\mathbf{x})} f_r(\omega_p \to \omega) \Phi_p$$

渐进式改进：
- 逐步缩小搜索半径：$r_{i+1} = r_i \sqrt{\frac{i}{i+1}}$
- 保证渐近无偏和一致性收敛

## 本章小结

**核心概念**：
1. 辐射度量学提供了描述光传输的物理框架
2. 渲染方程是全局光照的数学基础：
   $$L_o = L_e + \int_{\Omega} f_r L_i \cos\theta_i \, d\omega_i$$
3. 蒙特卡洛方法提供了数值求解高维积分的框架
4. 路径追踪通过递归采样求解渲染方程

**关键公式**：
- 辐射度定义：$L = \frac{d^2\Phi}{dA \cos\theta \, d\omega}$
- BRDF定义：$f_r = \frac{dL_o}{dE_i}$
- 蒙特卡洛估计：$\langle I \rangle = \frac{1}{N} \sum_{i=1}^{N} \frac{f(X_i)}{p(X_i)}$
- MIS权重：$w_i(x) = \frac{p_i^2(x)}{\sum_j p_j^2(x)}$

**算法要点**：
- 重要性采样降低方差
- 俄罗斯轮盘赌保证无偏性
- 直接光照采样提高效率
- 多重重要性采样组合多种策略

## 练习题

### 基础题

**练习8.1**：证明辐射度在真空中沿光线传播时保持不变。
<details>
<summary>提示</summary>
考虑两个微分面元之间的能量传输，使用立体角和面积的关系。
</details>

<details>
<summary>答案</summary>
设两点 $\mathbf{p}_1$ 和 $\mathbf{p}_2$ 之间的辐射传输，从 $\mathbf{p}_1$ 发出的功率：
$$d\Phi = L_1 dA_1 \cos\theta_1 d\omega_1$$
其中 $d\omega_1 = \frac{dA_2 \cos\theta_2}{r^2}$

在 $\mathbf{p}_2$ 接收的功率：
$$d\Phi = L_2 dA_2 \cos\theta_2 d\omega_2$$
其中 $d\omega_2 = \frac{dA_1 \cos\theta_1}{r^2}$

由能量守恒：$L_1 = L_2$
</details>

**练习8.2**：推导Lambertian BRDF $f_r = \rho/\pi$ 满足能量守恒的条件。
<details>
<summary>提示</summary>
计算半球积分 $\int_{\Omega} f_r \cos\theta \, d\omega$。
</details>

<details>
<summary>答案</summary>
$$\int_{\Omega} \frac{\rho}{\pi} \cos\theta \, d\omega = \frac{\rho}{\pi} \int_0^{2\pi} \int_0^{\pi/2} \cos\theta \sin\theta \, d\theta \, d\phi = \frac{\rho}{\pi} \cdot 2\pi \cdot \frac{1}{2} = \rho$$
因此能量守恒要求 $\rho \leq 1$。
</details>

**练习8.3**：计算均匀采样半球方向时，估计器 $\frac{\cos\theta}{p(\omega)}$ 的方差。
<details>
<summary>提示</summary>
均匀采样时 $p(\omega) = 1/(2\pi)$，计算 $E[(\cos\theta)^2]$。
</details>

<details>
<summary>答案</summary>
$$V = E\left[\left(\frac{\cos\theta}{1/(2\pi)}\right)^2\right] - \left(E\left[\frac{\cos\theta}{1/(2\pi)}\right]\right)^2$$
$$= (2\pi)^2 \int_{\Omega} \cos^2\theta \frac{1}{2\pi} d\omega - \pi^2$$
$$= 2\pi \cdot \frac{\pi}{2} - \pi^2 = 0$$
</details>

### 挑战题

**练习8.4**：设计一个采样策略，同时考虑BRDF和入射光照分布，并推导相应的pdf。
<details>
<summary>提示</summary>
考虑乘积 $f_r(\omega) L_i(\omega) \cos\theta$ 的重要性采样。
</details>

<details>
<summary>答案</summary>
理想的pdf应该正比于被积函数：
$$p(\omega) \propto f_r(\omega_i \to \omega_o) L_i(\omega_i) \cos\theta_i$$

实践中可以使用多重重要性采样组合BRDF采样和光源采样：
- BRDF采样：$p_1(\omega) \propto f_r(\omega)$
- 光源采样：$p_2(\omega) \propto L_i(\omega) \cos\theta$

使用MIS权重组合两种策略。
</details>

**练习8.5**：分析路径追踪中不同长度路径的相对贡献，假设场景中材质的平均反射率为 $\bar{\rho}$。
<details>
<summary>提示</summary>
考虑几何级数和路径的概率。
</details>

<details>
<summary>答案</summary>
长度为 $n$ 的路径的平均贡献大约为 $\bar{\rho}^n$。总贡献为几何级数：
$$\sum_{n=0}^{\infty} \bar{\rho}^n = \frac{1}{1-\bar{\rho}}$$

例如，$\bar{\rho} = 0.5$ 时，直接光照贡献50%，一次反射贡献25%，以此类推。这解释了为什么大多数场景中3-5次反射就足够了。
</details>

**练习8.6**：推导双向路径追踪中，长度为 $k$ 的路径有多少种不同的采样策略，并分析每种策略的效率。
<details>
<summary>提示</summary>
考虑从相机追踪 $s$ 步，从光源追踪 $t$ 步，其中 $s + t = k$。
</details>

<details>
<summary>答案</summary>
对于长度为 $k$ 的路径（$k+1$ 个顶点），有 $k+1$ 种采样策略：
- $(s,t) = (0,k)$：纯光源追踪（光子映射）
- $(s,t) = (1,k-1)$：直接光照
- ...
- $(s,t) = (k,0)$：纯相机追踪（路径追踪）

不同策略的效率取决于场景特征：
- 焦散：光源追踪更有效
- 间接漫反射：双向连接更有效
- 镜面反射：相机追踪更有效
</details>

## 常见陷阱与错误

1. **辐射度单位混淆**
   - 错误：将辐照度（W/m²）与辐射度（W/(m²·sr)）混淆
   - 正确：明确区分不同辐射度量及其物理意义

2. **BRDF归一化**
   - 错误：假设 $\int f_r d\omega = 1$
   - 正确：能量守恒要求 $\int f_r \cos\theta_o d\omega_o \leq 1$

3. **概率密度转换**
   - 错误：直接使用立体角域的pdf进行面积采样
   - 正确：$p_A = p_\omega \frac{r^2}{\cos\theta}$

4. **俄罗斯轮盘赌偏差**
   - 错误：不除以继续概率 $P_{rr}$
   - 正确：必须除以 $P_{rr}$ 保证无偏

5. **MIS权重计算**
   - 错误：只考虑当前采样策略的pdf
   - 正确：需要评估所有可能策略的pdf

6. **数值精度问题**
   - 错误：在接近零的pdf处评估
   - 正确：添加小的epsilon避免除零

## 最佳实践检查清单

### 算法实现
- [ ] 验证BRDF满足互易性和能量守恒
- [ ] 正确实现概率密度的测度转换
- [ ] 使用俄罗斯轮盘赌时保证无偏性
- [ ] 实现鲁棒的光线-物体相交，处理数值误差

### 采样策略
- [ ] 针对不同BRDF类型实现专门的采样方法
- [ ] 实现直接光照采样提高收敛速度
- [ ] 使用多重重要性采样组合多种策略
- [ ] 根据材质属性自适应调整俄罗斯轮盘赌概率

### 优化技巧
- [ ] 实现高效的光源采样（如别名方法）
- [ ] 使用空间数据结构加速光线求交
- [ ] 实现路径重用技术（如光子映射）
- [ ] 考虑使用准蒙特卡洛序列

### 验证与调试
- [ ] 实现白炉测试验证能量守恒
- [ ] 使用解析解验证简单场景
- [ ] 实现可视化工具调试采样分布
- [ ] 记录方差和收敛速度统计