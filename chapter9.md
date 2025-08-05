# 第9章：材质与外观

材质与外观建模是计算机图形学中连接物理世界与虚拟渲染的关键桥梁。本章将深入探讨如何用数学模型描述真实世界中光与物质的相互作用，从基础的BRDF理论到复杂的外观现象。对于AI科学家而言，理解这些模型不仅有助于传统渲染，更是神经渲染、可微渲染等新兴领域的理论基础。

本章学习目标：
- 掌握BRDF/BSDF的数学框架和物理约束
- 理解微表面理论及其在现代渲染中的应用
- 学习复杂外观现象的建模方法
- 了解材质表示在机器学习中的应用前景

## 9.1 BRDF/BSDF理论

### 9.1.1 辐射度量学回顾与BRDF定义

双向反射分布函数（Bidirectional Reflectance Distribution Function, BRDF）描述了光线如何从物体表面反射。回顾辐射度量学基础：

**辐射通量（Radiant Flux）**：$\Phi$，单位时间内的辐射能量，单位：瓦特(W)

**辐照度（Irradiance）**：$E = \frac{d\Phi}{dA}$，单位面积接收的辐射通量，单位：$\text{W/m}^2$

**辐射度（Radiance）**：$L = \frac{d^2\Phi}{dA \cos\theta d\omega}$，单位面积、单位立体角的辐射通量，单位：$\text{W/(m}^2\cdot\text{sr)}$

BRDF的正式定义：

$$f_r(\omega_i, \omega_o) = \frac{dL_o(\omega_o)}{dE_i(\omega_i)} = \frac{dL_o(\omega_o)}{L_i(\omega_i)\cos\theta_i d\omega_i}$$

其中：
- $\omega_i$：入射方向
- $\omega_o$：出射方向
- $L_i$：入射辐射度
- $L_o$：出射辐射度
- $\theta_i$：入射角（与法线夹角）

### 9.1.2 BRDF的数学性质

BRDF必须满足以下物理约束：

**1. 非负性（Non-negativity）**：
$$f_r(\omega_i, \omega_o) \geq 0$$

**2. 互易性（Helmholtz Reciprocity）**：
$$f_r(\omega_i, \omega_o) = f_r(\omega_o, \omega_i)$$

这一性质源于光路可逆原理，在实际渲染中可用于优化计算。

**3. 能量守恒（Energy Conservation）**：
$$\int_{\Omega} f_r(\omega_i, \omega_o) \cos\theta_o d\omega_o \leq 1$$

对于任意入射方向，反射的总能量不能超过入射能量。

**4. 线性性（Linearity）**：
BRDF对入射光强度是线性的，这使得我们可以将复杂光照分解为简单光源的叠加。

### 9.1.3 从BRDF到BSDF：透射与散射

双向散射分布函数（Bidirectional Scattering Distribution Function, BSDF）是BRDF的推广，包含了反射和透射：

$$f_s(\omega_i, \omega_o) = f_r(\omega_i, \omega_o) + f_t(\omega_i, \omega_o)$$

其中$f_t$是双向透射分布函数（BTDF）。

对于透射，需要考虑折射定律（Snell's Law）：
$$\eta_i \sin\theta_i = \eta_t \sin\theta_t$$

其中$\eta_i$和$\eta_t$分别是入射和透射介质的折射率。

**临界角与全内反射**：
当光从光密介质射向光疏介质时，存在临界角：
$$\theta_c = \arcsin\left(\frac{\eta_t}{\eta_i}\right)$$

### 9.1.4 球谐函数与BRDF表示

球谐函数（Spherical Harmonics, SH）提供了在球面上表示函数的正交基：

$$Y_l^m(\theta, \phi) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|}(\cos\theta) e^{im\phi}$$

其中$P_l^m$是关联勒让德多项式。

BRDF的球谐展开：
$$f_r(\omega_i, \omega_o) = \sum_{l=0}^{\infty} \sum_{m=-l}^{l} a_{lm} Y_l^m(\omega_o)$$

实践中通常截断到较低阶（如$l \leq 4$），这对于漫反射和低频光照效果良好。

### 9.1.5 BRDF的测量与拟合

**测量设备**：测角光度计（Gonioreflectometer）通过机械臂控制光源和检测器位置，系统地测量不同角度组合下的反射率。

**数据表示**：
1. **表格形式**：直接存储测量数据，使用时插值
2. **解析模型拟合**：将测量数据拟合到参数化模型
3. **基函数分解**：使用球谐函数、小波或其他基函数

**拟合误差度量**：
$$E = \int_{\Omega_i} \int_{\Omega_o} w(\omega_i, \omega_o) |f_r^{\text{measured}} - f_r^{\text{model}}|^2 d\omega_i d\omega_o$$

其中$w$是权重函数，通常强调掠射角等重要区域。

## 9.2 高级材质模型

### 9.2.1 微表面理论（Microfacet Theory）

微表面理论将粗糙表面建模为大量微小镜面的统计分布。核心思想：宏观BRDF是微观几何的统计平均。

**Cook-Torrance模型**：
$$f_r(\omega_i, \omega_o) = \frac{D(\omega_h) G(\omega_i, \omega_o) F(\omega_i, \omega_h)}{4 \cos\theta_i \cos\theta_o}$$

其中：
- $D$：法线分布函数（Normal Distribution Function）
- $G$：几何遮蔽函数（Geometry Function）
- $F$：菲涅尔项（Fresnel Term）
- $\omega_h$：半程向量（Half Vector），$\omega_h = \frac{\omega_i + \omega_o}{|\omega_i + \omega_o|}$

### 9.2.2 法线分布函数（D项）

法线分布函数$D(\omega_h)$描述了微表面法线的统计分布，必须满足归一化条件：
$$\int_{\Omega} D(\omega_h) \cos\theta_h d\omega_h = 1$$

**Beckmann分布**：
$$D_{\text{Beckmann}}(\omega_h) = \frac{1}{\pi \alpha^2 \cos^4\theta_h} \exp\left(-\frac{\tan^2\theta_h}{\alpha^2}\right)$$

其中$\alpha$是粗糙度参数。

**GGX/Trowbridge-Reitz分布**：
$$D_{\text{GGX}}(\omega_h) = \frac{\alpha^2}{\pi ((\alpha^2 - 1)\cos^2\theta_h + 1)^2}$$

GGX分布具有更长的尾部，能更好地模拟真实材质的高光。

**各向异性扩展**：
对于各向异性材质，使用两个粗糙度参数$\alpha_x$和$\alpha_y$：
$$D_{\text{aniso}}(\omega_h) = \frac{1}{\pi \alpha_x \alpha_y} \frac{1}{(\frac{h_x^2}{\alpha_x^2} + \frac{h_y^2}{\alpha_y^2} + h_z^2)^2}$$

### 9.2.3 几何遮蔽函数（G项）

几何函数描述了微表面间的相互遮蔽和阴影效应。

**Smith模型**：
$$G(\omega_i, \omega_o) = G_1(\omega_i) G_1(\omega_o)$$

其中$G_1$是单向遮蔽函数。

**GGX的Smith-G1**：
$$G_1(\omega) = \frac{2\cos\theta}{1 + \sqrt{1 + \alpha^2 \tan^2\theta}}$$

**高度相关遮蔽（Height-Correlated Masking-Shadowing）**：
考虑入射和出射方向的相关性：
$$G(\omega_i, \omega_o) = \frac{1}{1 + \Lambda(\omega_i) + \Lambda(\omega_o)}$$

其中$\Lambda$是辅助函数：
$$\Lambda(\omega) = \frac{-1 + \sqrt{1 + \alpha^2 \tan^2\theta}}{2}$$

### 9.2.4 菲涅尔反射（F项）与复杂IOR

菲涅尔方程描述了反射率随入射角的变化。

**Schlick近似**：
$$F_{\text{Schlick}}(\omega_i, \omega_h) = F_0 + (1 - F_0)(1 - \cos\theta_i)^5$$

其中$F_0$是垂直入射时的反射率：
$$F_0 = \left(\frac{\eta_1 - \eta_2}{\eta_1 + \eta_2}\right)^2$$

**复数折射率**：
对于导体，折射率是复数$\tilde{\eta} = \eta + i\kappa$，其中$\kappa$是消光系数。

完整的菲涅尔方程（非极化光）：
$$F = \frac{1}{2}(F_s + F_p)$$

其中$F_s$和$F_p$分别是s偏振和p偏振的反射率。

**色散效应**：
折射率随波长变化，可用Cauchy方程或Sellmeier方程建模：
$$\eta(\lambda) = A + \frac{B}{\lambda^2} + \frac{C}{\lambda^4}$$

### 9.2.5 分层材质与BSDF组合

现实材质常由多层组成，如汽车漆（底漆+金属漆+清漆）。

**层间多次反射**：
设上层BSDF为$f_1$，透射率为$T_1$，下层BSDF为$f_2$，则组合BSDF：
$$f_{\text{combined}} = f_1 + \frac{T_1^{\downarrow} f_2 T_1^{\uparrow}}{1 - f_2 * R_1}$$

其中$*$表示卷积，$R_1$是上层的反射算子。

**能量补偿**：
微表面模型可能丢失多次散射能量，需要补偿项：
$$f_r^{\text{ms}} = \frac{(1-E(\omega_i))(1-E(\omega_o))}{\pi(1-E_{\text{avg}})} F_{\text{avg}}$$

其中$E(\omega)$是方向反照率，$E_{\text{avg}}$是平均反照率。

## 9.3 复杂外观建模

### 9.3.1 各向异性材质

各向异性材质的反射特性依赖于表面的切向量方向。

**拉丝金属**：
微表面沿特定方向排列，使用各向异性法线分布：
$$D_{\text{aniso}}(\omega_h) = \frac{1}{\pi \alpha_t \alpha_b} \frac{1}{\left(\frac{(\omega_h \cdot \mathbf{t})^2}{\alpha_t^2} + \frac{(\omega_h \cdot \mathbf{b})^2}{\alpha_b^2} + (\omega_h \cdot \mathbf{n})^2\right)^2}$$

其中$\mathbf{t}$和$\mathbf{b}$是切向量和副切向量。

**织物材质**：
纤维结构导致复杂的散射行为，常用专门的BRDF模型如Ashikhmin-Shirley或微圆柱模型。

**头发渲染**：
Marschner模型将头发散射分解为三个分量：
- R：表面反射
- TT：穿透-穿透
- TRT：穿透-反射-穿透

每个分量有不同的纵向和方位散射函数。

### 9.3.2 次表面散射（皮肤、玉石、牛奶）

次表面散射（Subsurface Scattering, SSS）描述光线进入半透明材质内部，经多次散射后从不同位置射出的现象。

**BSSRDF**（双向次表面散射反射分布函数）：
$$L_o(\mathbf{x}_o, \omega_o) = \int_A \int_{\Omega} S(\mathbf{x}_i, \omega_i, \mathbf{x}_o, \omega_o) L_i(\mathbf{x}_i, \omega_i) \cos\theta_i d\omega_i dA$$

其中$S$是BSSRDF，描述从位置$\mathbf{x}_i$入射的光如何影响位置$\mathbf{x}_o$的出射。

**偶极子近似（Dipole Approximation）**：
$$R_d(r) = \frac{\alpha'}{4\pi} \left[\frac{z_r(1+\sigma_{tr}d_r)e^{-\sigma_{tr}d_r}}{d_r^3} + \frac{z_v(1+\sigma_{tr}d_v)e^{-\sigma_{tr}d_v}}{d_v^3}\right]$$

其中：
- $r$：表面距离
- $\sigma_{tr} = \sqrt{3\sigma_a(\sigma_a + \sigma_s')}$：有效传输系数
- $\alpha' = \sigma_s' / \sigma_{tr}$：约化反照率
- $d_r, d_v$：到实源和虚源的距离

**多极子方法**：
对于薄层材质，需要考虑多次内反射，使用多个镜像源。

**量化散射（Quantized Diffusion）**：
将连续扩散过程离散化，适合GPU实现：
$$u(r,t+\Delta t) = \sum_{i} G(r-r_i, \Delta t) u(r_i, t)$$

### 9.3.3 参与介质与体积渲染

参与介质（烟雾、云、大气）中的光传输由辐射传输方程（RTE）描述：

$$(\omega \cdot \nabla) L(\mathbf{x}, \omega) = -\sigma_t(\mathbf{x}) L(\mathbf{x}, \omega) + \sigma_s(\mathbf{x}) \int_{\Omega} p(\omega', \omega) L(\mathbf{x}, \omega') d\omega' + \sigma_a(\mathbf{x}) L_e(\mathbf{x}, \omega)$$

其中：
- $\sigma_t = \sigma_a + \sigma_s$：消光系数
- $\sigma_a$：吸收系数
- $\sigma_s$：散射系数
- $p(\omega', \omega)$：相函数

**相函数**：
Henyey-Greenstein相函数广泛用于各向异性散射：
$$p_{HG}(\cos\theta) = \frac{1}{4\pi} \frac{1 - g^2}{(1 + g^2 - 2g\cos\theta)^{3/2}}$$

其中$g \in [-1, 1]$是各向异性参数。

**体积渲染积分**：
$$L(\mathbf{x}, \omega) = \int_0^d T(0,t) \sigma_s(\mathbf{x}+t\omega) L_i(\mathbf{x}+t\omega, \omega) dt + T(0,d) L_{\text{background}}$$

其中透射率$T(a,b) = \exp\left(-\int_a^b \sigma_t(\mathbf{x}+s\omega) ds\right)$

### 9.3.4 波动光学效应

当几何特征接近光波长时，需考虑波动光学效应。

**薄膜干涉**：
对于厚度$d$的薄膜，反射率：
$$R(\lambda) = \frac{r_1^2 + r_2^2 + 2r_1r_2\cos(2\beta)}{1 + r_1^2r_2^2 + 2r_1r_2\cos(2\beta)}$$

其中$\beta = 2\pi nd\cos\theta_t / \lambda$，$r_1, r_2$是界面反射系数。

**衍射光栅**：
光栅方程：
$$d(\sin\theta_i + \sin\theta_o) = m\lambda$$

其中$d$是光栅周期，$m$是衍射级次。

**结构色**：
蝴蝶翅膀、CD等的虹彩效应源于微结构的衍射和干涉。可用FDTD（时域有限差分）或严格耦合波分析（RCWA）模拟。

### 9.3.5 数据驱动的外观建模

**BTF（双向纹理函数）**：
$$\text{BTF}(u,v,\omega_i,\omega_o) \in \mathbb{R}^{m \times n \times k^2 \times 3}$$

捕获空间变化的材质属性，需要大量存储和采样。

**神经材质表示**：
使用神经网络表示BRDF：
$$f_r = \text{MLP}(\omega_i, \omega_o; \theta)$$

优势：
- 紧凑表示
- 连续可微
- 易于优化和编辑

**可微渲染与材质反演**：
通过梯度下降优化材质参数：
$$\mathcal{L} = \|I_{\text{rendered}}(\theta) - I_{\text{target}}\|^2$$

使用自动微分框架计算$\frac{\partial \mathcal{L}}{\partial \theta}$。

## 本章小结

本章深入探讨了计算机图形学中材质与外观的数学建模：

**核心概念**：
- BRDF定义：$f_r(\omega_i, \omega_o) = \frac{dL_o(\omega_o)}{L_i(\omega_i)\cos\theta_i d\omega_i}$
- 物理约束：非负性、互易性、能量守恒
- 微表面模型：$f_r = \frac{D \cdot G \cdot F}{4 \cos\theta_i \cos\theta_o}$

**关键公式**：
- GGX法线分布：$D_{GGX} = \frac{\alpha^2}{\pi ((\alpha^2 - 1)\cos^2\theta_h + 1)^2}$
- Smith遮蔽函数：$G = G_1(\omega_i) G_1(\omega_o)$
- Schlick菲涅尔近似：$F = F_0 + (1 - F_0)(1 - \cos\theta)^5$
- 辐射传输方程：$(\omega \cdot \nabla) L = -\sigma_t L + \sigma_s \int p(\omega', \omega) L d\omega' + \sigma_a L_e$

**高级技术**：
- 各向异性材质建模（拉丝金属、织物、头发）
- 次表面散射（偶极子近似、量化扩散）
- 参与介质渲染（相函数、体积积分）
- 波动光学效应（薄膜干涉、衍射光栅）
- 数据驱动方法（BTF、神经材质）

**实践要点**：
- 选择合适的BRDF模型需权衡物理准确性和计算效率
- 能量守恒是保证渲染正确性的关键
- 复杂材质常需组合多种模型
- 可微渲染为材质获取和编辑提供新途径

## 练习题

### 基础题

**习题9.1** 证明BRDF的互易性原理。从麦克斯韦方程组出发，说明为什么$f_r(\omega_i, \omega_o) = f_r(\omega_o, \omega_i)$。

*提示：考虑时间反演对称性和洛伦兹互易定理。*

<details>
<summary>答案</summary>

从麦克斯韦方程组的时间反演对称性出发：
1. 电磁场满足：$\mathbf{E}(-t) = \mathbf{E}(t)$，$\mathbf{B}(-t) = -\mathbf{B}(t)$
2. 坡印廷矢量：$\mathbf{S} = \mathbf{E} \times \mathbf{H}$在时间反演下改变方向
3. 光路可逆：若光线可从A到B，则必可从B到A
4. 在线性、非磁性介质中，散射矩阵满足互易性
5. BRDF作为散射的宏观描述继承此性质

因此$f_r(\omega_i, \omega_o) = f_r(\omega_o, \omega_i)$。
</details>

**习题9.2** 给定粗糙度$\alpha = 0.3$的GGX分布，计算法线与半程向量夹角为30°时的$D$值。

*提示：直接代入GGX公式计算。*

<details>
<summary>答案</summary>

$\theta_h = 30° = \pi/6$，$\cos\theta_h = \sqrt{3}/2$

代入GGX公式：
$$D = \frac{\alpha^2}{\pi ((\alpha^2 - 1)\cos^2\theta_h + 1)^2}$$
$$= \frac{0.09}{\pi ((0.09 - 1) \cdot 0.75 + 1)^2}$$
$$= \frac{0.09}{\pi ((-0.91) \cdot 0.75 + 1)^2}$$
$$= \frac{0.09}{\pi (0.3175)^2}$$
$$≈ 0.284$$
</details>

**习题9.3** 计算空气（$\eta_1 = 1$）到玻璃（$\eta_2 = 1.5$）界面的菲涅尔反射率$F_0$。使用此值，计算入射角为60°时的Schlick近似值。

*提示：先计算$F_0$，再使用Schlick公式。*

<details>
<summary>答案</summary>

垂直入射反射率：
$$F_0 = \left(\frac{\eta_1 - \eta_2}{\eta_1 + \eta_2}\right)^2 = \left(\frac{1 - 1.5}{1 + 1.5}\right)^2 = \left(\frac{-0.5}{2.5}\right)^2 = 0.04$$

入射角60°时，$\cos\theta_i = 0.5$：
$$F = F_0 + (1 - F_0)(1 - \cos\theta_i)^5$$
$$= 0.04 + 0.96 \times 0.5^5$$
$$= 0.04 + 0.96 \times 0.03125$$
$$= 0.04 + 0.03 = 0.07$$
</details>

**习题9.4** 对于Henyey-Greenstein相函数，当$g = 0.8$时，计算前向散射（$\theta = 0°$）和后向散射（$\theta = 180°$）的相函数值之比。

*提示：代入HG相函数公式，注意$\cos(0°) = 1$，$\cos(180°) = -1$。*

<details>
<summary>答案</summary>

HG相函数：$p_{HG}(\cos\theta) = \frac{1}{4\pi} \frac{1 - g^2}{(1 + g^2 - 2g\cos\theta)^{3/2}}$

前向散射（$\cos\theta = 1$）：
$$p_{HG}(1) = \frac{1}{4\pi} \frac{1 - 0.64}{(1 + 0.64 - 1.6)^{3/2}} = \frac{1}{4\pi} \frac{0.36}{0.04^{3/2}} = \frac{1}{4\pi} \frac{0.36}{0.008} = \frac{45}{4\pi}$$

后向散射（$\cos\theta = -1$）：
$$p_{HG}(-1) = \frac{1}{4\pi} \frac{0.36}{(1 + 0.64 + 1.6)^{3/2}} = \frac{1}{4\pi} \frac{0.36}{3.24^{3/2}} = \frac{1}{4\pi} \frac{0.36}{5.832} ≈ \frac{0.0617}{4\pi}$$

比值：$\frac{p_{HG}(1)}{p_{HG}(-1)} = \frac{45}{0.0617} ≈ 729$
</details>

### 挑战题

**习题9.5** 推导各向异性GGX分布的归一化条件。证明：
$$\int_{\Omega} D_{aniso}(\omega_h) \cos\theta_h d\omega_h = 1$$

*提示：使用球坐标系，将$\omega_h$分解为切空间坐标。*

<details>
<summary>答案</summary>

设切空间坐标系，$\omega_h = (h_x, h_y, h_z)$，其中$h_x = \sin\theta\cos\phi$，$h_y = \sin\theta\sin\phi$，$h_z = \cos\theta$。

各向异性GGX：
$$D_{aniso} = \frac{1}{\pi \alpha_x \alpha_y} \frac{1}{(\frac{h_x^2}{\alpha_x^2} + \frac{h_y^2}{\alpha_y^2} + h_z^2)^2}$$

积分：
$$\int_0^{2\pi} \int_0^{\pi/2} D_{aniso} \cos\theta \sin\theta d\theta d\phi$$

令$u = \sin^2\theta$，则：
$$= \frac{1}{\pi \alpha_x \alpha_y} \int_0^{2\pi} \int_0^1 \frac{1-u}{2(a(\phi)u + (1-u))^2} du d\phi$$

其中$a(\phi) = \frac{\cos^2\phi}{\alpha_x^2} + \frac{\sin^2\phi}{\alpha_y^2} - 1$。

经过复杂积分（使用留数定理或数值验证），结果为1。
</details>

**习题9.6** 设计一个BRDF模型，同时具有漫反射、镜面反射和逆反射（retroreflection）特性。给出数学表达式并证明满足能量守恒。

*提示：考虑将多个BRDF分量线性组合，逆反射可用窄高斯分布在$\omega_o = -\omega_i$附近建模。*

<details>
<summary>答案</summary>

组合BRDF：
$$f_r = k_d f_d + k_s f_s + k_r f_r$$

其中：
- 漫反射：$f_d = \frac{1}{\pi}$
- 镜面反射：$f_s = \frac{D \cdot G \cdot F}{4 \cos\theta_i \cos\theta_o}$（微表面模型）
- 逆反射：$f_r = \frac{A}{\cos\theta_i} \exp\left(-\frac{|\omega_o + \omega_i|^2}{2\sigma^2}\right)$

能量守恒要求：
1. $k_d + k_s + k_r \leq 1$
2. 对逆反射项，需要：$A = \frac{1}{2\pi\sigma^2(1-e^{-2/\sigma^2})}$

验证：
$$\int_{\Omega} f_r \cos\theta_o d\omega_o = k_d + k_s + k_r \int_{\Omega} \frac{A}{\cos\theta_i} \exp(...) \cos\theta_o d\omega_o$$

通过适当选择$\sigma$和归一化常数$A$，可保证总反射率不超过1。
</details>

**习题9.7** 推导薄膜干涉的反射光谱，考虑多次内反射。解释为什么肥皂泡呈现彩虹色。

*提示：使用光程差和复振幅叠加，考虑相位关系。*

<details>
<summary>答案</summary>

设薄膜厚度$d$，折射率$n$，入射角$\theta_i$，折射角$\theta_t$。

光程差：$\Delta = 2nd\cos\theta_t$

考虑多次反射，总反射振幅：
$$r_{total} = r_1 + \frac{t_1 t_1' r_2 e^{i\delta}}{1 - r_1' r_2 e^{i\delta}}$$

其中$\delta = 4\pi nd\cos\theta_t / \lambda$。

反射率：
$$R = |r_{total}|^2 = \frac{r_1^2 + r_2^2 + 2r_1r_2\cos\delta}{1 + r_1^2r_2^2 + 2r_1r_2\cos\delta}$$

彩虹色原因：
1. 不同波长有不同的相位差$\delta$
2. 相长干涉条件：$\Delta = m\lambda$（亮纹）
3. 相消干涉条件：$\Delta = (m+1/2)\lambda$（暗纹）
4. 薄膜厚度变化导致不同位置强化不同颜色
5. 观察角度变化也改变光程差，产生颜色变化
</details>

**习题9.8** 分析神经网络表示BRDF的优缺点。设计一个保证物理约束（互易性、能量守恒）的网络架构。

*提示：考虑如何在网络结构或损失函数中编码物理约束。*

<details>
<summary>答案</summary>

神经BRDF优缺点：

优点：
- 紧凑表示复杂BRDF
- 连续可微，适合优化
- 可学习未知材质
- 易于插值和编辑

缺点：
- 难以保证物理约束
- 需要大量训练数据
- 推理计算成本
- 缺乏可解释性

物理约束架构设计：

1. **互易性**：使用对称输入编码
   ```
   input = sort([ω_i, ω_o]) 或 input = [ω_i + ω_o, |ω_i - ω_o|]
   ```

2. **能量守恒**：输出层归一化
   ```
   f_r_raw = MLP(input)
   albedo = ∫f_r_raw cosθ dω (预计算或学习)
   f_r = f_r_raw / max(1, albedo)
   ```

3. **非负性**：使用ReLU或Softplus激活

4. **损失函数**：
   ```
   L = L_data + λ_1 L_reciprocity + λ_2 L_energy
   L_reciprocity = |f_r(ω_i,ω_o) - f_r(ω_o,ω_i)|²
   L_energy = max(0, ∫f_r cosθ dω - 1)²
   ```

5. **架构示例**：
   - 输入：6D（两个方向）→ 4D（对称编码）
   - 隐藏层：[128, 256, 256, 128]，使用SIREN激活
   - 输出：RGB值，经过物理约束处理
</details>

## 常见陷阱与错误 (Gotchas)

### 1. BRDF实现错误

**错误**：忘记半程向量归一化
```
half = wi + wo  // 错误！
```
**正确**：
```
half = normalize(wi + wo)
```

**错误**：混淆角度定义
- $\theta_h$是半程向量与法线夹角
- $\theta_i$是入射方向与法线夹角
- $\theta_d$是入射方向与半程向量夹角（用于菲涅尔）

### 2. 能量守恒违反

**陷阱**：微表面BRDF可能超过1
- 原因：掠射角时G项不足
- 解决：使用改进的G函数或能量补偿

**陷阱**：多层材质能量累加错误
- 错误：直接相加各层BRDF
- 正确：考虑层间多次反射

### 3. 数值稳定性问题

**GGX在$\alpha = 0$时的奇异性**：
```
// 错误：直接使用alpha
D = alpha^2 / (pi * denominator^2)

// 正确：钳制最小值
alpha = max(0.001, roughness^2)
```

**菲涅尔计算的边界情况**：
- 全内反射时需特殊处理
- 复数IOR的平方根需选择正确分支

### 4. 采样效率问题

**错误**：使用均匀采样计算BRDF积分
- 问题：收敛极慢，噪声大
- 解决：使用重要性采样

**错误**：忽略多重重要性采样（MIS）
- 问题：某些情况下方差很大
- 解决：组合BRDF采样和光源采样

### 5. 物理准确性错误

**混淆辐射度量单位**：
- Radiance: W/(m²·sr)
- Irradiance: W/m²
- BRDF: 1/sr

**忽略色散**：
- 钻石、棱镜等材质必须考虑色散
- 至少对RGB三通道使用不同IOR

### 6. 性能优化陷阱

**预计算错误**：
- LUT采样不足导致失真
- 插值方法选择不当

**着色器精度问题**：
- 半精度浮点可能导致明显误差
- 特别是在计算very rough材质时

## 最佳实践检查清单

### 设计阶段
- [ ] 明确材质的物理特性（金属/电介质、各向同性/异性）
- [ ] 确定所需的外观特征（光泽度、颜色、纹理）
- [ ] 评估性能需求vs质量需求
- [ ] 考虑是否需要艺术可控性

### 实现阶段
- [ ] 验证BRDF满足互易性
- [ ] 检查能量守恒（白炉测试）
- [ ] 实现重要性采样
- [ ] 处理数值稳定性边界情况
- [ ] 添加能量补偿项（如需要）

### 验证阶段
- [ ] 与参考实现对比（如Mitsuba、PBRT）
- [ ] 测试极限参数（roughness=0/1）
- [ ] 验证不同光照条件
- [ ] 检查时间相关性（动画中的闪烁）

### 优化阶段
- [ ] 分析性能瓶颈
- [ ] 考虑近似方法（如Schlick vs完整菲涅尔）
- [ ] 实现LOD策略
- [ ] 预计算可重用项

### 艺术家工作流
- [ ] 提供直观的参数（避免物理单位）
- [ ] 实现实时预览
- [ ] 提供预设材质库
- [ ] 支持纹理映射所有参数

### 扩展性考虑
- [ ] 模块化设计（易于添加新BRDF）
- [ ] 支持材质分层
- [ ] 考虑未来功能（如光谱渲染）
- [ ] 文档化所有假设和限制
