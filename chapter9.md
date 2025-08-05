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

**物理直觉**：BRDF描述了一个微分光束如何被表面"重新分配"到各个出射方向。它本质上是一个概率密度函数的光学类比，但带有$\cos\theta_i$项的几何因子。

**单位分析**：BRDF的单位是$\text{sr}^{-1}$（球面度的倒数），这可以从定义式推导：
$$[f_r] = \frac{[\text{W/(m}^2\cdot\text{sr)}]}{[\text{W/(m}^2\cdot\text{sr)}] \cdot 1 \cdot [\text{sr}]} = [\text{sr}^{-1}]$$

**反射方程**：给定入射光场，出射辐射度通过反射方程计算：
$$L_o(\omega_o) = \int_{\Omega} f_r(\omega_i, \omega_o) L_i(\omega_i) \cos\theta_i d\omega_i$$

这是渲染方程的核心组成部分，也是所有基于物理的渲染算法的基础。

### 9.1.2 BRDF的数学性质

BRDF必须满足以下物理约束：

**1. 非负性（Non-negativity）**：
$$f_r(\omega_i, \omega_o) \geq 0$$

物理意义：能量不能为负，表面不能"吸收"后发出负光。

**2. 互易性（Helmholtz Reciprocity）**：
$$f_r(\omega_i, \omega_o) = f_r(\omega_o, \omega_i)$$

这一性质源于光路可逆原理，在实际渲染中可用于优化计算。互易性的深层含义来自于时间反演对称性和洛伦兹互易定理。

**实际应用**：在双向路径追踪中，互易性允许我们交换光线方向而不改变贡献，这对于连接光路至关重要。

**3. 能量守恒（Energy Conservation）**：
$$\int_{\Omega} f_r(\omega_i, \omega_o) \cos\theta_o d\omega_o \leq 1$$

对于任意入射方向，反射的总能量不能超过入射能量。

**方向反照率（Directional Albedo）**：
$$\rho(\omega_i) = \int_{\Omega} f_r(\omega_i, \omega_o) \cos\theta_o d\omega_o$$

这表示从方向$\omega_i$入射的光被反射的总比例。

**半球反照率（Hemispherical Albedo）**：
$$\rho_{hh} = \frac{1}{\pi} \int_{\Omega_i} \int_{\Omega_o} f_r(\omega_i, \omega_o) \cos\theta_i \cos\theta_o d\omega_i d\omega_o$$

**4. 线性性（Linearity）**：
BRDF对入射光强度是线性的，这使得我们可以将复杂光照分解为简单光源的叠加。数学表述：
$$L_o = \int f_r L_i \cos\theta_i d\omega_i = \int f_r (L_{i1} + L_{i2}) \cos\theta_i d\omega_i = L_{o1} + L_{o2}$$

**5. 平滑性与连续性**：
物理BRDF通常是连续的（除了镜面反射的狄拉克δ函数），这保证了渲染结果的视觉连续性。

### 9.1.3 从BRDF到BSDF：透射与散射

双向散射分布函数（Bidirectional Scattering Distribution Function, BSDF）是BRDF的推广，包含了反射和透射：

$$f_s(\omega_i, \omega_o) = f_r(\omega_i, \omega_o) + f_t(\omega_i, \omega_o)$$

其中$f_t$是双向透射分布函数（BTDF）。

**注意**：这里$\omega_o$可能在表面的另一侧（透射情况），需要仔细定义坐标系。

对于透射，需要考虑折射定律（Snell's Law）：
$$\eta_i \sin\theta_i = \eta_t \sin\theta_t$$

其中$\eta_i$和$\eta_t$分别是入射和透射介质的折射率。

**折射方向计算**：
给定入射方向$\omega_i$和表面法线$\mathbf{n}$，折射方向$\omega_t$为：
$$\omega_t = \frac{\eta_i}{\eta_t}\omega_i + \left(\frac{\eta_i}{\eta_t}\cos\theta_i - \sqrt{1 - \left(\frac{\eta_i}{\eta_t}\right)^2(1 - \cos^2\theta_i)}\right)\mathbf{n}$$

**临界角与全内反射**：
当光从光密介质射向光疏介质时，存在临界角：
$$\theta_c = \arcsin\left(\frac{\eta_t}{\eta_i}\right)$$

当$\theta_i > \theta_c$时，发生全内反射，此时BTDF为0，所有能量都被反射。

**BTDF的特殊性质**：
1. **非对称性**：由于折射率不同，BTDF不满足简单的互易性，而是：
   $$\eta_i^2 f_t(\omega_i, \omega_o) = \eta_t^2 f_t(\omega_o, \omega_i)$$

2. **雅可比行列式**：从立体角到立体角的变换需要考虑折射导致的立体角压缩/扩展：
   $$\frac{d\omega_t}{d\omega_i} = \frac{\eta_t^2}{\eta_i^2} \frac{\cos\theta_t}{\cos\theta_i}$$

**理想折射BTDF**：
对于理想光滑表面：
$$f_t(\omega_i, \omega_o) = \frac{\eta_t^2}{\eta_i^2} \frac{|n \cdot \omega_i|}{|n \cdot \omega_t|} (1 - F(\omega_i)) \delta(\omega_o - \omega_t)$$

其中$F(\omega_i)$是菲涅尔反射率，$(1-F)$是透射率。

### 9.1.4 球谐函数与BRDF表示

球谐函数（Spherical Harmonics, SH）提供了在球面上表示函数的正交基：

$$Y_l^m(\theta, \phi) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-|m|)!}{(l+|m|)!}} P_l^{|m|}(\cos\theta) e^{im\phi}$$

其中$P_l^m$是关联勒让德多项式。

**前几阶球谐函数**：
- $Y_0^0 = \frac{1}{2\sqrt{\pi}}$ （常数项）
- $Y_1^{-1} = \frac{1}{2}\sqrt{\frac{3}{\pi}}\sin\theta e^{-i\phi}$
- $Y_1^0 = \frac{1}{2}\sqrt{\frac{3}{\pi}}\cos\theta$
- $Y_1^1 = -\frac{1}{2}\sqrt{\frac{3}{\pi}}\sin\theta e^{i\phi}$

**BRDF的球谐展开**：
由于BRDF是四维函数，需要更复杂的表示。一种方法是固定入射方向：
$$f_r(\omega_i, \omega_o) = \sum_{l=0}^{\infty} \sum_{m=-l}^{l} a_{lm}(\omega_i) Y_l^m(\omega_o)$$

**旋转不变性**：
对于各向同性BRDF，可以利用旋转不变性简化为只依赖于$\omega_i \cdot \omega_o$的函数，使用Legendre多项式展开：
$$f_r(\omega_i \cdot \omega_o) = \sum_{l=0}^{\infty} a_l P_l(\omega_i \cdot \omega_o)$$

**Zonal Harmonics**：
当BRDF具有旋转对称性时，可以使用Zonal Harmonics（$m=0$的球谐函数）：
$$f_r(\theta) = \sum_{l=0}^{\infty} a_l \sqrt{\frac{2l+1}{4\pi}} P_l(\cos\theta)$$

实践中通常截断到较低阶（如$l \leq 4$），这对于漫反射和低频光照效果良好。

**优势与局限**：
- 优势：紧凑表示、快速评估、易于旋转和卷积
- 局限：难以表示高频特征（如锐利高光）、负值问题（需要额外处理）

**应用场景**：
- 预计算辐射传输（Precomputed Radiance Transfer, PRT）
- 环境光照的快速近似
- 实时全局光照中的间接光近似

### 9.1.5 BRDF的测量与拟合

**测量设备**：测角光度计（Gonioreflectometer）通过机械臂控制光源和检测器位置，系统地测量不同角度组合下的反射率。

**测量挑战**：
- **高维度**：4D函数需要大量采样（典型需要$10^6$以上测量点）
- **动态范围**：从漫反射到镜面反射跨越多个数量级
- **掠射角**：接近90°时测量困难但视觉重要
- **时间成本**：完整测量可能需要数小时到数天

**数据表示**：
1. **表格形式**：直接存储测量数据，使用时插值
   - 优点：保真度高
   - 缺点：存储量大、不连续

2. **解析模型拟合**：将测量数据拟合到参数化模型
   - 优点：紧凑、可微、物理意义明确
   - 缺点：可能丢失细节

3. **基函数分解**：使用球谐函数、小波或其他基函数
   - 优点：多分辨率、渐进传输
   - 缺点：高频重建困难

**数据库与标准**：
- MERL BRDF数据库：100种材质的密集测量
- UTIA数据库：包含各向异性材质
- Disney BRDF数据集：生产级材质参考

**拟合方法**：

**1. 最小二乘拟合**：
$$\min_{\theta} \sum_{i} w_i |f_r^{\text{measured}}(\omega_{i,in}, \omega_{i,out}) - f_r^{\text{model}}(\omega_{i,in}, \omega_{i,out}; \theta)|^2$$

**2. 最大似然估计**：
考虑测量噪声模型：
$$\mathcal{L}(\theta) = \prod_i p(y_i | f_r(\omega_{i,in}, \omega_{i,out}; \theta))$$

**3. 贝叶斯方法**：
加入先验知识：
$$p(\theta | \mathcal{D}) \propto p(\mathcal{D} | \theta) p(\theta)$$

**拟合误差度量**：
$$E = \int_{\Omega_i} \int_{\Omega_o} w(\omega_i, \omega_o) |f_r^{\text{measured}} - f_r^{\text{model}}|^2 d\omega_i d\omega_o$$

其中$w$是权重函数，通常强调掠射角等重要区域。

**感知误差度量**：
考虑人眼感知特性：
$$E_{\text{perceptual}} = \int \int w(\omega_i, \omega_o) \Delta E_{CIE}(L_o^{\text{measured}}, L_o^{\text{model}}) d\omega_i d\omega_o$$

其中$\Delta E_{CIE}$是CIE色差公式。

## 9.2 高级材质模型

### 9.2.1 微表面理论（Microfacet Theory）

微表面理论将粗糙表面建模为大量微小镜面的统计分布。核心思想：宏观BRDF是微观几何的统计平均。

**基本假设**：
1. 表面由微小的理想镜面组成
2. 微表面尺度远大于光波长（几何光学有效）
3. 微表面尺度远小于像素（统计平均有意义）
4. 每个微表面是完美镜面（菲涅尔反射）

**微表面BRDF推导**：
从微观到宏观的统计过程：
1. 只有法线为$\omega_h$的微表面能将光从$\omega_i$反射到$\omega_o$
2. 需要考虑微表面的可见性（遮蔽和阴影）
3. 积分所有贡献的微表面

**Cook-Torrance模型**：
$$f_r(\omega_i, \omega_o) = \frac{D(\omega_h) G(\omega_i, \omega_o) F(\omega_i, \omega_h)}{4 \cos\theta_i \cos\theta_o}$$

其中：
- $D$：法线分布函数（Normal Distribution Function）
- $G$：几何遮蔽函数（Geometry Function）
- $F$：菲涅尔项（Fresnel Term）
- $\omega_h$：半程向量（Half Vector），$\omega_h = \frac{\omega_i + \omega_o}{|\omega_i + \omega_o|}$

**分母的物理意义**：
$4 \cos\theta_i \cos\theta_o$项来自于：
- 微表面积分到宏观表面积的雅可比变换
- 将微分立体角转换为微分面积
- 确保能量守恒

**微表面模型的优势**：
1. 物理基础扎实
2. 参数直观（粗糙度、折射率）
3. 能量守恒（正确实现时）
4. 可扩展（各向异性、多次散射等）

### 9.2.2 法线分布函数（D项）

法线分布函数$D(\omega_h)$描述了微表面法线的统计分布，必须满足归一化条件：
$$\int_{\Omega} D(\omega_h) \cos\theta_h d\omega_h = 1$$

这确保了微表面的投影面积等于宏观表面面积。

**物理意义**：
- $D(\omega_h)$表示法线为$\omega_h$的微表面的面积密度
- $D(\omega_h) \cos\theta_h$是投影到宏观表面的面积密度
- 积分为1保证了面积守恒

**Beckmann分布**：
$$D_{\text{Beckmann}}(\omega_h) = \frac{1}{\pi \alpha^2 \cos^4\theta_h} \exp\left(-\frac{\tan^2\theta_h}{\alpha^2}\right)$$

其中$\alpha$是粗糙度参数。

**特性**：
- 基于高斯分布的斜率统计
- 在$\alpha \to 0$时收敛到δ函数（完美镜面）
- 计算成本相对较高（指数函数）

**GGX/Trowbridge-Reitz分布**：
$$D_{\text{GGX}}(\omega_h) = \frac{\alpha^2}{\pi ((\alpha^2 - 1)\cos^2\theta_h + 1)^2}$$

GGX分布具有更长的尾部，能更好地模拟真实材质的高光。

**GGX vs Beckmann**：
- GGX有更长的尾部（重尾分布）
- 更真实的高光形状
- 更好的掠射角表现
- 计算效率更高

**各向异性扩展**：
对于各向异性材质，使用两个粗糙度参数$\alpha_x$和$\alpha_y$：
$$D_{\text{aniso}}(\omega_h) = \frac{1}{\pi \alpha_x \alpha_y} \frac{1}{(\frac{h_x^2}{\alpha_x^2} + \frac{h_y^2}{\alpha_y^2} + h_z^2)^2}$$

其中$(h_x, h_y, h_z)$是半程向量在切空间的坐标。

**其他分布**：
1. **Phong分布**（已过时）：$D = \frac{n+2}{2\pi}\cos^n\theta_h$
2. **Berry分布**：考虑自相关长度的物理模型
3. **ABC分布**：统一框架，可调节分布形状

**粗糙度参数映射**：
实践中常用感知线性的粗糙度：
$$\alpha = \text{roughness}^2$$

这使得参数调节更直观。

### 9.2.3 几何遮蔽函数（G项）

几何函数描述了微表面间的相互遮蔽和阴影效应。

**物理意义**：
- **遮蔽（Masking）**：入射光被其他微表面挡住
- **阴影（Shadowing）**：出射光被其他微表面挡住
- $G \in [0,1]$表示未被遮挡的比例

**Smith模型**：
$$G(\omega_i, \omega_o) = G_1(\omega_i) G_1(\omega_o)$$

其中$G_1$是单向遮蔽函数。

**假设**：
- 遮蔽和阴影统计独立（近似）
- 微表面高度与法线方向无关
- 适用于各向同性粗糙表面

**GGX的Smith-G1**：
$$G_1(\omega) = \frac{2\cos\theta}{1 + \sqrt{1 + \alpha^2 \tan^2\theta}}$$

**推导思路**：
1. 假设微表面高度服从正态分布
2. 使用射线与随机粗糙表面相交的统计理论
3. 得出可见性的解析表达式

**高度相关遮蔽（Height-Correlated Masking-Shadowing）**：
考虑入射和出射方向的相关性：
$$G(\omega_i, \omega_o) = \frac{1}{1 + \Lambda(\omega_i) + \Lambda(\omega_o)}$$

其中$\Lambda$是辅助函数：
$$\Lambda(\omega) = \frac{-1 + \sqrt{1 + \alpha^2 \tan^2\theta}}{2}$$

**相关性的重要性**：
- 独立假设（分离式）在某些情况下过于保守
- 高度相关模型更准确，特别是在掠射角
- 计算成本相近

**其他G函数**：
1. **Cook-Torrance G**：$G = \min(1, \frac{2n \cdot h \cdot n \cdot v}{v \cdot h}, \frac{2n \cdot h \cdot n \cdot l}{v \cdot h})$
2. **Kelemen G**：$G = \frac{n \cdot l \cdot n \cdot v}{v \cdot h}$（快速近似）
3. **V-cavity模型**：基于V形凹槽的简化模型

**实用近似**：
Schlick-GGX近似（用于实时渲染）：
$$G_{\text{Schlick}}(\omega) = \frac{\cos\theta}{\cos\theta(1-k) + k}$$

其中$k = \frac{(\text{roughness}+1)^2}{8}$（直接光照）或$k = \frac{\text{roughness}^2}{2}$（IBL）。

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
