# 第1章：计算机图形学概述与数学基础

计算机图形学是研究如何通过计算生成、操作和显示视觉内容的学科。本章将从宏观视角介绍图形学的核心概念，并建立必要的数学基础。我们将深入探讨向量空间、线性变换以及四元数等高级数学工具，这些内容构成了后续章节的理论基石。通过本章学习，读者将掌握图形学中的基本数学语言，并理解各种数学工具在实际应用中的优势与局限。

## 1.1 计算机图形学概述

### 1.1.1 图形学的本质与应用领域

计算机图形学的核心任务是将三维世界的几何信息转换为二维图像。这个过程涉及几何建模、光照计算、投影变换和图像合成等多个环节。在现代应用中，图形学已经渗透到游戏开发、电影特效、科学可视化、虚拟现实、计算机辅助设计等众多领域。

从数学角度看，图形学本质上是一系列函数映射：
$$f: \mathbb{R}^3 \times \mathcal{L} \times \mathcal{M} \rightarrow \mathbb{R}^2 \times \mathcal{C}$$

其中$\mathbb{R}^3$表示三维空间坐标，$\mathcal{L}$表示光照条件，$\mathcal{M}$表示材质属性，$\mathbb{R}^2$表示图像平面坐标，$\mathcal{C}$表示颜色空间。

### 1.1.2 实时渲染与离线渲染

图形学应用可分为两大类：

**实时渲染**：追求在有限时间内（通常16-33毫秒）生成图像，常用于游戏和交互应用。主要技术包括：
- 光栅化管线（Z-buffer、着色器）
- 层次细节（LOD）优化
- 时间相关性利用（时序抗锯齿、运动模糊）

**离线渲染**：追求物理准确性和视觉真实感，计算时间可达数小时。核心技术包括：
- 光线追踪与路径追踪
- 全局光照算法
- 基于物理的渲染（PBR）

### 1.1.3 图形管线概述

现代图形管线可以抽象为以下阶段：

1. **应用阶段**：场景管理、碰撞检测、动画更新
2. **几何处理**：顶点变换、裁剪、图元装配
3. **光栅化**：将几何图元转换为像素片段
4. **片段处理**：着色计算、深度测试、混合

每个阶段都涉及特定的数学变换和算法优化，我们将在后续章节详细探讨。

### 1.1.4 现代图形学的挑战与机遇

随着AI技术的发展，图形学面临新的机遇：

- **神经渲染**：使用深度学习加速或替代传统渲染算法
- **可微渲染**：将渲染过程设计为可微函数，支持基于梯度的优化
- **实时光线追踪**：硬件加速使得实时光线追踪成为可能
- **虚拟现实与增强现实**：对延迟和分辨率提出更高要求

## 1.2 向量与线性代数

### 1.2.1 向量空间与基本运算

向量是图形学中最基本的数学对象。在$n$维欧几里得空间$\mathbb{R}^n$中，向量$\mathbf{v} = (v_1, v_2, ..., v_n)^T$表示一个有向线段。

**向量加法与标量乘法**：
$$\mathbf{u} + \mathbf{v} = (u_1 + v_1, u_2 + v_2, ..., u_n + v_n)^T$$
$$k\mathbf{v} = (kv_1, kv_2, ..., kv_n)^T$$

**内积（点积）**：
$$\mathbf{u} \cdot \mathbf{v} = \sum_{i=1}^{n} u_i v_i = |\mathbf{u}||\mathbf{v}|\cos\theta$$

内积的几何意义：
- 计算投影长度：$\text{proj}_{\mathbf{v}}\mathbf{u} = \frac{\mathbf{u} \cdot \mathbf{v}}{|\mathbf{v}|^2}\mathbf{v}$
- 判断方向关系：$\mathbf{u} \cdot \mathbf{v} > 0$表示夹角小于90°
- 计算功和能量：在物理模拟中广泛应用

**外积（叉积）**：仅在三维空间定义
$$\mathbf{u} \times \mathbf{v} = \begin{pmatrix}
u_2v_3 - u_3v_2 \\
u_3v_1 - u_1v_3 \\
u_1v_2 - u_2v_1
\end{pmatrix}$$

外积的性质：
- $|\mathbf{u} \times \mathbf{v}| = |\mathbf{u}||\mathbf{v}|\sin\theta$
- 方向遵循右手定则
- 反交换律：$\mathbf{u} \times \mathbf{v} = -\mathbf{v} \times \mathbf{u}$

### 1.2.2 矩阵与线性变换

线性变换$T: \mathbb{R}^n \rightarrow \mathbb{R}^m$可以用$m \times n$矩阵$\mathbf{A}$表示：
$$T(\mathbf{v}) = \mathbf{A}\mathbf{v}$$

**基本变换矩阵**：

缩放矩阵：
$$\mathbf{S} = \begin{pmatrix}
s_x & 0 & 0 \\
0 & s_y & 0 \\
0 & 0 & s_z
\end{pmatrix}$$

绕$z$轴旋转$\theta$角：
$$\mathbf{R}_z(\theta) = \begin{pmatrix}
\cos\theta & -\sin\theta & 0 \\
\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{pmatrix}$$

**矩阵运算的几何意义**：
- 矩阵乘法对应变换的复合
- 逆矩阵对应逆变换
- 行列式表示体积缩放因子
- 特征值和特征向量揭示变换的不变方向

### 1.2.3 齐次坐标与仿射变换

为了统一表示线性变换和平移，引入齐次坐标。三维点$(x, y, z)$的齐次坐标为$(x, y, z, 1)$，向量的齐次坐标为$(x, y, z, 0)$。

仿射变换的一般形式：
$$\mathbf{T} = \begin{pmatrix}
\mathbf{A}_{3×3} & \mathbf{t}_{3×1} \\
\mathbf{0}_{1×3} & 1
\end{pmatrix}$$

其中$\mathbf{A}$是线性变换部分，$\mathbf{t}$是平移向量。

**变换的组合与分解**：
- 任意刚体变换可分解为旋转和平移
- 任意仿射变换可通过SVD分解为旋转、缩放和剪切
- 变换顺序很重要：$\mathbf{T}_1\mathbf{T}_2 \neq \mathbf{T}_2\mathbf{T}_1$

### 1.2.4 正交性与投影

**正交矩阵**满足$\mathbf{Q}^T\mathbf{Q} = \mathbf{I}$，具有以下性质：
- 保持长度和角度
- 行列式为±1
- 逆矩阵等于转置矩阵

**投影矩阵**：
正交投影到子空间$\text{span}(\mathbf{v}_1, ..., \mathbf{v}_k)$：
$$\mathbf{P} = \mathbf{V}(\mathbf{V}^T\mathbf{V})^{-1}\mathbf{V}^T$$

其中$\mathbf{V} = [\mathbf{v}_1 | ... | \mathbf{v}_k]$。

透视投影矩阵（后续章节详述）将三维点映射到规范化设备坐标。

## 1.3 四元数与复数在图形学中的应用

### 1.3.1 复数与二维旋转

复数提供了优雅的二维旋转表示。复数$z = a + bi$可以写成极坐标形式：
$$z = r(\cos\theta + i\sin\theta) = re^{i\theta}$$

二维点$(x, y)$对应复数$x + yi$，绕原点旋转$\phi$角等价于乘以$e^{i\phi}$：
$$(x + yi)e^{i\phi} = (x + yi)(\cos\phi + i\sin\phi)$$
$$= (x\cos\phi - y\sin\phi) + i(x\sin\phi + y\cos\phi)$$

这恰好对应旋转矩阵的作用结果。复数乘法的几何意义：
- 模相乘：$|z_1z_2| = |z_1||z_2|$
- 幅角相加：$\arg(z_1z_2) = \arg(z_1) + \arg(z_2)$

### 1.3.2 四元数基础

四元数是复数的推广，形式为：
$$q = w + xi + yj + zk$$

其中$i^2 = j^2 = k^2 = ijk = -1$，且：
- $ij = k, jk = i, ki = j$
- $ji = -k, kj = -i, ik = -j$

四元数可以写成标量-向量形式：
$$q = (s, \mathbf{v}) = s + v_xi + v_yj + v_zk$$

**四元数运算**：

加法：$(s_1, \mathbf{v}_1) + (s_2, \mathbf{v}_2) = (s_1 + s_2, \mathbf{v}_1 + \mathbf{v}_2)$

乘法（Hamilton积）：
$$q_1q_2 = (s_1s_2 - \mathbf{v}_1 \cdot \mathbf{v}_2, s_1\mathbf{v}_2 + s_2\mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2)$$

共轭：$q^* = (s, -\mathbf{v})$

模：$|q| = \sqrt{qq^*} = \sqrt{s^2 + |\mathbf{v}|^2}$

### 1.3.3 四元数与三维旋转

单位四元数可以表示三维旋转。绕单位轴$\mathbf{n}$旋转$\theta$角对应的四元数为：
$$q = \cos\frac{\theta}{2} + \sin\frac{\theta}{2}(n_xi + n_yj + n_zk)$$

**旋转公式**：点$\mathbf{p}$经四元数$q$旋转后得到：
$$\mathbf{p}' = q\mathbf{p}q^*$$

其中$\mathbf{p}$被视为纯四元数$(0, \mathbf{p})$。

**四元数转旋转矩阵**：
设$q = w + xi + yj + zk$，对应的旋转矩阵为：
$$\mathbf{R} = \begin{pmatrix}
1-2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\
2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\
2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2)
\end{pmatrix}$$

### 1.3.4 四元数的优势与应用

**相比欧拉角的优势**：
1. 无万向锁问题
2. 插值更平滑（球面线性插值SLERP）
3. 数值稳定性更好
4. 存储更紧凑（4个数vs9个数）

**球面线性插值（SLERP）**：
$$\text{slerp}(q_1, q_2, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}q_1 + \frac{\sin(t\Omega)}{\sin\Omega}q_2$$

其中$\cos\Omega = q_1 \cdot q_2$。

**在图形学中的应用**：
- 骨骼动画的关节旋转
- 相机控制（避免万向锁）
- 物理模拟中的刚体方向
- 法线贴图的切线空间变换

### 1.3.5 对偶四元数与刚体变换

对偶四元数将旋转和平移统一表示：
$$\hat{q} = q_r + \epsilon q_d$$

其中$\epsilon^2 = 0$（对偶单位），$q_r$表示旋转，$q_d = \frac{1}{2}tq_r$（$t$是平移四元数）。

对偶四元数的优势：
- 统一表示刚体变换
- 支持平滑插值（包括旋转和平移）
- 避免矩阵乘法的数值误差累积

## 本章小结

本章建立了计算机图形学的数学基础：

**核心概念**：
1. 图形学管线：应用→几何→光栅化→片段处理
2. 向量运算：点积（投影、角度）、叉积（法向、面积）
3. 矩阵变换：线性变换、仿射变换、齐次坐标
4. 四元数：紧凑的旋转表示，避免万向锁

**关键公式**：
- 点积：$\mathbf{u} \cdot \mathbf{v} = |\mathbf{u}||\mathbf{v}|\cos\theta$
- 叉积：$|\mathbf{u} \times \mathbf{v}| = |\mathbf{u}||\mathbf{v}|\sin\theta$
- 四元数旋转：$\mathbf{p}' = q\mathbf{p}q^*$
- SLERP：$\text{slerp}(q_1, q_2, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}q_1 + \frac{\sin(t\Omega)}{\sin\Omega}q_2$

**重要性质**：
- 正交矩阵保持长度和角度
- 单位四元数形成3-球面$S^3$
- 旋转的复合对应四元数乘法
- 仿射变换保持直线和平行关系

## 练习题

### 基础题

**1.1 向量运算基础**
给定向量$\mathbf{a} = (1, 2, 3)$，$\mathbf{b} = (4, -1, 2)$，计算：
- (a) $\mathbf{a} \cdot \mathbf{b}$
- (b) $\mathbf{a} \times \mathbf{b}$
- (c) $\mathbf{a}$在$\mathbf{b}$上的投影向量

*Hint*: 投影公式为$\text{proj}_{\mathbf{b}}\mathbf{a} = \frac{\mathbf{a} \cdot \mathbf{b}}{|\mathbf{b}|^2}\mathbf{b}$

<details>
<summary>答案</summary>

(a) $\mathbf{a} \cdot \mathbf{b} = 1×4 + 2×(-1) + 3×2 = 8$

(b) $\mathbf{a} \times \mathbf{b} = \begin{pmatrix} 2×2-3×(-1) \\ 3×4-1×2 \\ 1×(-1)-2×4 \end{pmatrix} = \begin{pmatrix} 7 \\ 10 \\ -9 \end{pmatrix}$

(c) $|\mathbf{b}|^2 = 16 + 1 + 4 = 21$，所以投影向量为$\frac{8}{21}(4, -1, 2) = (\frac{32}{21}, -\frac{8}{21}, \frac{16}{21})$

</details>

**1.2 矩阵变换组合**
一个物体先绕$z$轴旋转45°，再沿$x$轴方向缩放2倍，最后平移$(3, 1, 0)$。写出复合变换矩阵。

*Hint*: 注意变换顺序，使用齐次坐标

<details>
<summary>答案</summary>

变换顺序（从右到左）：旋转→缩放→平移

$$\mathbf{T} = \begin{pmatrix}
1 & 0 & 0 & 3 \\
0 & 1 & 0 & 1 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
2 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
\cos45° & -\sin45° & 0 & 0 \\
\sin45° & \cos45° & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}$$

$$= \begin{pmatrix}
\sqrt{2} & -\frac{\sqrt{2}}{2} & 0 & 3 \\
\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2} & 0 & 1 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}$$

</details>

**1.3 四元数旋转**
将四元数$q = \frac{1}{2} + \frac{1}{2}i + \frac{1}{2}j + \frac{1}{2}k$转换为旋转矩阵。

*Hint*: 首先验证这是单位四元数

<details>
<summary>答案</summary>

验证：$|q| = \sqrt{(\frac{1}{2})^2 + (\frac{1}{2})^2 + (\frac{1}{2})^2 + (\frac{1}{2})^2} = 1$ ✓

设$w = \frac{1}{2}, x = y = z = \frac{1}{2}$，代入公式：

$$\mathbf{R} = \begin{pmatrix}
1-2(\frac{1}{4}+\frac{1}{4}) & 2(\frac{1}{4}-\frac{1}{4}) & 2(\frac{1}{4}+\frac{1}{4}) \\
2(\frac{1}{4}+\frac{1}{4}) & 1-2(\frac{1}{4}+\frac{1}{4}) & 2(\frac{1}{4}-\frac{1}{4}) \\
2(\frac{1}{4}-\frac{1}{4}) & 2(\frac{1}{4}+\frac{1}{4}) & 1-2(\frac{1}{4}+\frac{1}{4})
\end{pmatrix}$$

$$= \begin{pmatrix}
0 & 0 & 1 \\
1 & 0 & 0 \\
0 & 1 & 0
\end{pmatrix}$$

这是一个120°旋转，轴为$(1,1,1)$方向。

</details>

### 挑战题

**1.4 最小二乘投影**
推导将向量$\mathbf{b}$投影到由列向量$\mathbf{a}_1, \mathbf{a}_2$张成的平面上的投影矩阵$\mathbf{P}$，并证明$\mathbf{P}^2 = \mathbf{P}$。

*Hint*: 使用正规方程$\mathbf{A}^T\mathbf{A}\mathbf{x} = \mathbf{A}^T\mathbf{b}$

<details>
<summary>答案</summary>

设$\mathbf{A} = [\mathbf{a}_1 | \mathbf{a}_2]$，投影为$\mathbf{p} = \mathbf{A}\mathbf{x}$。

最小化$|\mathbf{b} - \mathbf{A}\mathbf{x}|^2$，得正规方程：
$$\mathbf{A}^T\mathbf{A}\mathbf{x} = \mathbf{A}^T\mathbf{b}$$

解得：$\mathbf{x} = (\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\mathbf{b}$

投影：$\mathbf{p} = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\mathbf{b}$

因此投影矩阵：$\mathbf{P} = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T$

证明幂等性：
$$\mathbf{P}^2 = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T = \mathbf{P}$$

</details>

**1.5 四元数插值优化**
给定两个单位四元数$q_1$和$q_2$，当$q_1 \cdot q_2 < 0$时，SLERP会选择"长路径"。如何修正这个问题？证明你的方法保持插值的连续性。

*Hint*: 考虑四元数的双覆盖性质

<details>
<summary>答案</summary>

四元数$q$和$-q$表示相同的旋转（双覆盖）。当$q_1 \cdot q_2 < 0$时，$q_1$和$q_2$在4D球面上的夹角大于90°。

修正方法：若$q_1 \cdot q_2 < 0$，将$q_2$替换为$-q_2$。

证明连续性：
设原始插值路径为$\gamma(t) = \text{slerp}(q_1, q_2, t)$，修正后为$\gamma'(t) = \text{slerp}(q_1, -q_2, t)$。

在$t=0$：$\gamma'(0) = q_1 = \gamma(0)$
在$t=1$：$\gamma'(1) = -q_2$，但$-q_2$和$q_2$表示相同旋转

修正后的夹角$\Omega' = \arccos(-q_1 \cdot q_2) = \pi - \Omega < \pi/2$，确保选择短路径。

</details>

**1.6 仿射变换的不动点**
证明：任意二维仿射变换$T(\mathbf{x}) = \mathbf{A}\mathbf{x} + \mathbf{b}$（其中$\det(\mathbf{A}-\mathbf{I}) \neq 0$）有唯一不动点。找出该不动点的表达式。

*Hint*: 不动点满足$T(\mathbf{x}) = \mathbf{x}$

<details>
<summary>答案</summary>

设不动点为$\mathbf{x}_0$，则：
$$\mathbf{A}\mathbf{x}_0 + \mathbf{b} = \mathbf{x}_0$$
$$(\mathbf{A} - \mathbf{I})\mathbf{x}_0 = -\mathbf{b}$$

由于$\det(\mathbf{A}-\mathbf{I}) \neq 0$，矩阵$\mathbf{A}-\mathbf{I}$可逆，因此：
$$\mathbf{x}_0 = -(\mathbf{A} - \mathbf{I})^{-1}\mathbf{b}$$

这是唯一解。

几何意义：不动点是变换的"中心"，所有点围绕它进行变换。

</details>

**1.7 复数与四元数的联系**
将二维旋转的复数表示推广到三维。具体地，如何用两个复数表示一个三维旋转？这与四元数有什么关系？

*Hint*: 考虑Cayley-Dickson构造

<details>
<summary>答案</summary>

四元数可以看作"复数对"：$q = z_1 + z_2j$，其中$z_1 = a + bi$，$z_2 = c + di$。

这给出：$q = a + bi + cj + dk$

关系：
- 当$c = d = 0$时，退化为复数（二维旋转）
- 四元数乘法规则可从复数乘法推广得出

另一种表示：使用两个复数$(z_1, z_2)$表示三维旋转，满足$|z_1|^2 + |z_2|^2 = 1$。这实际上是$SU(2)$群的表示，与单位四元数同构。

对应关系：
$$\begin{pmatrix} z_1 & -\bar{z}_2 \\ z_2 & \bar{z}_1 \end{pmatrix} \leftrightarrow a + bi + cj + dk$$

其中$z_1 = a + bi$，$z_2 = c + di$。

</details>

**1.8 图形管线的能量守恒**
在物理真实的渲染中，为什么需要确保BRDF满足能量守恒？给出数学表述并解释其在实时渲染中的简化方法。

*Hint*: 考虑反射率的积分

<details>
<summary>答案</summary>

能量守恒要求出射能量不超过入射能量。对于BRDF $f_r(\omega_i, \omega_o)$：

$$\int_{\Omega} f_r(\omega_i, \omega_o) \cos\theta_i \, d\omega_i \leq 1$$

对所有出射方向$\omega_o$成立。

物理意义：反射的总能量不能超过入射能量。

实时渲染简化：
1. 预计算积分表（如UE4的环境BRDF查找表）
2. 使用归一化的解析BRDF（如归一化Phong）
3. 分离漫反射和镜面反射，分别确保守恒

例如，对于简化的BRDF：
$$f_r = \frac{k_d}{\pi} + k_s \cdot f_{specular}$$

确保$k_d + k_s \leq 1$（菲涅尔项会进一步调整）。

</details>

## 常见陷阱与错误

### 向量与点的混淆
- **错误**：将向量当作点进行平移变换
- **正确**：向量的齐次坐标第四分量为0，不受平移影响
- **调试技巧**：检查齐次坐标的w分量

### 矩阵乘法顺序
- **错误**：`translate * rotate * scale`（错误顺序）
- **正确**：`translate * rotate * scale * vertex`（从右向左）
- **记忆方法**：变换从右向左应用，最右边的先执行

### 四元数归一化
- **错误**：忘记归一化导致缩放
- **正确**：每次运算后重新归一化
- **优化**：使用平方根的快速近似算法

### 叉积的手性
- **错误**：混淆左手系和右手系
- **正确**：OpenGL用右手系，DirectX传统上用左手系
- **验证**：$\mathbf{x} \times \mathbf{y} = \mathbf{z}$（右手系）

### 浮点精度问题
- **错误**：直接比较浮点数相等
- **正确**：使用epsilon容差：`abs(a - b) < epsilon`
- **选择epsilon**：通常用`1e-6`，但需根据场景调整

### Gimbal Lock（万向锁）
- **症状**：欧拉角在某些角度失去一个自由度
- **解决**：使用四元数或轴角表示
- **转换时机**：仅在最终输出时转换为欧拉角

## 最佳实践检查清单

### 数学库设计
- [ ] 向量和点使用不同的类型或明确区分
- [ ] 提供行主序和列主序的明确说明
- [ ] 实现SIMD优化的向量运算
- [ ] 包含数值稳定的归一化函数
- [ ] 提供便捷的调试输出功能

### 变换系统
- [ ] 使用齐次坐标统一处理变换
- [ ] 缓存常用的变换矩阵
- [ ] 提供矩阵分解功能（提取旋转、缩放等）
- [ ] 实现视锥体裁剪的快速路径
- [ ] 支持不同坐标系之间的转换

### 四元数使用
- [ ] 始终保持四元数归一化
- [ ] 实现SLERP和NLERP两种插值
- [ ] 处理四元数点积为负的情况
- [ ] 提供与欧拉角、轴角的转换函数
- [ ] 考虑使用对偶四元数处理刚体变换

### 性能优化
- [ ] 避免不必要的矩阵求逆
- [ ] 使用快速平方根倒数近似
- [ ] 预计算能复用的中间结果
- [ ] 考虑内存对齐和缓存友好性
- [ ] 使用查表法加速三角函数

### 数值稳定性
- [ ] 使用稳定的正交化算法（如改进的Gram-Schmidt）
- [ ] 避免在接近奇异的情况下求逆
- [ ] 使用适当的epsilon进行浮点比较
- [ ] 定期重新正交化旋转矩阵
- [ ] 考虑使用双精度进行关键计算