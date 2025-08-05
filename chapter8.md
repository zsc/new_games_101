# 第8章：全局光照

全局光照（Global Illumination）是计算机图形学中最具挑战性的问题之一。与局部光照不同，全局光照考虑了光线在场景中的多次反射、折射和散射，能够生成逼真的视觉效果，如软阴影、颜色渗透（color bleeding）、焦散（caustics）等。本章将从辐射度量学的基础出发，推导渲染方程，并介绍蒙特卡洛方法在求解全局光照中的应用。

## 8.1 辐射度量学

辐射度量学（Radiometry）是描述光能传播的物理框架。理解这些概念对于正确实现基于物理的渲染至关重要。

### 8.1.1 基本概念

**辐射能量（Radiant Energy）**
- 符号：$Q$
- 单位：焦耳（J）
- 光源发出的总能量

**辐射功率/辐射通量（Radiant Power/Flux）**
- 符号：$\Phi = \frac{dQ}{dt}$
- 单位：瓦特（W）或流明（lm）
- 单位时间内通过某个表面的能量

**辐射强度（Radiant Intensity）**
- 符号：$I = \frac{d\Phi}{d\omega}$
- 单位：W/sr（瓦特/球面度）
- 点光源在某个方向上的功率密度

**辐照度（Irradiance）**
- 符号：$E = \frac{d\Phi}{dA}$
- 单位：W/m²
- 单位面积接收到的功率

**辐射度（Radiance）**
- 符号：$L = \frac{d^2\Phi}{dA \cos\theta \, d\omega}$
- 单位：W/(m²·sr)
- 单位面积、单位立体角上的功率

### 8.1.2 辐射度的性质

辐射度是渲染中最重要的量，具有以下关键性质：

1. **方向性**：$L(\mathbf{p}, \omega)$ 依赖于位置 $\mathbf{p}$ 和方向 $\omega$

2. **不变性**：在真空中传播时，辐射度沿光线保持不变
   $$L(\mathbf{p}, \omega) = L(\mathbf{p} + t\omega, \omega)$$

3. **可加性**：多个光源的贡献可以线性叠加

### 8.1.3 立体角与投影立体角

**立体角**的微分形式：
$$d\omega = \sin\theta \, d\theta \, d\phi$$

在球坐标系中，对整个半球积分：
$$\int_{\Omega} d\omega = \int_0^{2\pi} \int_0^{\pi/2} \sin\theta \, d\theta \, d\phi = 2\pi$$

**投影立体角**：
$$d\omega^\perp = \cos\theta \, d\omega$$

这在计算辐照度时特别有用：
$$E(\mathbf{p}) = \int_{\Omega} L_i(\mathbf{p}, \omega) \cos\theta \, d\omega$$

### 8.1.4 BRDF与反射方程

双向反射分布函数（BRDF）定义了表面如何反射光线：

$$f_r(\omega_i \to \omega_o) = \frac{dL_o(\omega_o)}{dE(\omega_i)} = \frac{dL_o(\omega_o)}{L_i(\omega_i) \cos\theta_i \, d\omega_i}$$

BRDF的性质：
1. **非负性**：$f_r \geq 0$
2. **互易性**（Helmholtz互易原理）：$f_r(\omega_i \to \omega_o) = f_r(\omega_o \to \omega_i)$
3. **能量守恒**：$\int_{\Omega} f_r(\omega_i \to \omega_o) \cos\theta_o \, d\omega_o \leq 1$

反射方程描述了出射辐射度：
$$L_o(\mathbf{p}, \omega_o) = L_e(\mathbf{p}, \omega_o) + \int_{\Omega} f_r(\mathbf{p}, \omega_i \to \omega_o) L_i(\mathbf{p}, \omega_i) \cos\theta_i \, d\omega_i$$

## 8.2 渲染方程

### 8.2.1 渲染方程的推导

Kajiya在1986年提出的渲染方程是全局光照的数学基础。从反射方程出发，考虑入射辐射度来自其他表面的反射：

$$L_i(\mathbf{p}, \omega_i) = L_o(\mathbf{p}', -\omega_i)$$

其中 $\mathbf{p}'$ 是从 $\mathbf{p}$ 沿 $\omega_i$ 方向看到的最近表面点。

将此代入反射方程，得到完整的渲染方程：
$$L_o(\mathbf{p}, \omega_o) = L_e(\mathbf{p}, \omega_o) + \int_{\Omega} f_r(\mathbf{p}, \omega_i \to \omega_o) L_o(\mathbf{p}', -\omega_i) \cos\theta_i \, d\omega_i$$

### 8.2.2 积分算子形式

定义光传输算子 $\mathcal{K}$：
$$(\mathcal{K}L)(\mathbf{p}, \omega) = \int_{\Omega} f_r(\mathbf{p}, \omega_i \to \omega) L(\mathbf{p}', -\omega_i) \cos\theta_i \, d\omega_i$$

渲染方程可写为：
$$L = L_e + \mathcal{K}L$$

形式解为：
$$L = (I - \mathcal{K})^{-1}L_e = \sum_{n=0}^{\infty} \mathcal{K}^n L_e$$

这个级数展开具有物理意义：
- $\mathcal{K}^0 L_e = L_e$：直接光照
- $\mathcal{K}^1 L_e$：一次反射
- $\mathcal{K}^2 L_e$：二次反射
- ...

### 8.2.3 路径积分形式

考虑长度为 $n$ 的光路径 $\mathbf{x}_0 \mathbf{x}_1 ... \mathbf{x}_n$，其中 $\mathbf{x}_0$ 在相机，$\mathbf{x}_n$ 在光源。路径的贡献为：

$$L(\mathbf{x}_0 \to \mathbf{x}_1) = \int L_e(\mathbf{x}_n \to \mathbf{x}_{n-1}) \prod_{i=1}^{n-1} f_r(\mathbf{x}_{i+1} \to \mathbf{x}_i \to \mathbf{x}_{i-1}) G(\mathbf{x}_i \leftrightarrow \mathbf{x}_{i+1}) \, dA(\mathbf{x}_2) ... dA(\mathbf{x}_n)$$

其中几何因子：
$$G(\mathbf{x} \leftrightarrow \mathbf{y}) = \frac{\cos\theta_x \cos\theta_y}{||\mathbf{x} - \mathbf{y}||^2} V(\mathbf{x} \leftrightarrow \mathbf{y})$$

$V$ 是可见性函数（0或1）。

### 8.2.4 测度变换

从立体角测度到面积测度的转换：
$$d\omega = \frac{\cos\theta' \, dA'}{||\mathbf{p} - \mathbf{p}'||^2}$$

这使得渲染方程可以在面积域上积分：
$$L_o(\mathbf{p}, \omega_o) = L_e(\mathbf{p}, \omega_o) + \int_{\mathcal{M}} f_r(\mathbf{p}, \mathbf{p}' \to \omega_o) L_o(\mathbf{p}', \mathbf{p} - \mathbf{p}') G(\mathbf{p} \leftrightarrow \mathbf{p}') \, dA'$$

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