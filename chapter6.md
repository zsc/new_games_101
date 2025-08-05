# 第6章：几何表示

本章深入探讨计算机图形学中几何对象的各种表示方法。我们将从基本的几何表示开始，逐步深入到参数化曲线曲面理论，最后讨论现代图形学中的网格处理技术和阴影图算法。对于AI科学家而言，理解这些几何表示不仅有助于3D视觉任务，还能为神经隐式表示（Neural Implicit Representations）等前沿技术打下坚实基础。

## 6.1 基本几何表示方法

### 6.1.1 隐式表示与显式表示

几何对象的数学表示主要分为两大类：

**隐式表示（Implicit Representation）**：通过方程 $f(\mathbf{x}) = 0$ 定义几何对象，其中 $\mathbf{x} \in \mathbb{R}^n$。

优点：
- 易于判断点的内外关系：$f(\mathbf{x}) < 0$ 表示内部，$f(\mathbf{x}) > 0$ 表示外部
- 便于进行布尔运算（并、交、差）
- 自然支持实体建模

缺点：
- 难以参数化采样
- 不易控制局部形状

**显式表示（Explicit Representation）**：直接给出几何对象上点的坐标，如参数化表示 $\mathbf{x} = \mathbf{g}(u, v)$。

优点：
- 易于采样和遍历
- 便于纹理映射
- 直观的形状控制

缺点：
- 难以判断点的内外关系
- 拓扑变化处理复杂

### 6.1.2 代数表面

代数表面是由多项式方程定义的隐式表面。对于度数为 $d$ 的代数表面：

$$f(\mathbf{x}) = \sum_{i+j+k \leq d} a_{ijk} x^i y^j z^k = 0$$

**二次曲面（Quadrics）** 是最常用的代数表面（$d=2$）：

$$\mathbf{x}^T \mathbf{Q} \mathbf{x} + 2\mathbf{b}^T \mathbf{x} + c = 0$$

其中 $\mathbf{Q}$ 是 $3 \times 3$ 对称矩阵。通过特征值分解 $\mathbf{Q} = \mathbf{V}\mathbf{\Lambda}\mathbf{V}^T$，可以分类二次曲面：
- 椭球体：三个正特征值
- 单叶双曲面：两正一负
- 双叶双曲面：一正两负
- 椭圆抛物面：两个非零特征值
- 圆柱面/圆锥面：一个零特征值

### 6.1.3 构造实体几何（CSG）

CSG通过布尔运算组合简单几何体构建复杂形状。对于两个隐式表示的对象 $A: f_A(\mathbf{x}) \leq 0$ 和 $B: f_B(\mathbf{x}) \leq 0$：

- 并集：$A \cup B: \min(f_A(\mathbf{x}), f_B(\mathbf{x})) \leq 0$
- 交集：$A \cap B: \max(f_A(\mathbf{x}), f_B(\mathbf{x})) \leq 0$
- 差集：$A - B: \max(f_A(\mathbf{x}), -f_B(\mathbf{x})) \leq 0$

**平滑布尔运算**：为避免尖锐边缘，可使用平滑最小/最大函数：

$$\text{smin}(a, b, k) = -\frac{1}{k}\log(e^{-ka} + e^{-kb})$$

### 6.1.4 符号距离场（SDF）

SDF是一种特殊的隐式表示，其中 $f(\mathbf{x})$ 表示点 $\mathbf{x}$ 到表面的有符号距离：

$$f(\mathbf{x}) = \text{sign}(\mathbf{x}) \cdot \min_{\mathbf{y} \in \partial S} \|\mathbf{x} - \mathbf{y}\|$$

SDF的关键性质：
- 梯度为单位向量：$\|\nabla f\| = 1$（几乎处处成立）
- 法向量：$\mathbf{n} = \nabla f$
- 最近点投影：$\mathbf{p}_{\text{closest}} = \mathbf{x} - f(\mathbf{x})\nabla f(\mathbf{x})$

**Eikonal方程**：SDF满足偏微分方程 $\|\nabla f\| = 1$，可通过快速行进法（Fast Marching Method）或快速扫描法（Fast Sweeping Method）求解。

### 6.1.5 点云与体素表示

**点云表示**：
- 无连接信息的3D点集合 $\{\mathbf{p}_i\}_{i=1}^N$
- 常用于3D扫描和激光雷达数据
- 表面重建方法：泊松重建、移动最小二乘（MLS）

**体素表示**：
- 3D网格中的占用信息
- 层次结构：八叉树（Octree）表示
- 稀疏体素八叉树（SVO）用于大规模场景

### 6.1.6 网格表示

多边形网格（主要是三角网格）是实时图形学的主流表示：

**半边数据结构（Half-Edge）**：
```
struct HalfEdge {
    Vertex* vertex;      // 指向的顶点
    HalfEdge* twin;      // 对偶半边
    HalfEdge* next;      // 下一条半边
    Face* face;          // 所属面片
};
```

**欧拉公式**：对于连通的多面体网格：
$$V - E + F = 2(1 - g)$$
其中 $V$ 是顶点数，$E$ 是边数，$F$ 是面数，$g$ 是亏格（genus）。

## 6.2 曲线与曲面

### 6.2.1 参数曲线基础

参数曲线 $\mathbf{c}(t): [a, b] \rightarrow \mathbb{R}^3$ 的基本微分几何：

**切向量**：$\mathbf{t}(t) = \frac{d\mathbf{c}}{dt}$

**弧长参数化**：$s(t) = \int_a^t \|\mathbf{c}'(\tau)\| d\tau$

**曲率**：$\kappa = \frac{\|\mathbf{t}' \times \mathbf{t}\|}{\|\mathbf{t}\|^3}$

**Frenet-Serret标架**：
- 切向量：$\mathbf{T} = \frac{\mathbf{c}'}{|\mathbf{c}'|}$
- 法向量：$\mathbf{N} = \frac{\mathbf{T}'}{|\mathbf{T}'|}$
- 副法向量：$\mathbf{B} = \mathbf{T} \times \mathbf{N}$

Frenet-Serret公式：
$$\begin{pmatrix} \mathbf{T}' \\ \mathbf{N}' \\ \mathbf{B}' \end{pmatrix} = \begin{pmatrix} 0 & \kappa & 0 \\ -\kappa & 0 & \tau \\ 0 & -\tau & 0 \end{pmatrix} \begin{pmatrix} \mathbf{T} \\ \mathbf{N} \\ \mathbf{B} \end{pmatrix}$$

### 6.2.2 贝塞尔曲线

$n$ 次贝塞尔曲线由 $n+1$ 个控制点 $\{\mathbf{P}_i\}_{i=0}^n$ 定义：

$$\mathbf{B}(t) = \sum_{i=0}^n B_i^n(t) \mathbf{P}_i$$

其中 $B_i^n(t) = \binom{n}{i}t^i(1-t)^{n-i}$ 是Bernstein基函数。

**性质**：
- 端点插值：$\mathbf{B}(0) = \mathbf{P}_0$，$\mathbf{B}(1) = \mathbf{P}_n$
- 凸包性质：曲线完全位于控制点的凸包内
- 仿射不变性
- 变差减少性质

**de Casteljau算法**：递归计算 $\mathbf{B}(t)$：
$$\mathbf{P}_i^{(k)}(t) = (1-t)\mathbf{P}_i^{(k-1)} + t\mathbf{P}_{i+1}^{(k-1)}$$

**导数**：
$$\mathbf{B}'(t) = n\sum_{i=0}^{n-1} B_i^{n-1}(t)(\mathbf{P}_{i+1} - \mathbf{P}_i)$$

### 6.2.3 B样条曲线

B样条通过局部基函数实现局部控制：

$$\mathbf{C}(t) = \sum_{i=0}^{n} N_{i,p}(t) \mathbf{P}_i$$

其中 $N_{i,p}(t)$ 是 $p$ 次B样条基函数，由Cox-de Boor递归定义：

$$N_{i,0}(t) = \begin{cases} 1 & t_i \leq t < t_{i+1} \\ 0 & \text{否则} \end{cases}$$

$$N_{i,p}(t) = \frac{t - t_i}{t_{i+p} - t_i}N_{i,p-1}(t) + \frac{t_{i+p+1} - t}{t_{i+p+1} - t_{i+1}}N_{i+1,p-1}(t)$$

**节点向量**：$\{t_0, t_1, ..., t_{m}\}$，其中 $m = n + p + 1$

**性质**：
- 局部支撑：$N_{i,p}(t) = 0$ 当 $t \notin [t_i, t_{i+p+1})$
- $C^{p-1}$ 连续性（简单节点情况）
- 分段多项式

### 6.2.4 NURBS（非均匀有理B样条）

NURBS通过引入权重和有理形式扩展B样条：

$$\mathbf{C}(t) = \frac{\sum_{i=0}^{n} w_i N_{i,p}(t) \mathbf{P}_i}{\sum_{i=0}^{n} w_i N_{i,p}(t)}$$

**优势**：
- 精确表示二次曲线（圆、椭圆、抛物线）
- 投影不变性
- 统一的数学框架

### 6.2.5 参数曲面

双参数曲面 $\mathbf{S}(u,v): D \subset \mathbb{R}^2 \rightarrow \mathbb{R}^3$

**第一基本形式**：度量曲面上的距离
$$I = E du^2 + 2F du dv + G dv^2$$
其中：
- $E = \mathbf{S}_u \cdot \mathbf{S}_u$
- $F = \mathbf{S}_u \cdot \mathbf{S}_v$
- $G = \mathbf{S}_v \cdot \mathbf{S}_v$

**法向量**：
$$\mathbf{n} = \frac{\mathbf{S}_u \times \mathbf{S}_v}{|\mathbf{S}_u \times \mathbf{S}_v|}$$

**第二基本形式**：描述曲面的弯曲
$$II = L du^2 + 2M du dv + N dv^2$$
其中：
- $L = \mathbf{S}_{uu} \cdot \mathbf{n}$
- $M = \mathbf{S}_{uv} \cdot \mathbf{n}$
- $N = \mathbf{S}_{vv} \cdot \mathbf{n}$

**主曲率**：Weingarten矩阵的特征值
$$\mathbf{W} = \begin{pmatrix} L & M \\ M & N \end{pmatrix} \begin{pmatrix} E & F \\ F & G \end{pmatrix}^{-1}$$

**高斯曲率**：$K = \kappa_1 \kappa_2 = \frac{LN - M^2}{EG - F^2}$

**平均曲率**：$H = \frac{\kappa_1 + \kappa_2}{2} = \frac{EN - 2FM + GL}{2(EG - F^2)}$

### 6.2.6 细分曲面

细分方案通过递归细化网格生成光滑曲面。

**Loop细分**（三角网格）：
- 新顶点位置（边中点）：
  $$\mathbf{v}_{\text{new}} = \frac{3}{8}(\mathbf{v}_1 + \mathbf{v}_2) + \frac{1}{8}(\mathbf{v}_3 + \mathbf{v}_4)$$
- 旧顶点更新：
  $$\mathbf{v}_{\text{updated}} = (1 - n\beta)\mathbf{v}_{\text{old}} + \beta\sum_{i=1}^n \mathbf{v}_i$$
  其中 $\beta = \frac{1}{n}(\frac{5}{8} - (\frac{3}{8} + \frac{1}{4}\cos\frac{2\pi}{n})^2)$

**Catmull-Clark细分**（四边形网格）：
- 面点：面内顶点的平均
- 边点：相邻两个顶点和两个面点的平均
- 顶点更新：基于价数的加权平均

极限曲面在正则顶点处为 $C^2$ 连续，在非正则顶点处为 $C^1$ 连续。

## 6.3 网格处理与阴影图

### 6.3.1 网格简化

网格简化旨在减少多边形数量同时保持视觉质量。

**边收缩（Edge Collapse）**：
将边 $(v_i, v_j)$ 收缩到新顶点 $\bar{v}$。使用二次误差度量（QEM）：

每个顶点关联误差二次型：
$$Q_v = \sum_{p \in \text{planes}(v)} K_p$$

其中 $K_p = \mathbf{n}_p\mathbf{n}_p^T$ 对于平面 $\mathbf{n}_p^T\mathbf{x} + d = 0$。

边收缩代价：
$$\text{cost}(v_i, v_j) = \min_{\bar{v}} \bar{v}^T(Q_i + Q_j)\bar{v}$$

最优位置通过求解线性系统获得：
$$\bar{v} = -(Q_i + Q_j)^{-1}\mathbf{b}$$

**顶点聚类**：
将空间划分为网格，每个单元内的顶点合并为一个代表顶点。

### 6.3.2 网格参数化

将3D曲面映射到2D参数域 $\phi: M \rightarrow \mathbb{R}^2$。

**调和映射（Harmonic Mapping）**：
最小化Dirichlet能量：
$$E_D = \frac{1}{2}\int_M |\nabla \phi|^2 dA$$

离散化后得到线性系统：
$$\sum_{j \in N(i)} w_{ij}(u_j - u_i) = 0$$

权重选择：
- 均匀权重：$w_{ij} = 1$
- 余切权重：$w_{ij} = \cot\alpha_{ij} + \cot\beta_{ij}$
- 平均值坐标：考虑角度和距离

**保角参数化（Conformal Parameterization）**：
最小化角度畸变，使用复数表示和Cauchy-Riemann方程。

**LSCM（最小二乘保角映射）**：
$$\min_{\phi} \sum_{T} A_T ||\nabla \phi - R_{90°}\nabla \phi||^2$$

### 6.3.3 网格平滑

**拉普拉斯平滑**：
$$\mathbf{v}_i^{new} = \mathbf{v}_i + \lambda \mathbf{L}(\mathbf{v}_i)$$

其中离散拉普拉斯算子：
$$\mathbf{L}(\mathbf{v}_i) = \frac{1}{|N(i)|}\sum_{j \in N(i)} (\mathbf{v}_j - \mathbf{v}_i)$$

**双边滤波扩展**：
$$\mathbf{v}_i^{new} = \mathbf{v}_i + \lambda \sum_{j \in N(i)} w_s(\||\mathbf{v}_j - \mathbf{v}_i||) w_n(\mathbf{n}_i \cdot \mathbf{n}_j) (\mathbf{v}_j - \mathbf{v}_i)$$

### 6.3.4 网格修复

**孔洞填充**：
1. 检测边界环
2. 三角化边界多边形（最小角度、最小面积）
3. 细分和平滑新增三角形

**自相交检测与修复**：
- 使用空间划分结构加速相交测试
- 局部重新网格化解决相交

### 6.3.5 阴影图（Shadow Mapping）

阴影图是实时渲染中生成阴影的标准技术。

**基本算法**：
1. 从光源视角渲染场景，存储深度
2. 从相机视角渲染，测试每个片元是否在阴影中

**深度比较**：
对于世界空间点 $\mathbf{p}$：
1. 变换到光源空间：$\mathbf{p}_{light} = \mathbf{M}_{light}\mathbf{p}$
2. 采样阴影图：$z_{stored} = \text{ShadowMap}(p_{light}.xy)$
3. 比较深度：$\text{inShadow} = (p_{light}.z > z_{stored} + \text{bias})$

**深度偏移（Bias）问题**：
- 恒定偏移：$\text{bias} = \epsilon$
- 斜率偏移：$\text{bias} = \epsilon \cdot \tan(\theta)$
- 法线偏移：沿法线方向偏移采样点

### 6.3.6 级联阴影图（CSM）

将视锥体分割为多个子区域，每个使用独立的阴影图。

**分割策略**：
- 均匀分割：$z_i = z_{near} + i \cdot \frac{z_{far} - z_{near}}{N}$
- 对数分割：$z_i = z_{near} \cdot (z_{far}/z_{near})^{i/N}$
- 实用分割：两者的混合

**层级选择**：
```glsl
float depth = length(viewPos);
int cascade = 0;
for (int i = 0; i < NUM_CASCADES - 1; i++) {
    if (depth > cascadeSplits[i]) cascade = i + 1;
}
```

### 6.3.7 软阴影技术

**PCF（Percentage Closer Filtering）**：
$$\text{shadow} = \frac{1}{N}\sum_{i=1}^N \text{ShadowTest}(\mathbf{p} + \delta_i)$$

**PCSS（Percentage Closer Soft Shadows）**：
1. 阻挡物搜索：估计平均阻挡物深度
2. 半影估计：$w_{penumbra} = (d_{receiver} - d_{blocker}) \cdot w_{light} / d_{blocker}$
3. PCF过滤：使用估计的半影大小

**方差阴影图（VSM）**：
存储深度和深度平方，使用Chebyshev不等式估计阴影概率：
$$P(z \geq t) \leq \frac{\sigma^2}{\sigma^2 + (t - \mu)^2}$$

### 6.3.8 光线追踪阴影

**硬阴影**：
发射阴影光线从着色点到光源：
$$\text{visible} = \begin{cases} 1 & \text{无遮挡} \\ 0 & \text{有遮挡} \end{cases}$$

**软阴影**：
对面光源采样多条光线：
$$L = \frac{1}{N}\sum_{i=1}^N V(\mathbf{p}, \mathbf{x}_i) \cdot L_e(\mathbf{x}_i) \cdot G(\mathbf{p}, \mathbf{x}_i)$$

其中 $V$ 是可见性函数，$G$ 是几何项。

### 6.3.9 环境光遮蔽（AO）

**屏幕空间环境光遮蔽（SSAO）**：
$$AO(\mathbf{p}) = \frac{1}{N}\sum_{i=1}^N \text{visibility}(\mathbf{p}, \mathbf{p} + \mathbf{s}_i)$$

采样点在半球内随机分布，使用深度缓冲判断遮挡。

**预计算AO**：
- 每顶点存储遮蔽值
- 使用光线投射计算
- 考虑法线加权：$AO = \frac{1}{\pi}\int_{\Omega} V(\omega) \cdot (\mathbf{n} \cdot \omega)^+ d\omega$

## 本章小结

本章系统介绍了计算机图形学中的几何表示方法：

**基本几何表示**：
- 隐式表示：$f(\mathbf{x}) = 0$，便于内外判断和布尔运算
- 显式表示：参数化形式，便于采样和纹理映射
- 符号距离场（SDF）：特殊的隐式表示，梯度为单位向量
- CSG：通过布尔运算构建复杂几何

**曲线与曲面理论**：
- 贝塞尔曲线：$\mathbf{B}(t) = \sum_{i=0}^n B_i^n(t) \mathbf{P}_i$
- B样条：局部控制，$C^{p-1}$ 连续性
- NURBS：有理形式，精确表示二次曲线
- 微分几何：第一/第二基本形式，高斯曲率 $K$，平均曲率 $H$
- 细分曲面：Loop（三角网格）、Catmull-Clark（四边形网格）

**网格处理与阴影技术**：
- 网格简化：QEM边收缩算法
- 网格参数化：调和映射、LSCM
- 阴影图：深度比较 + 偏移处理
- 级联阴影图：多分辨率解决透视走样
- 软阴影：PCF、PCSS、VSM
- 环境光遮蔽：SSAO、预计算AO

## 练习题

### 基础题

**练习 6.1**：给定一个椭球的隐式方程 $\frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} = 1$，推导其在点 $(x_0, y_0, z_0)$ 处的法向量。

*提示*：法向量是梯度方向 $\nabla f$。

<details>
<summary>答案</summary>

对隐式函数 $f(x,y,z) = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} - 1$ 求梯度：

$$\nabla f = \left(\frac{2x}{a^2}, \frac{2y}{b^2}, \frac{2z}{c^2}\right)$$

在点 $(x_0, y_0, z_0)$ 处的单位法向量：

$$\mathbf{n} = \frac{1}{\sqrt{\frac{4x_0^2}{a^4} + \frac{4y_0^2}{b^4} + \frac{4z_0^2}{c^4}}} \left(\frac{2x_0}{a^2}, \frac{2y_0}{b^2}, \frac{2z_0}{c^2}\right)$$

</details>

**练习 6.2**：证明三次贝塞尔曲线的中点 $\mathbf{B}(0.5)$ 可以表示为：
$$\mathbf{B}(0.5) = \frac{1}{8}(\mathbf{P}_0 + 3\mathbf{P}_1 + 3\mathbf{P}_2 + \mathbf{P}_3)$$

*提示*：代入 $t = 0.5$ 到贝塞尔曲线公式。

<details>
<summary>答案</summary>

三次贝塞尔曲线：
$$\mathbf{B}(t) = (1-t)^3\mathbf{P}_0 + 3t(1-t)^2\mathbf{P}_1 + 3t^2(1-t)\mathbf{P}_2 + t^3\mathbf{P}_3$$

代入 $t = 0.5$：
- $(1-0.5)^3 = 0.125 = \frac{1}{8}$
- $3 \cdot 0.5 \cdot (0.5)^2 = 3 \cdot 0.5 \cdot 0.25 = 0.375 = \frac{3}{8}$
- $3 \cdot (0.5)^2 \cdot 0.5 = 3 \cdot 0.25 \cdot 0.5 = 0.375 = \frac{3}{8}$
- $(0.5)^3 = 0.125 = \frac{1}{8}$

因此：
$$\mathbf{B}(0.5) = \frac{1}{8}\mathbf{P}_0 + \frac{3}{8}\mathbf{P}_1 + \frac{3}{8}\mathbf{P}_2 + \frac{1}{8}\mathbf{P}_3 = \frac{1}{8}(\mathbf{P}_0 + 3\mathbf{P}_1 + 3\mathbf{P}_2 + \mathbf{P}_3)$$

</details>

**练习 6.3**：对于一个三角形网格，已知顶点数 $V = 100$，面数 $F = 196$，使用欧拉公式计算边数 $E$。假设网格是单连通的（genus = 0）。

*提示*：欧拉公式 $V - E + F = 2(1 - g)$。

<details>
<summary>答案</summary>

对于单连通网格（$g = 0$）：
$$V - E + F = 2(1 - 0) = 2$$

代入已知值：
$$100 - E + 196 = 2$$
$$296 - E = 2$$
$$E = 294$$

验证：对于三角网格，每个面有3条边，每条边被2个面共享，所以 $E = \frac{3F}{2} = \frac{3 \times 196}{2} = 294$ ✓

</details>

**练习 6.4**：给定参数曲面 $\mathbf{S}(u,v) = (u\cos v, u\sin v, u^2)$，计算在点 $(u,v) = (1, 0)$ 处的法向量。

*提示*：计算 $\mathbf{S}_u \times \mathbf{S}_v$。

<details>
<summary>答案</summary>

首先计算偏导数：
$$\mathbf{S}_u = (\cos v, \sin v, 2u)$$
$$\mathbf{S}_v = (-u\sin v, u\cos v, 0)$$

在 $(u,v) = (1,0)$ 处：
$$\mathbf{S}_u(1,0) = (1, 0, 2)$$
$$\mathbf{S}_v(1,0) = (0, 1, 0)$$

叉积：
$$\mathbf{n} = \mathbf{S}_u \times \mathbf{S}_v = \begin{vmatrix} \mathbf{i} & \mathbf{j} & \mathbf{k} \\ 1 & 0 & 2 \\ 0 & 1 & 0 \end{vmatrix} = (-2, 0, 1)$$

单位法向量：
$$\hat{\mathbf{n}} = \frac{(-2, 0, 1)}{\sqrt{4 + 0 + 1}} = \frac{1}{\sqrt{5}}(-2, 0, 1)$$

</details>

### 挑战题

**练习 6.5**：证明对于符号距离场（SDF），如果 $f$ 是到曲面 $S$ 的符号距离函数，则 $\|\nabla f\| = 1$ 几乎处处成立。

*提示*：考虑最近点投影和距离函数的定义。

<details>
<summary>答案</summary>

设 $\mathbf{x}$ 是空间中一点，$\mathbf{p}$ 是 $\mathbf{x}$ 在曲面 $S$ 上的最近点。

定义：$f(\mathbf{x}) = \text{sign}(\mathbf{x}) \cdot \|\mathbf{x} - \mathbf{p}\|$

对于曲面外的点（$f > 0$），沿任意方向 $\mathbf{v}$ 的方向导数：
$$\frac{\partial f}{\partial \mathbf{v}} = \lim_{h \to 0} \frac{f(\mathbf{x} + h\mathbf{v}) - f(\mathbf{x})}{h}$$

由于 $\mathbf{p}$ 是最近点，当 $h$ 足够小时，$\mathbf{x} + h\mathbf{v}$ 的最近点仍是 $\mathbf{p}$（除非 $\mathbf{x}$ 在中轴上）。

因此：
$$f(\mathbf{x} + h\mathbf{v}) \approx \|\mathbf{x} + h\mathbf{v} - \mathbf{p}\|$$

使用泰勒展开：
$$\|\mathbf{x} + h\mathbf{v} - \mathbf{p}\| = \|\mathbf{x} - \mathbf{p}\| + h\frac{(\mathbf{x} - \mathbf{p}) \cdot \mathbf{v}}{\|\mathbf{x} - \mathbf{p}\|} + O(h^2)$$

所以：
$$\frac{\partial f}{\partial \mathbf{v}} = \frac{(\mathbf{x} - \mathbf{p}) \cdot \mathbf{v}}{\|\mathbf{x} - \mathbf{p}\|}$$

梯度为：
$$\nabla f = \frac{\mathbf{x} - \mathbf{p}}{\|\mathbf{x} - \mathbf{p}\|}$$

因此 $\|\nabla f\| = 1$。

注意：在中轴（medial axis）上，最近点不唯一，梯度不存在。

</details>

**练习 6.6**：推导Loop细分中 $\beta$ 值的选择，使得细分曲面在正则顶点处达到 $C^2$ 连续性。

*提示*：分析细分矩阵的特征值。

<details>
<summary>答案</summary>

对于价数为 $n$ 的顶点，Loop细分的更新规则：
$$\mathbf{v}_{\text{new}} = (1 - n\beta)\mathbf{v}_{\text{old}} + \beta\sum_{i=1}^n \mathbf{v}_i$$

细分矩阵 $\mathbf{S}$ 的特征值决定了曲面的连续性。

对于 $C^2$ 连续性，需要：
1. 最大特征值 $\lambda_0 = 1$（保持仿射不变性）
2. 次大特征值 $\lambda_1 = \lambda_2 < 1$（收敛性）
3. $\lambda_1/\lambda_0 > \lambda_3/\lambda_1$（$C^2$ 条件）

通过傅里叶分析，特征值为：
$$\lambda_k = \frac{1}{8}(5 - 3\cos\frac{2\pi k}{n}) + \beta(1 - \cos\frac{2\pi k}{n})$$

要使 $\lambda_1 = \lambda_2$，需要旋转对称性。

最优的 $\beta$ 值（Warren's choice）：
$$\beta = \frac{1}{n}\left(\frac{5}{8} - \left(\frac{3}{8} + \frac{1}{4}\cos\frac{2\pi}{n}\right)^2\right)$$

这保证了：
- $n = 3$：$\beta = 3/16$
- $n = 6$：$\beta = 1/16$（正则顶点）
- 大 $n$：$\beta \approx 3/(8n)$

</details>

**练习 6.7**：设计一个算法，将任意多边形网格转换为符号距离场表示。讨论算法的时间复杂度和精度权衡。

*提示*：考虑快速行进法（FMM）或快速扫描法。

<details>
<summary>答案</summary>

**算法设计**：

1. **初始化**：
   - 在规则网格上采样
   - 标记网格点的内外关系（射线法）
   - 计算到最近三角形的距离

2. **精确距离计算**（边界附近）：
   - 对每个网格点，找最近的三角形
   - 计算点到三角形的距离（考虑顶点、边、面）
   - 使用空间划分加速（KD-tree或BVH）

3. **快速行进法**（远离边界）：
   - 从已知距离的点开始传播
   - 求解Eikonal方程：$\|\nabla f\| = 1$
   - 使用堆维护波前

4. **符号确定**：
   - 奇偶规则或缠绕数
   - 法向量一致性检查

**时间复杂度**：
- 初始化：$O(N \log M)$，$N$ 是网格点数，$M$ 是三角形数
- FMM传播：$O(N \log N)$
- 总体：$O(N \log N + N \log M)$

**精度权衡**：
- 网格分辨率 vs. 内存消耗
- 窄带宽度 vs. 计算时间
- 插值阶数 vs. 平滑度

**优化策略**：
- 自适应网格细化
- GPU并行化（特别适合快速扫描法）
- 层次化表示（稀疏体素八叉树）

</details>

**练习 6.8**：分析PCSS（Percentage Closer Soft Shadows）算法中，光源大小、遮挡物距离和接收器距离如何影响半影大小。推导半影估计公式。

*提示*：使用相似三角形原理。

<details>
<summary>答案</summary>

**几何分析**：

设置：
- 光源宽度：$w_{light}$
- 遮挡物到光源距离：$d_{blocker}$
- 接收器到光源距离：$d_{receiver}$
- 遮挡物到接收器距离：$d_{receiver} - d_{blocker}$

**半影大小推导**：

使用相似三角形：
$$\frac{w_{penumbra}}{d_{receiver} - d_{blocker}} = \frac{w_{light}}{d_{blocker}}$$

解得：
$$w_{penumbra} = \frac{(d_{receiver} - d_{blocker}) \cdot w_{light}}{d_{blocker}}$$

**影响因素分析**：

1. **光源大小** $w_{light}$：
   - 线性关系：光源越大，半影越宽
   - 点光源（$w_{light} \to 0$）：硬阴影

2. **遮挡物距离** $d_{blocker}$：
   - 反比关系：遮挡物越近光源，半影越大
   - 极限情况：$d_{blocker} \to 0$，$w_{penumbra} \to \infty$

3. **接收器距离** $d_{receiver}$：
   - 正相关：接收器越远，半影越大
   - 相对距离 $(d_{receiver} - d_{blocker})$ 是关键

**PCSS算法步骤**：

1. **搜索遮挡物**：
   ```glsl
   float searchRadius = lightSize * (receiverDepth - nearPlane) / receiverDepth;
   float avgBlockerDepth = findAverageBlockerDepth(shadowCoord, searchRadius);
   ```

2. **估计半影**：
   ```glsl
   float penumbraWidth = (receiverDepth - avgBlockerDepth) * lightSize / avgBlockerDepth;
   ```

3. **PCF过滤**：
   ```glsl
   float shadow = PCF(shadowCoord, penumbraWidth);
   ```

**误差来源**：
- 单一平均遮挡物深度的近似
- 离散采样导致的噪声
- 透视投影的非线性

</details>

## 常见陷阱与错误

1. **符号距离场的梯度假设**
   - 错误：认为SDF的梯度处处为1
   - 正确：在中轴（medial axis）上梯度不连续

2. **贝塞尔曲线的参数化**
   - 错误：假设参数 $t$ 对应弧长
   - 正确：贝塞尔曲线通常不是弧长参数化的

3. **网格简化的拓扑保持**
   - 错误：边收缩可能改变拓扑（创建非流形结构）
   - 正确：需要检查收缩操作的有效性

4. **阴影图的深度精度**
   - 错误：忽略深度缓冲的有限精度
   - 正确：使用适当的偏移和高精度深度格式

5. **细分曲面的边界处理**
   - 错误：对边界顶点使用内部规则
   - 正确：边界需要特殊的细分规则

6. **NURBS权重的影响**
   - 错误：认为权重只是简单的缩放
   - 正确：权重影响曲线的参数化速度

7. **CSG运算的数值稳定性**
   - 错误：直接使用 min/max 进行布尔运算
   - 正确：在边界附近需要特殊处理

8. **环境光遮蔽的法线偏移**
   - 错误：不考虑自遮挡问题
   - 正确：沿法线偏移采样起点

## 最佳实践检查清单

### 几何表示选择
- [ ] 根据应用需求选择合适的表示（隐式/显式）
- [ ] 考虑内存占用和查询效率的平衡
- [ ] 评估修改操作的频率和复杂度
- [ ] 确保表示支持所需的查询操作

### 曲线曲面设计
- [ ] 选择适当的连续性要求（$C^0$/$C^1$/$C^2$）
- [ ] 验证控制点数量与曲线复杂度的匹配
- [ ] 检查参数域的有效范围
- [ ] 考虑数值精度和稳定性

### 网格处理
- [ ] 保持网格的流形性质
- [ ] 验证欧拉特征数的一致性
- [ ] 使用半边结构进行拓扑操作
- [ ] 实现鲁棒的相交测试

### 阴影实现
- [ ] 选择合适的阴影图分辨率
- [ ] 正确设置深度偏移参数
- [ ] 实现阴影图的过滤（PCF等）
- [ ] 考虑级联阴影图的必要性

### 性能优化
- [ ] 使用空间数据结构加速查询
- [ ] 实现细节层次（LOD）系统
- [ ] 利用GPU并行计算
- [ ] 缓存重复计算的结果

### 数值稳定性
- [ ] 避免除零和接近奇异的情况
- [ ] 使用稳定的算法（如QR分解代替法线方程）
- [ ] 选择合适的浮点精度
- [ ] 实现容错机制
