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

#### 参数化曲线与弧长参数化

在实际应用中，我们经常需要控制物体沿曲线的运动速度。考虑参数曲线 $\mathbf{r}(u)$，其中 $u \in [0,1]$。弧长定义为：

$$s(u) = \int_0^u \|\mathbf{r}'(\tau)\| d\tau$$

理想的弧长参数化满足 $\|\mathbf{r}'(s)\| = 1$，这保证了匀速运动。然而，对于大多数曲线，弧长积分没有解析解。实践中使用数值方法：

1. **查找表方法**：预计算 $(u_i, s_i)$ 对，使用插值
2. **Newton-Raphson迭代**：求解 $s(u) - s_{target} = 0$

#### 时间扭曲与缓动函数

为了创建更自然的动画，我们使用时间扭曲函数 $\tau = w(t)$：

$$\mathbf{p}(t) = \mathbf{q}(w(t))$$

常用的缓动函数包括：

**缓入（Ease-in）**：$w(t) = t^n$，$n > 1$

**缓出（Ease-out）**：$w(t) = 1 - (1-t)^n$

**缓入缓出（Ease-in-out）**：
$$w(t) = \begin{cases}
\frac{1}{2}(2t)^n & t < 0.5 \\
1 - \frac{1}{2}(2(1-t))^n & t \geq 0.5
\end{cases}$$

**Hermite基缓动**：使用三次Hermite曲线，控制端点导数：
$$w(t) = h_{00}(t) + h_{10}(t)v_0 + h_{01}(t) + h_{11}(t)v_1$$

其中Hermite基函数：
$$h_{00}(t) = 2t^3 - 3t^2 + 1, \quad h_{10}(t) = t^3 - 2t^2 + t$$
$$h_{01}(t) = -2t^3 + 3t^2, \quad h_{11}(t) = t^3 - t^2$$

#### 高阶运动分析

除了位置、速度和加速度，高阶导数也有物理意义：

**急动度（Jerk）**：$\mathbf{j}(t) = \frac{d^3\mathbf{x}}{dt^3}$

**急动度变化率（Snap）**：$\mathbf{s}(t) = \frac{d^4\mathbf{x}}{dt^4}$

在运动规划中，最小化急动度可以产生更平滑的轨迹：

$$\min \int_0^T \|\mathbf{j}(t)\|^2 dt$$

这导致七次多项式轨迹，满足位置、速度、加速度和急动度的边界条件。

### 11.1.2 关键帧动画与插值

关键帧动画是最基础也是最常用的动画技术。给定一系列关键帧 $\{(\mathbf{p}_i, t_i)\}_{i=0}^n$，我们需要在这些关键帧之间进行插值。

**线性插值**：
$$\mathbf{p}(t) = \mathbf{p}_i + \frac{t - t_i}{t_{i+1} - t_i}(\mathbf{p}_{i+1} - \mathbf{p}_i), \quad t \in [t_i, t_{i+1}]$$

线性插值虽然简单，但会在关键帧处产生速度不连续。参数 $s = \frac{t - t_i}{t_{i+1} - t_i}$ 称为插值参数。

#### 旋转插值的深入分析

**球面线性插值（SLERP）**：
对于四元数 $\mathbf{q}_i$ 和 $\mathbf{q}_{i+1}$：

$$\text{slerp}(\mathbf{q}_i, \mathbf{q}_{i+1}, s) = \frac{\sin((1-s)\theta)}{\sin\theta}\mathbf{q}_i + \frac{\sin(s\theta)}{\sin\theta}\mathbf{q}_{i+1}$$

其中 $\cos\theta = \mathbf{q}_i \cdot \mathbf{q}_{i+1}$，$s = \frac{t - t_i}{t_{i+1} - t_i}$。

**数值稳定性考虑**：
当 $\theta \approx 0$ 时，$\sin\theta \approx 0$，导致数值不稳定。此时使用线性插值近似：

$$\text{slerp}(\mathbf{q}_i, \mathbf{q}_{i+1}, s) \approx (1-s)\mathbf{q}_i + s\mathbf{q}_{i+1}, \quad \text{当} \theta < \epsilon$$

**最短路径问题**：
四元数 $\mathbf{q}$ 和 $-\mathbf{q}$ 表示相同旋转。为确保最短路径插值：

$$\text{如果} \quad \mathbf{q}_i \cdot \mathbf{q}_{i+1} < 0, \quad \text{则} \quad \mathbf{q}_{i+1} \leftarrow -\mathbf{q}_{i+1}$$

#### 高阶插值方法

**立方Hermite插值**：
给定端点值和导数：
$$\mathbf{p}(t) = h_{00}(s)\mathbf{p}_i + h_{10}(s)(t_{i+1}-t_i)\mathbf{m}_i + h_{01}(s)\mathbf{p}_{i+1} + h_{11}(s)(t_{i+1}-t_i)\mathbf{m}_{i+1}$$

其中 $\mathbf{m}_i$ 是点 $i$ 处的切向量（速度）。

**导数估计方法**：

1. **有限差分**：
   $$\mathbf{m}_i = \frac{\mathbf{p}_{i+1} - \mathbf{p}_{i-1}}{t_{i+1} - t_{i-1}}$$

2. **Catmull-Rom切线**：
   $$\mathbf{m}_i = \frac{1}{2}\left(\frac{\mathbf{p}_{i+1} - \mathbf{p}_i}{t_{i+1} - t_i} + \frac{\mathbf{p}_i - \mathbf{p}_{i-1}}{t_i - t_{i-1}}\right)$$

3. **Cardinal样条切线**：
   $$\mathbf{m}_i = (1-c)\frac{\mathbf{p}_{i+1} - \mathbf{p}_{i-1}}{t_{i+1} - t_{i-1}}$$
   
   其中 $c \in [0,1]$ 是张力参数。

#### 四元数样条

**Squad（球面四次插值）**：
$$\text{squad}(\mathbf{q}_i, \mathbf{q}_{i+1}, \mathbf{s}_i, \mathbf{s}_{i+1}, t) = \text{slerp}(\text{slerp}(\mathbf{q}_i, \mathbf{q}_{i+1}, t), \text{slerp}(\mathbf{s}_i, \mathbf{s}_{i+1}, t), 2t(1-t))$$

其中辅助四元数：
$$\mathbf{s}_i = \mathbf{q}_i \exp\left(-\frac{\log(\mathbf{q}_i^{-1}\mathbf{q}_{i-1}) + \log(\mathbf{q}_i^{-1}\mathbf{q}_{i+1})}{4}\right)$$

**Bezier四元数曲线**：
$$\mathbf{q}(t) = \text{slerp}(\text{slerp}(\text{slerp}(\mathbf{q}_0, \mathbf{q}_1, t), \text{slerp}(\mathbf{q}_1, \mathbf{q}_2, t), t), \text{slerp}(\text{slerp}(\mathbf{q}_1, \mathbf{q}_2, t), \text{slerp}(\mathbf{q}_2, \mathbf{q}_3, t), t), t)$$

#### 多维插值的耦合问题

当同时插值位置和旋转时，需要考虑它们的耦合关系：

**螺旋运动（Screw Motion）**：
结合平移和旋转的最自然方式：
$$\mathbf{T}(t) = \begin{bmatrix}
\mathbf{R}(t) & \mathbf{d}(t) \\
\mathbf{0} & 1
\end{bmatrix}$$

其中 $\mathbf{R}(t)$ 是旋转矩阵，$\mathbf{d}(t)$ 是位移向量。

**对偶四元数插值**：
使用对偶四元数 $\hat{\mathbf{q}} = \mathbf{q}_r + \epsilon\mathbf{q}_d$ 统一表示旋转和平移：
$$\hat{\mathbf{q}}(t) = \text{DLB}(\hat{\mathbf{q}}_0, \hat{\mathbf{q}}_1, t)$$

其中DLB是对偶四元数线性混合。

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

#### 非均匀参数化

标准Catmull-Rom假设均匀参数化，但实际关键帧时间可能不均匀。非均匀Catmull-Rom使用弦长或向心参数化：

**弦长参数化**：
$$t_{i+1} - t_i = \|\mathbf{p}_{i+1} - \mathbf{p}_i\|$$

**向心参数化**：
$$t_{i+1} - t_i = \|\mathbf{p}_{i+1} - \mathbf{p}_i\|^{0.5}$$

向心参数化避免了尖点和自交，产生更自然的曲线。

#### B样条在动画中的应用

B样条提供局部控制和更高的连续性：

**均匀三次B样条**：
$$\mathbf{p}(t) = \sum_{i} N_{i,3}(t)\mathbf{p}_i$$

基函数：
$$N_{i,3}(t) = \frac{1}{6}\begin{cases}
(1-t)^3 & t \in [0,1] \\
3t^3 - 12t^2 + 12t - 3 & t \in [1,2] \\
-3t^3 + 12t^2 - 12t + 1 & t \in [2,3] \\
t^3 & t \in [3,4]
\end{cases}$$

**NURBS（非均匀有理B样条）**：
$$\mathbf{p}(t) = \frac{\sum_{i} w_i N_{i,k}(t)\mathbf{p}_i}{\sum_{i} w_i N_{i,k}(t)}$$

权重$w_i$提供额外的形状控制，特别适合表示圆锥曲线。

#### 样条的连续性分析

**参数连续性**：
- $C^0$：位置连续
- $C^1$：切向量连续
- $C^2$：曲率连续

**几何连续性**：
- $G^0$：位置连续
- $G^1$：切向连续（方向相同，大小可不同）
- $G^2$：曲率连续

对于动画，$G^1$连续通常足够，但相机路径可能需要$G^2$以避免视觉跳动。

#### 张力-连续性-偏置（TCB）样条

Kochanek-Bartels样条提供三个参数控制：

$$\mathbf{m}_i^{in} = \frac{(1-t)(1+c)(1+b)}{2}(\mathbf{p}_i - \mathbf{p}_{i-1}) + \frac{(1-t)(1-c)(1-b)}{2}(\mathbf{p}_{i+1} - \mathbf{p}_i)$$

$$\mathbf{m}_i^{out} = \frac{(1-t)(1-c)(1+b)}{2}(\mathbf{p}_i - \mathbf{p}_{i-1}) + \frac{(1-t)(1+c)(1-b)}{2}(\mathbf{p}_{i+1} - \mathbf{p}_i)$$

其中：
- $t$：张力（tension），控制曲线的"紧度"
- $c$：连续性（continuity），控制拐角的尖锐度
- $b$：偏置（bias），控制曲线偏向前后控制点的程度

#### 样条曲线的实时评估优化

**前向差分法**：
对于均匀参数的三次曲线，使用前向差分避免重复计算：

$$\Delta^0 \mathbf{p} = \mathbf{p}(0)$$
$$\Delta^1 \mathbf{p} = a\mathbf{v}_0 + b\mathbf{v}_1 + c\mathbf{p}_0 + d\mathbf{p}_1$$
$$\Delta^2 \mathbf{p} = 6ah^2\mathbf{v}_0 + 6bh^2\mathbf{v}_1$$
$$\Delta^3 \mathbf{p} = 6h^3(\mathbf{v}_0 + \mathbf{v}_1)$$

递推关系：
$$\mathbf{p}_{i+1} = \mathbf{p}_i + \Delta^1\mathbf{p}_i$$
$$\Delta^1\mathbf{p}_{i+1} = \Delta^1\mathbf{p}_i + \Delta^2\mathbf{p}_i$$
$$\Delta^2\mathbf{p}_{i+1} = \Delta^2\mathbf{p}_i + \Delta^3\mathbf{p}_i$$

**de Casteljau算法的并行化**：
贝塞尔曲线的de Casteljau算法天然适合SIMD并行：
$$\mathbf{b}_i^{(k)} = (1-t)\mathbf{b}_i^{(k-1)} + t\mathbf{b}_{i+1}^{(k-1)}$$

可以同时计算多个$t$值或多条曲线。

### 11.1.4 动画压缩与优化

实际应用中，动画数据可能非常庞大。考虑一个30fps的动画，每秒需要存储30帧数据，一分钟的动画就需要1800帧。对于复杂的角色动画，每帧可能包含数百个骨骼的变换矩阵。

#### 存储需求分析

考虑一个典型的人形角色：
- 骨骼数量：$n = 60$
- 每骨骼存储：位置(3) + 四元数(4) + 缩放(3) = 10个浮点数
- 每帧数据：$60 \times 10 \times 4 = 2400$ 字节
- 30fps动画每秒：$2400 \times 30 = 72$ KB
- 一分钟动画：$72 \times 60 = 4.32$ MB

这还不包括面部动画、手指细节等。

#### 关键帧约简

**Douglas-Peucker算法**用于移除冗余关键帧：

1. 连接首尾关键帧形成直线
2. 找到离直线最远的中间帧
3. 如果距离超过阈值$\epsilon$，保留该帧并递归处理两段
4. 否则移除所有中间帧

时间复杂度：$O(n\log n)$（平均情况）

**自适应误差度量**：
不同属性使用不同阈值：
- 位置：$\epsilon_{pos} = 0.001$ 单位
- 旋转：$\epsilon_{rot} = 0.5°$
- 缩放：$\epsilon_{scale} = 0.01$

**感知重要性加权**：
$$d_{weighted} = w_{bone} \cdot w_{velocity} \cdot d_{geometric}$$

其中：
- $w_{bone}$：骨骼重要性（脊椎 > 四肢末端）
- $w_{velocity}$：速度权重（快速运动允许更大误差）

#### 曲线拟合压缩

使用**最小二乘法**拟合动画曲线：

对于$n$次多项式$p(t) = \sum_{i=0}^n a_i t^i$，求解：
$$\min_{\{a_i\}} \sum_{j} \|\mathbf{p}(t_j) - \mathbf{p}_{data}(t_j)\|^2$$

这导致法方程：
$$\mathbf{A}^T\mathbf{A}\mathbf{x} = \mathbf{A}^T\mathbf{b}$$

其中$\mathbf{A}_{ji} = t_j^i$。

**分段拟合策略**：
1. 自动检测断点（速度/加速度突变）
2. 每段独立拟合
3. 确保段间连续性

**贝塞尔曲线压缩**：
使用三次贝塞尔曲线，仅存储4个控制点：
$$\mathbf{p}(t) = (1-t)^3\mathbf{P}_0 + 3(1-t)^2t\mathbf{P}_1 + 3(1-t)t^2\mathbf{P}_2 + t^3\mathbf{P}_3$$

控制点优化：
$$\min_{\{\mathbf{P}_i\}} \sum_{j} \|\mathbf{p}(t_j) - \mathbf{p}_{data}(t_j)\|^2$$

受约束于：$\mathbf{P}_0 = \mathbf{p}_{data}(0)$，$\mathbf{P}_3 = \mathbf{p}_{data}(1)$

#### 主成分分析（PCA）压缩

对于$m$个相关的动画通道（如多个角色执行相似动作）：

1. 构建数据矩阵$\mathbf{X} \in \mathbb{R}^{m \times n}$
2. 计算协方差矩阵$\mathbf{C} = \frac{1}{n-1}\mathbf{X}\mathbf{X}^T$
3. 特征值分解$\mathbf{C} = \mathbf{U}\mathbf{\Lambda}\mathbf{U}^T$
4. 保留前$k$个主成分：$\mathbf{X}_{compressed} = \mathbf{U}_k^T\mathbf{X}$

压缩率：$\frac{k(m+n)}{mn}$

**增量PCA更新**：
新增动画时无需重新计算整个PCA：
$$\mathbf{C}_{new} = \frac{n-1}{n}\mathbf{C}_{old} + \frac{1}{n}\mathbf{x}_{new}\mathbf{x}_{new}^T$$

**局部PCA分解**：
1. 时间分段：将长动画分成短片段
2. 空间分组：将骨骼按运动相关性分组
3. 每组独立PCA

#### 小波压缩

使用多分辨率分析：
$$\mathbf{x}(t) = \sum_{k} c_{j_0,k}\phi_{j_0,k}(t) + \sum_{j=j_0}^{J} \sum_{k} d_{j,k}\psi_{j,k}(t)$$

其中：
- $\phi$：尺度函数
- $\psi$：小波函数
- $c_{j,k}$：尺度系数
- $d_{j,k}$：小波系数

**阈值处理**：
$$\tilde{d}_{j,k} = \begin{cases}
d_{j,k} & |d_{j,k}| > \epsilon_j \\
0 & \text{otherwise}
\end{cases}$$

自适应阈值：$\epsilon_j = \sigma\sqrt{2\log n}/2^j$

#### 四元数压缩

**最小三参数表示**：
利用单位四元数约束$\|\mathbf{q}\| = 1$：
1. 找到最大分量索引
2. 存储其他三个分量
3. 重建时：$q_{max} = \sqrt{1 - q_x^2 - q_y^2 - q_z^2}$

**对数映射压缩**：
$$\mathbf{v} = \log(\mathbf{q}) = \begin{cases}
\mathbf{0} & \theta = 0 \\
\frac{\theta}{\sin\theta}[q_x, q_y, q_z]^T & \text{otherwise}
\end{cases}$$

优势：线性空间，适合PCA

#### 误差度量与质量控制

除了简单的L2误差，还应考虑：

**感知误差**：
$$E_{perceptual} = \sum_t w(t)\|\mathbf{p}_{original}(t) - \mathbf{p}_{compressed}(t)\|^2$$

其中$w(t)$是基于运动速度的权重函数：
$$w(t) = 1 + \alpha\|\mathbf{v}(t)\|$$

**关节角度误差**（对于骨骼动画）：
$$E_{angle} = \sum_{joint} \sum_t \|\boldsymbol{\theta}_{original}(t) - \boldsymbol{\theta}_{compressed}(t)\|^2$$

**端点误差**（IK链）：
$$E_{endpoint} = \sum_t \|FK(\boldsymbol{\theta}_{original}(t)) - FK(\boldsymbol{\theta}_{compressed}(t))\|^2$$

#### 运行时解压优化

**层次细节（LOD）**：
- 远距离：仅主要骨骼
- 中距离：添加次要骨骼
- 近距离：完整骨骼+面部

**预测解码**：
$$\mathbf{x}(t) = \mathbf{x}_{base}(t) + \sum_{i=1}^{k} \alpha_i(t)\mathbf{e}_i$$

其中$\mathbf{e}_i$是预计算的误差模式。

**GPU友好格式**：
- 纹理存储：动画数据作为纹理
- 顶点纹理获取（VTF）
- 计算着色器解压

### 11.1.5 程序化动画与噪声函数

程序化动画通过数学函数生成动画，无需存储关键帧数据。

#### Perlin噪声

用于生成自然的随机运动：
$$noise(\mathbf{x}) = \sum_{i} fade(f_i) \cdot lerp(w_i, grad_i \cdot \mathbf{x}_i, grad_{i+1} \cdot \mathbf{x}_{i+1})$$

其中$fade(t) = 6t^5 - 15t^4 + 10t^3$保证C2连续性。

**梯度生成**：
使用伪随机数生成器，保证可重复性：
$$grad(\mathbf{i}) = normalize(hash(\mathbf{i}) \mod \mathbf{g}[hash(\mathbf{i}) \mod |\mathbf{g}|])$$

其中$\mathbf{g}$是预定义的梯度表。

**改进Perlin噪声**：
$$fade_{improved}(t) = t^3(t(t \cdot 6 - 15) + 10)$$

这保证了$fade'(0) = fade'(1) = 0$和$fade''(0) = fade''(1) = 0$。

#### Simplex噪声

Perlin的改进版本，计算效率更高：

**2D Simplex**：
$$F = \frac{\sqrt{3} - 1}{2}, \quad G = \frac{3 - \sqrt{3}}{6}$$

坐标变换：
$$s = (x + y) \cdot F$$
$$i = \lfloor x + s \rfloor, \quad j = \lfloor y + s \rfloor$$

**优势**：
- 更少的计算量（$2^n$顶点 vs $n!$顶点）
- 更好的各向同性
- 无方向伪影

#### 分形布朗运动（fBm）

通过叠加不同频率的噪声：
$$fBm(\mathbf{x}) = \sum_{i=0}^{n} \frac{noise(2^i \mathbf{x})}{2^{iH}}$$

其中$H$是Hurst指数，控制粗糙度。

**变体**：
- **Turbulence**：$turbulence(\mathbf{x}) = \sum_{i=0}^{n} \frac{|noise(2^i \mathbf{x})|}{2^i}$
- **Ridged noise**：$ridge(\mathbf{x}) = \sum_{i=0}^{n} \frac{(1-|noise(2^i \mathbf{x})|)^2}{2^i}$

#### 谐波叠加

用于周期性运动（如呼吸、摆动）：
$$\mathbf{p}(t) = \mathbf{p}_0 + \sum_{i=1}^n A_i \sin(\omega_i t + \phi_i)\mathbf{d}_i$$

通过调整振幅$A_i$、频率$\omega_i$和相位$\phi_i$可以创建复杂的周期运动。

**从参考运动提取参数**：
使用FFT分析：
$$A_k = \frac{2}{N}\left|\sum_{n=0}^{N-1} \mathbf{p}(t_n)e^{-2\pi ikn/N}\right|$$
$$\phi_k = \arg\left(\sum_{n=0}^{N-1} \mathbf{p}(t_n)e^{-2\pi ikn/N}\right)$$

#### Gabor噪声

结合频域和空域特性：
$$gabor(\mathbf{x}) = \sum_{i=1}^{n} w_i \cdot g(\mathbf{x} - \mathbf{x}_i) \cdot \cos(2\pi\mathbf{k}_i \cdot (\mathbf{x} - \mathbf{x}_i) + \phi_i)$$

其中$g$是高斯核：
$$g(\mathbf{x}) = K \exp\left(-\pi\|\mathbf{x}\|^2/a^2\right)$$

**参数控制**：
- $\mathbf{k}_i$：频率向量
- $a$：带宽
- $\phi_i$：相位

#### 流体噪声

时间相关的噪声，保持时间连续性：
$$flow(\mathbf{x}, t) = \sum_{i} \alpha_i(t) \cdot noise(\mathbf{x} + \mathbf{v}_i t)$$

其中：
$$\alpha_i(t) = \max(0, 1 - |t - t_i|/T)$$

这确保了平滑的时间过渡。

#### 向量场动画

使用向量场驱动粒子运动：
$$\frac{d\mathbf{x}}{dt} = \mathbf{v}(\mathbf{x}, t)$$

**无散度向量场**：
$$\mathbf{v} = \nabla \times \boldsymbol{\psi}$$

其中$\boldsymbol{\psi}$是流函数。

**流场的程序化生成**：
- 温流：$\mathbf{v} = \nabla \times (noise(\mathbf{x})\hat{\mathbf{z}})$
- 漩流：$\mathbf{v} = \omega \times (\mathbf{x} - \mathbf{c})$
- 汇/源：$\mathbf{v} = \pm\frac{\mathbf{x} - \mathbf{c}}{\|\mathbf{x} - \mathbf{c}\|^2}$

#### 群体行为模拟

**Boids模型**：
1. **分离**：$\mathbf{f}_{sep} = -\sum_{j \in neighbors} \frac{\mathbf{x}_j - \mathbf{x}_i}{\|\mathbf{x}_j - \mathbf{x}_i\|^2}$
2. **对齐**：$\mathbf{f}_{align} = \frac{1}{|neighbors|}\sum_{j \in neighbors} \mathbf{v}_j - \mathbf{v}_i$
3. **聚集**：$\mathbf{f}_{cohesion} = \frac{1}{|neighbors|}\sum_{j \in neighbors} \mathbf{x}_j - \mathbf{x}_i$

**势场方法**：
$$\mathbf{a}_i = -\nabla U(\mathbf{x}_i) + \sum_j \mathbf{f}_{interaction}(\mathbf{x}_i, \mathbf{x}_j)$$

其中$U$是环境势场。

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

#### 雅可比矩阵的几何意义

雅可比矩阵描述了关节空间速度到笛卡尔空间速度的线性映射：
$$\dot{\mathbf{x}} = \mathbf{J}(\boldsymbol{\theta})\dot{\boldsymbol{\theta}}$$

对于旋转关节$i$，其对末端执行器的贡献：
$$\mathbf{J}_i = \begin{cases}
[\mathbf{z}_i \times (\mathbf{p}_{end} - \mathbf{p}_i)]^T & \text{位置部分} \\
\mathbf{z}_i^T & \text{方向部分}
\end{cases}$$

其中$\mathbf{z}_i$是关节$i$的旋转轴，$\mathbf{p}_i$是关节位置。

#### 奇异性分析

**奇异配置**发生在：
1. **边界奇异性**：机械臂完全伸展或收缩
2. **内部奇异性**：两个或多个关节轴对齐

奇异性的数学特征：
- $\det(\mathbf{J}\mathbf{J}^T) = 0$
- 雅可比矩阵降秩
- 某些方向上失去自由度

**可操作性椭球**：
$$\mathcal{E} = \{\mathbf{v} : \mathbf{v}^T(\mathbf{J}\mathbf{J}^T)^{-1}\mathbf{v} \leq 1\}$$

椭球的主轴由$\mathbf{J}\mathbf{J}^T$的特征向量给出，轴长为特征值的平方根。

#### 处理奇异性的高级方法

1. **阻尼最小二乘（DLS）**：
$$\Delta\boldsymbol{\theta} = \mathbf{J}^T(\mathbf{J}\mathbf{J}^T + \lambda^2\mathbf{I})^{-1}\Delta\mathbf{x}$$

自适应阻尼因子：
$$\lambda = \begin{cases}
0 & \text{if } \sigma_{min} \geq \epsilon \\
\sqrt{\epsilon^2 - \sigma_{min}^2} & \text{otherwise}
\end{cases}$$

2. **奇异值分解（SVD）过滤**：
$$\mathbf{J} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

过滤伪逆：
$$\mathbf{J}^+ = \sum_{i=1}^r \frac{1}{\sigma_i} \mathbf{v}_i \mathbf{u}_i^T, \quad \text{仅当} \sigma_i > \epsilon$$

3. **梯度投影法**：

在零空间中优化次要目标：
$$\Delta\boldsymbol{\theta} = \mathbf{J}^+\Delta\mathbf{x} + (\mathbf{I} - \mathbf{J}^+\mathbf{J})\nabla h(\boldsymbol{\theta})$$

其中$h(\boldsymbol{\theta})$是次要目标函数（如避免关节限位）。

#### 任务优先级方法

处理多个任务时，使用优先级：

主任务：$\mathbf{J}_1\dot{\boldsymbol{\theta}} = \dot{\mathbf{x}}_1$

次任务在主任务零空间中：
$$\dot{\boldsymbol{\theta}} = \mathbf{J}_1^+\dot{\mathbf{x}}_1 + (\mathbf{I} - \mathbf{J}_1^+\mathbf{J}_1)\mathbf{J}_2^+\dot{\mathbf{x}}_2$$

### 11.2.6 高级约束处理

#### 位置基约束（PBD）

直接在位置级别满足约束：

1. 预测位置：$\mathbf{p}_i^* = \mathbf{x}_i + \Delta t \mathbf{v}_i$
2. 求解约束：
   $$\Delta\mathbf{p}_i = -\frac{C(\mathbf{p}^*)}{\sum_j w_j|\nabla_{\mathbf{p}_j}C|^2}w_i\nabla_{\mathbf{p}_i}C$$
3. 更新位置：$\mathbf{p}_i = \mathbf{p}_i^* + \Delta\mathbf{p}_i$
4. 更新速度：$\mathbf{v}_i = (\mathbf{p}_i - \mathbf{x}_i)/\Delta t$

其中$w_i = 1/m_i$是逆质量。

#### XPBD（扩展PBD）

引入拉格朗日乘数使约束更物理：
$$\Delta\lambda = \frac{-C - \tilde{\alpha}\lambda}{\sum_j w_j|\nabla_{\mathbf{p}_j}C|^2 + \tilde{\alpha}}$$

其中$\tilde{\alpha} = \alpha/(\Delta t)^2$是时间步长相关的柔度。

#### 分层求解器

对于复杂系统，使用多级求解：

1. **粗粒度**：解决全局约束（如整体形状保持）
2. **中粒度**：解决局部约束（如局部碰撞）  
3. **细粒度**：解决细节约束（如褶皱细节）

每级使用不同的迭代次数和精度要求。

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

#### 哈密顿系统与辛结构

哈密顿系统的相空间演化保持辛结构：
$$\frac{d}{dt}\begin{bmatrix} \mathbf{q} \\ \mathbf{p} \end{bmatrix} = \begin{bmatrix} \mathbf{0} & \mathbf{I} \\ -\mathbf{I} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \nabla_\mathbf{q} H \\ \nabla_\mathbf{p} H \end{bmatrix}$$

其中$H(\mathbf{q}, \mathbf{p})$是哈密顿量。

**Liouville定理**：相空间体积在哈密顿流下保持不变。

#### 辛积分器

**Störmer-Verlet方法**：
$$\mathbf{x}_{n+1} = 2\mathbf{x}_n - \mathbf{x}_{n-1} + h^2\mathbf{M}^{-1}\mathbf{f}_n$$

这是二阶精度的辛积分器。速度计算：
$$\mathbf{v}_n = \frac{\mathbf{x}_{n+1} - \mathbf{x}_{n-1}}{2h} + O(h^2)$$

**速度Verlet方法**（数值上更稳定）：
$$\mathbf{x}_{n+1} = \mathbf{x}_n + h\mathbf{v}_n + \frac{h^2}{2}\mathbf{a}_n$$
$$\mathbf{v}_{n+1} = \mathbf{v}_n + \frac{h}{2}(\mathbf{a}_n + \mathbf{a}_{n+1})$$

#### 能量误差分析

对于辛积分器，能量误差是有界的：
$$|H(t) - H(0)| \leq Ch^p$$

其中$p$是方法的阶数，$C$与时间$t$无关。

非辛方法（如RK4）的能量误差可能线性增长：
$$|H(t) - H(0)| \sim Cth^p$$

#### 后向误差分析

辛积分器精确求解一个扰动的哈密顿系统：
$$\tilde{H} = H + h^2H_2 + h^4H_4 + \cdots$$

这解释了为什么辛方法长期行为更好。

#### 约束系统的能量守恒

对于约束系统，使用**SHAKE算法**：

1. 无约束更新：$\mathbf{r}^* = \mathbf{r}_n + h\mathbf{v}_n + \frac{h^2}{2m}\mathbf{f}_n$
2. 求解约束：$g(\mathbf{r}_{n+1}) = 0$
3. 约束力：$\mathbf{f}_c = \frac{m}{h^2}(\mathbf{r}_{n+1} - \mathbf{r}^*)$

**RATTLE算法**同时约束位置和速度。

### 11.3.6 高阶积分方法与分裂算法

#### 高阶辛积分器

**Forest-Ruth方法**（4阶）：
$$\begin{align}
\mathbf{x}_{1/2} &= \mathbf{x}_n + \theta h \mathbf{v}_n/2 \\
\mathbf{v}_1 &= \mathbf{v}_n + \theta h \mathbf{a}(\mathbf{x}_{1/2}) \\
\mathbf{x}_1 &= \mathbf{x}_{1/2} + (1-\theta)h\mathbf{v}_1/2 \\
\mathbf{v}_{3/2} &= \mathbf{v}_1 + (1-2\theta)h\mathbf{a}(\mathbf{x}_1) \\
\mathbf{x}_{3/2} &= \mathbf{x}_1 + (1-\theta)h\mathbf{v}_{3/2}/2 \\
\mathbf{v}_2 &= \mathbf{v}_{3/2} + \theta h \mathbf{a}(\mathbf{x}_{3/2}) \\
\mathbf{x}_{n+1} &= \mathbf{x}_{3/2} + \theta h \mathbf{v}_2/2
\end{align}$$

其中$\theta = 2^{1/3}/(2-2^{1/3})$。

#### 算子分裂方法

对于$\frac{d\mathbf{y}}{dt} = \mathbf{A}(\mathbf{y}) + \mathbf{B}(\mathbf{y})$：

**Strang分裂**（二阶）：
$$\mathbf{y}_{n+1} = e^{h\mathbf{B}/2} e^{h\mathbf{A}} e^{h\mathbf{B}/2} \mathbf{y}_n$$

**Yoshida分裂**（4阶）：
$$\mathbf{y}_{n+1} = e^{w_1h\mathbf{A}} e^{w_1h\mathbf{B}} e^{w_2h\mathbf{A}} e^{w_2h\mathbf{B}} e^{w_3h\mathbf{A}} e^{w_3h\mathbf{B}} e^{w_2h\mathbf{A}} e^{w_2h\mathbf{B}} e^{w_1h\mathbf{A}} e^{w_1h\mathbf{B}} \mathbf{y}_n$$

其中权重经过精心选择以消除低阶误差项。

#### 多时间尺度方法

对于刚性系统，使用**IMEX方法**：
$$\mathbf{y}_{n+1} = \mathbf{y}_n + h[\mathbf{f}_{explicit}(\mathbf{y}_n) + \mathbf{f}_{implicit}(\mathbf{y}_{n+1})]$$

将快变（刚性）部分隐式处理，慢变部分显式处理。

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