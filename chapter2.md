# 第2章：几何变换

几何变换是计算机图形学的核心基础之一，它决定了我们如何在二维和三维空间中移动、旋转、缩放物体，以及如何将三维世界投影到二维屏幕上。本章将深入探讨变换的数学原理、实现细节和优化技巧，为后续的渲染管线学习打下坚实基础。

## 2.1 二维与三维变换

### 2.1.1 齐次坐标系

在计算机图形学中，我们使用齐次坐标（Homogeneous Coordinates）来统一表示点和向量，并使得仿射变换可以用矩阵乘法来表示。

对于二维空间：
- 点：$(x, y) \rightarrow (x, y, 1)$
- 向量：$(x, y) \rightarrow (x, y, 0)$

对于三维空间：
- 点：$(x, y, z) \rightarrow (x, y, z, 1)$
- 向量：$(x, y, z) \rightarrow (x, y, z, 0)$

齐次坐标的关键性质：
- $(x, y, z, w)$ 和 $(kx, ky, kz, kw)$ 表示同一个点（$k \neq 0$）
- 当 $w \neq 0$ 时，齐次坐标 $(x, y, z, w)$ 对应的笛卡尔坐标为 $(x/w, y/w, z/w)$

### 2.1.2 基本二维变换

**平移变换（Translation）**

$$\mathbf{T}(t_x, t_y) = \begin{bmatrix}
1 & 0 & t_x \\
0 & 1 & t_y \\
0 & 0 & 1
\end{bmatrix}$$

**旋转变换（Rotation）**

绕原点逆时针旋转 $\theta$ 角度：

$$\mathbf{R}(\theta) = \begin{bmatrix}
\cos\theta & -\sin\theta & 0 \\
\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{bmatrix}$$

**缩放变换（Scaling）**

$$\mathbf{S}(s_x, s_y) = \begin{bmatrix}
s_x & 0 & 0 \\
0 & s_y & 0 \\
0 & 0 & 1
\end{bmatrix}$$

**错切变换（Shearing）**

$$\mathbf{H}_x(s) = \begin{bmatrix}
1 & s & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}, \quad
\mathbf{H}_y(s) = \begin{bmatrix}
1 & 0 & 0 \\
s & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

### 2.1.3 基本三维变换

**三维平移**

$$\mathbf{T}(t_x, t_y, t_z) = \begin{bmatrix}
1 & 0 & 0 & t_x \\
0 & 1 & 0 & t_y \\
0 & 0 & 1 & t_z \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**三维缩放**

$$\mathbf{S}(s_x, s_y, s_z) = \begin{bmatrix}
s_x & 0 & 0 & 0 \\
0 & s_y & 0 & 0 \\
0 & 0 & s_z & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**三维旋转**

绕 x 轴旋转：
$$\mathbf{R}_x(\theta) = \begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & \cos\theta & -\sin\theta & 0 \\
0 & \sin\theta & \cos\theta & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

绕 y 轴旋转：
$$\mathbf{R}_y(\theta) = \begin{bmatrix}
\cos\theta & 0 & \sin\theta & 0 \\
0 & 1 & 0 & 0 \\
-\sin\theta & 0 & \cos\theta & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

绕 z 轴旋转：
$$\mathbf{R}_z(\theta) = \begin{bmatrix}
\cos\theta & -\sin\theta & 0 & 0 \\
\sin\theta & \cos\theta & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

### 2.1.4 任意轴旋转（Rodrigues' Rotation Formula）

绕任意单位向量 $\mathbf{n} = (n_x, n_y, n_z)$ 旋转 $\theta$ 角度：

$$\mathbf{R}(\mathbf{n}, \theta) = \cos\theta \mathbf{I} + (1-\cos\theta)\mathbf{n}\mathbf{n}^T + \sin\theta \mathbf{N}$$

其中 $\mathbf{N}$ 是 $\mathbf{n}$ 的反对称矩阵：

$$\mathbf{N} = \begin{bmatrix}
0 & -n_z & n_y \\
n_z & 0 & -n_x \\
-n_y & n_x & 0
\end{bmatrix}$$

### 2.1.5 变换的分解与组合

任何仿射变换都可以分解为：平移 × 旋转 × 缩放 × 错切

对于刚体变换（保持形状和大小）：
$$\mathbf{M} = \mathbf{T} \cdot \mathbf{R}$$

变换的组合遵循矩阵乘法规则，注意顺序很重要：
$$\mathbf{M}_{\text{combined}} = \mathbf{M}_n \cdot \mathbf{M}_{n-1} \cdot ... \cdot \mathbf{M}_2 \cdot \mathbf{M}_1$$

## 2.2 模型、视图、投影变换

### 2.2.1 坐标系统与变换管线

图形渲染管线中的坐标系统：
1. **模型空间（Model Space）**：物体的局部坐标系
2. **世界空间（World Space）**：场景的全局坐标系
3. **视图空间（View/Camera Space）**：以相机为原点的坐标系
4. **裁剪空间（Clip Space）**：投影后的齐次坐标系
5. **NDC空间（Normalized Device Coordinates）**：归一化的设备坐标
6. **屏幕空间（Screen Space）**：最终的像素坐标

### 2.2.2 模型变换（Model Transform）

模型变换将物体从模型空间变换到世界空间：
$$\mathbf{M}_{\text{model}} = \mathbf{T}_{\text{world}} \cdot \mathbf{R}_{\text{world}} \cdot \mathbf{S}_{\text{local}}$$

实践中常用的优化：
- 预计算静态物体的模型矩阵
- 使用实例化（Instancing）减少重复计算
- 层次化场景图进行变换继承

### 2.2.3 视图变换（View Transform）

视图变换将世界空间变换到视图空间。给定相机参数：
- 位置：$\mathbf{e}$ (eye position)
- 观察方向：$\mathbf{g}$ (gaze direction)
- 上方向：$\mathbf{t}$ (up direction)

构建相机坐标系：
$$\mathbf{w} = -\frac{\mathbf{g}}{||\mathbf{g}||}, \quad
\mathbf{u} = \frac{\mathbf{t} \times \mathbf{w}}{||\mathbf{t} \times \mathbf{w}||}, \quad
\mathbf{v} = \mathbf{w} \times \mathbf{u}$$

视图矩阵：
$$\mathbf{V} = \begin{bmatrix}
\mathbf{u}^T & -\mathbf{u} \cdot \mathbf{e} \\
\mathbf{v}^T & -\mathbf{v} \cdot \mathbf{e} \\
\mathbf{w}^T & -\mathbf{w} \cdot \mathbf{e} \\
0 & 1
\end{bmatrix}$$

### 2.2.4 投影变换（Projection Transform）

**透视投影（Perspective Projection）**

给定视锥体参数：
- 垂直视场角：$\text{fov}_y$
- 宽高比：$\text{aspect}$
- 近平面：$n$
- 远平面：$f$

透视投影矩阵：
$$\mathbf{P}_{\text{persp}} = \begin{bmatrix}
\frac{1}{\text{aspect} \cdot \tan(\text{fov}_y/2)} & 0 & 0 & 0 \\
0 & \frac{1}{\tan(\text{fov}_y/2)} & 0 & 0 \\
0 & 0 & \frac{f+n}{n-f} & \frac{2fn}{n-f} \\
0 & 0 & -1 & 0
\end{bmatrix}$$

**正交投影（Orthographic Projection）**

给定包围盒 $[l,r] \times [b,t] \times [n,f]$：

$$\mathbf{P}_{\text{ortho}} = \begin{bmatrix}
\frac{2}{r-l} & 0 & 0 & -\frac{r+l}{r-l} \\
0 & \frac{2}{t-b} & 0 & -\frac{t+b}{t-b} \\
0 & 0 & \frac{2}{n-f} & -\frac{n+f}{n-f} \\
0 & 0 & 0 & 1
\end{bmatrix}$$

### 2.2.5 视口变换（Viewport Transform）

将NDC坐标 $[-1,1]^3$ 映射到屏幕坐标 $[0,w] \times [0,h]$：

$$\mathbf{M}_{\text{viewport}} = \begin{bmatrix}
\frac{w}{2} & 0 & 0 & \frac{w}{2} \\
0 & \frac{h}{2} & 0 & \frac{h}{2} \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

## 2.3 变换矩阵的组合与优化

### 2.3.1 矩阵组合的数学性质

**结合律**：$(AB)C = A(BC)$
- 允许预计算变换序列
- 减少实时计算量

**不满足交换律**：$AB \neq BA$（一般情况）
- 变换顺序至关重要
- 常见错误：先旋转后平移 vs 先平移后旋转

**逆变换**：
- 平移的逆：$\mathbf{T}^{-1}(t_x, t_y, t_z) = \mathbf{T}(-t_x, -t_y, -t_z)$
- 旋转的逆：$\mathbf{R}^{-1} = \mathbf{R}^T$（正交矩阵性质）
- 缩放的逆：$\mathbf{S}^{-1}(s_x, s_y, s_z) = \mathbf{S}(1/s_x, 1/s_y, 1/s_z)$
- 组合的逆：$(AB)^{-1} = B^{-1}A^{-1}$

### 2.3.2 变换矩阵的优化技巧

**1. 矩阵预计算**
```
MVP = Projection × View × Model  // 预计算，每个物体只算一次
```

**2. 矩阵分解存储**
- 存储位置、旋转（四元数）、缩放分量
- 仅在需要时构建完整矩阵
- 节省内存，提高缓存效率

**3. 特殊矩阵的快速计算**
- 利用正交矩阵性质：$\mathbf{R}^{-1} = \mathbf{R}^T$
- 对角矩阵乘法优化
- 稀疏矩阵表示

**4. SIMD优化**
- 使用 SSE/AVX 指令集
- 矩阵-向量乘法的并行化
- 批量变换处理

### 2.3.3 法向量变换

法向量的变换不能直接使用模型矩阵，需要使用逆转置矩阵：

设模型变换矩阵为 $\mathbf{M}$，法向量变换矩阵为：
$$\mathbf{G} = (\mathbf{M}^{-1})^T$$

推导：保持法向量与切向量垂直
- 原始：$\mathbf{n} \cdot \mathbf{t} = 0$
- 变换后：$\mathbf{n}' \cdot \mathbf{t}' = 0$
- 其中 $\mathbf{t}' = \mathbf{M}\mathbf{t}$

对于仅包含旋转的变换，$\mathbf{G} = \mathbf{M}$（因为 $\mathbf{R}^{-1} = \mathbf{R}^T$）

### 2.3.4 精度问题与数值稳定性

**浮点精度损失**
- 连续变换累积误差
- 远离原点的大坐标值
- 解决：定期正交化、双精度关键计算

**深度精度（Z-fighting）**
- 近远平面比例过大
- 解决：对数深度缓冲、反向Z

**矩阵正交化**
使用 Gram-Schmidt 过程保持旋转矩阵的正交性：
```
重新正交化旋转矩阵的前三列
```

### 2.3.5 变换缓存与层次结构

**场景图优化**
- 父子关系的变换继承
- 脏标记（Dirty Flag）机制
- 惰性求值（Lazy Evaluation）

**实例化渲染**
- 基础模型矩阵 + 实例偏移
- GPU实例化缓冲区
- 减少 Draw Call

## 本章小结

### 核心概念
1. **齐次坐标**：统一点和向量表示，使仿射变换可用矩阵乘法表示
2. **基本变换**：平移、旋转、缩放、错切
3. **变换管线**：模型→世界→视图→裁剪→NDC→屏幕
4. **投影类型**：透视投影（真实感）vs 正交投影（工程图）
5. **优化技巧**：预计算、SIMD、实例化

### 关键公式
- Rodrigues旋转公式：$\mathbf{R}(\mathbf{n}, \theta) = \cos\theta \mathbf{I} + (1-\cos\theta)\mathbf{n}\mathbf{n}^T + \sin\theta \mathbf{N}$
- 视图矩阵：基于相机位置和朝向构建
- 透视投影：$w' = -z$，产生透视除法
- 法向量变换：$\mathbf{G} = (\mathbf{M}^{-1})^T$

### 性能考虑
- 矩阵乘法顺序影响计算量
- 预计算静态变换
- 利用特殊矩阵性质（正交、对称等）
- 批量处理和GPU并行化

## 练习题

### 基础题

**练习 2.1** 证明二维旋转矩阵 $\mathbf{R}(\theta)$ 是正交矩阵，即 $\mathbf{R}^T\mathbf{R} = \mathbf{I}$。

*Hint*: 计算 $\mathbf{R}^T(\theta)\mathbf{R}(\theta)$ 并利用三角恒等式。

<details>
<summary>答案</summary>

$$\mathbf{R}^T(\theta)\mathbf{R}(\theta) = \begin{bmatrix}
\cos\theta & \sin\theta \\
-\sin\theta & \cos\theta
\end{bmatrix}
\begin{bmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{bmatrix}$$

$$= \begin{bmatrix}
\cos^2\theta + \sin^2\theta & 0 \\
0 & \cos^2\theta + \sin^2\theta
\end{bmatrix} = \begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix} = \mathbf{I}$$

</details>

**练习 2.2** 给定点 $P = (2, 3, 1)$，先绕原点旋转 $45°$，再平移 $(1, -1, 0)$。写出组合变换矩阵并计算变换后的点坐标。

*Hint*: 注意变换顺序，先旋转后平移意味着 $\mathbf{T} \cdot \mathbf{R}$。

<details>
<summary>答案</summary>

组合变换矩阵：
$$\mathbf{M} = \mathbf{T}(1,-1,0) \cdot \mathbf{R}_z(45°) = \begin{bmatrix}
\frac{\sqrt{2}}{2} & -\frac{\sqrt{2}}{2} & 0 & 1 \\
\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2} & 0 & -1 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

变换后的点：
$$P' = \mathbf{M} \cdot [2, 3, 1, 1]^T = [-\frac{\sqrt{2}}{2} + 1, \frac{5\sqrt{2}}{2} - 1, 1, 1]^T$$

</details>

**练习 2.3** 推导透视投影后，为什么需要进行透视除法（除以 $w$ 分量）？

*Hint*: 考虑相似三角形原理和齐次坐标的意义。

<details>
<summary>答案</summary>

透视投影模拟了小孔成像原理。对于点 $(x, y, z)$，投影到近平面 $z = n$ 上的坐标应该是：
$$x' = \frac{nx}{z}, \quad y' = \frac{ny}{z}$$

透视投影矩阵将 $z$ 存储在 $w$ 分量中（通常是 $-z$），所以齐次坐标变为 $(x", y", z", -z)$。
透视除法 $(x"/w, y"/w, z"/w)$ 正好实现了上述的除以 $z$ 操作，完成了透视投影。

</details>

### 挑战题

**练习 2.4** 设计一个算法，将任意 $4 \times 4$ 仿射变换矩阵分解为平移、旋转和缩放的组合。假设没有错切变换。

*Hint*: 先提取平移部分，然后对 $3 \times 3$ 子矩阵进行极分解或QR分解。

<details>
<summary>答案</summary>

1. 提取平移：$\mathbf{t} = [M_{14}, M_{24}, M_{34}]^T$
2. 提取 $3 \times 3$ 子矩阵 $\mathbf{A}$
3. 计算缩放：
   - $s_x = ||\mathbf{a}_1||$, $s_y = ||\mathbf{a}_2||$, $s_z = ||\mathbf{a}_3||$
   - 其中 $\mathbf{a}_i$ 是 $\mathbf{A}$ 的第 $i$ 列
4. 计算旋转：
   - $\mathbf{R} = \mathbf{A} \cdot \text{diag}(1/s_x, 1/s_y, 1/s_z)$
   - 如果 $\det(\mathbf{R}) < 0$，则有反射，需要调整一个缩放因子的符号
5. 验证：$\mathbf{M} = \mathbf{T}(\mathbf{t}) \cdot \mathbf{R} \cdot \mathbf{S}(s_x, s_y, s_z)$

</details>

**练习 2.5** 证明 Rodrigues 旋转公式，并说明当 $\theta = 180°$ 时的特殊情况。

*Hint*: 将向量分解为平行和垂直于旋转轴的分量。

<details>
<summary>答案</summary>

设向量 $\mathbf{v}$ 绕单位轴 $\mathbf{n}$ 旋转 $\theta$ 得到 $\mathbf{v}'$。

分解：$\mathbf{v} = \mathbf{v}_\parallel + \mathbf{v}_\perp$
- $\mathbf{v}_\parallel = (\mathbf{n} \cdot \mathbf{v})\mathbf{n}$
- $\mathbf{v}_\perp = \mathbf{v} - \mathbf{v}_\parallel$

旋转后：
- $\mathbf{v}'_\parallel = \mathbf{v}_\parallel$（不变）
- $\mathbf{v}'_\perp = \cos\theta \mathbf{v}_\perp + \sin\theta (\mathbf{n} \times \mathbf{v}_\perp)$

合并得到 Rodrigues 公式。

当 $\theta = 180°$：
$$\mathbf{R}(\mathbf{n}, 180°) = -\mathbf{I} + 2\mathbf{n}\mathbf{n}^T$$
这是关于轴 $\mathbf{n}$ 的反射。

</details>

**练习 2.6** 给定视锥体的六个平面，推导如何在齐次裁剪空间进行视锥体裁剪。

*Hint*: 在裁剪空间中，视锥体变成了 $[-w, w]$ 的立方体。

<details>
<summary>答案</summary>

在齐次裁剪空间，点 $(x, y, z, w)$ 在视锥体内当且仅当：
- $-w \leq x \leq w$（左右平面）
- $-w \leq y \leq w$（上下平面）
- $0 \leq z \leq w$（近远平面，假设使用 $[0,1]$ 深度范围）

裁剪测试：
1. 如果所有顶点都在某个平面外侧，则完全剔除
2. 如果所有顶点都在内侧，则完全保留
3. 否则需要裁剪，使用 Sutherland-Hodgman 算法

优化：使用 6 位掩码编码顶点相对于 6 个平面的位置。

</details>

**练习 2.7**（开放题）讨论在虚拟现实（VR）渲染中，如何修改投影矩阵以适应双目渲染和镜头畸变校正。

*Hint*: 考虑非对称视锥体和逆畸变预处理。

<details>
<summary>答案</summary>

VR渲染的特殊考虑：

1. **非对称投影**：
   - 左右眼的视锥体不是中心对称的
   - 需要偏移投影中心：修改投影矩阵的 $(1,3)$ 和 $(2,3)$ 元素

2. **镜头畸变**：
   - 桶形畸变：$r' = r(1 + k_1r^2 + k_2r^4)$
   - 预畸变：在投影后应用反向畸变
   - 或使用顶点着色器中的非线性变换

3. **性能优化**：
   - 固定注视点渲染（Foveated Rendering）
   - 多分辨率着色（Multi-Resolution Shading）
   - 单通道立体渲染（Single Pass Stereo）

4. **时间扭曲**（Timewarp）：
   - 使用最新的头部姿态更新已渲染图像
   - 需要深度信息进行重投影

</details>

## 常见陷阱与错误

### 1. 变换顺序错误
```
错误：Model × View × Projection
正确：Projection × View × Model
```
记忆技巧：从右到左读，"先模型变换，再视图变换，最后投影"

### 2. 列主序 vs 行主序
- OpenGL：列主序，$\mathbf{v}' = \mathbf{M}\mathbf{v}$
- DirectX：行主序，$\mathbf{v}' = \mathbf{v}\mathbf{M}^T$
- 混淆会导致变换完全错误

### 3. 法向量变换错误
```
错误：n' = M × n
正确：n' = (M^(-1))^T × n
```
只有在仅包含旋转和均匀缩放时才能直接使用模型矩阵

### 4. 深度精度问题
- 近远平面比例不当（如 0.01 到 10000）
- 解决：使用对数深度或反向 Z

### 5. 万向锁（Gimbal Lock）
- 使用欧拉角时的奇异性
- 解决：使用四元数或轴角表示

### 6. 浮点累积误差
- 连续小角度旋转导致矩阵不再正交
- 解决：定期重新正交化

## 最佳实践检查清单

### 设计阶段
- [ ] 选择合适的坐标系（左手/右手）
- [ ] 确定变换顺序和矩阵存储格式
- [ ] 考虑数值精度需求
- [ ] 设计变换缓存策略

### 实现阶段
- [ ] 使用齐次坐标统一表示
- [ ] 预计算静态变换矩阵
- [ ] 正确处理法向量变换
- [ ] 实现矩阵正交化机制
- [ ] 利用SIMD指令优化

### 优化阶段
- [ ] 最小化矩阵乘法次数
- [ ] 使用场景图减少重复计算
- [ ] 实现视锥体裁剪
- [ ] 考虑GPU实例化
- [ ] profile矩阵运算热点

### 调试阶段
- [ ] 验证矩阵正交性
- [ ] 检查变换顺序
- [ ] 测试极端情况（大坐标值、小角度）
- [ ] 验证投影矩阵的近远平面
- [ ] 检查深度精度问题
