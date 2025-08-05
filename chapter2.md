# 第2章：几何变换

几何变换是计算机图形学的核心基础之一，它决定了我们如何在二维和三维空间中移动、旋转、缩放物体，以及如何将三维世界投影到二维屏幕上。本章将深入探讨变换的数学原理、实现细节和优化技巧，为后续的渲染管线学习打下坚实基础。

## 2.1 二维与三维变换

### 2.1.1 齐次坐标系

在计算机图形学中，我们使用齐次坐标（Homogeneous Coordinates）来统一表示点和向量，并使得仿射变换可以用矩阵乘法来表示。齐次坐标是投影几何的核心概念，它优雅地解决了平移变换无法用矩阵乘法表示的问题。

#### 从欧几里得坐标到齐次坐标

对于二维空间：
- 点：$(x, y) \rightarrow (x, y, 1)$ 或更一般地 $(wx, wy, w)$，其中 $w \neq 0$
- 向量：$(x, y) \rightarrow (x, y, 0)$

对于三维空间：
- 点：$(x, y, z) \rightarrow (x, y, z, 1)$ 或更一般地 $(wx, wy, wz, w)$，其中 $w \neq 0$
- 向量：$(x, y, z) \rightarrow (x, y, z, 0)$

#### 齐次坐标的关键性质

1. **等价性**：$(x, y, z, w)$ 和 $(kx, ky, kz, kw)$ 表示同一个点（$k \neq 0$）
   - 这定义了一个等价类，所有满足比例关系的坐标表示同一点
   - 规范化形式：当 $w \neq 0$ 时，可以规范化为 $(x/w, y/w, z/w, 1)$

2. **笛卡尔坐标恢复**：
   - 当 $w \neq 0$ 时：$(x, y, z, w) \rightarrow (x/w, y/w, z/w)$
   - 当 $w = 0$ 时：表示无穷远处的点或方向向量

3. **代数封闭性**：
   - 两条平行线在投影空间中相交于无穷远点
   - 使得某些几何定理（如对偶原理）更加统一和优雅

#### 齐次坐标的几何解释

齐次坐标本质上是将 $n$ 维空间嵌入到 $n+1$ 维投影空间中：

1. **投影模型**：
   - 想象三维空间中通过原点的所有直线
   - 每条直线（除原点外）代表二维平面上的一个点
   - 直线与 $w=1$ 平面的交点就是对应的笛卡尔坐标

2. **无穷远元素**：
   - $w=0$ 的超平面包含所有无穷远点
   - 平行线在此平面上相交，实现了投影几何的完备性
   - 例如：$(1, 0, 0)$ 表示 $x$ 轴方向的无穷远点

3. **对偶性**：
   - 点和线在投影平面上具有对偶关系
   - 线的方程 $ax + by + c = 0$ 可表示为 $[a, b, c]$
   - 点 $(x, y, 1)$ 在线 $[a, b, c]$ 上当且仅当 $ax + by + c = 0$

#### 为什么区分点和向量

这种区分不仅是数学上的优雅，更具有深刻的几何意义：

1. **仿射组合**：
   - 点的仿射组合：$\sum_i \alpha_i \mathbf{p}_i$，当 $\sum_i \alpha_i = 1$ 时结果是点
   - 向量的线性组合：$\sum_i \beta_i \mathbf{v}_i$ 总是向量
   - 这保证了几何运算的类型安全性

2. **变换行为**：
   - 平移对点有效：$\mathbf{p}' = \mathbf{p} + \mathbf{t}$
   - 平移对向量无效：$\mathbf{v}' = \mathbf{v}$（方向不因位置改变）
   - 旋转和缩放对两者都有效

3. **物理意义**：
   - 点表示位置（position）
   - 向量表示位移（displacement）或方向（direction）
   - 两点之差自然产生向量：$\mathbf{v} = \mathbf{p}_2 - \mathbf{p}_1$

#### 齐次坐标的运算规则

1. **加法运算**：
   - $(x_1, y_1, w_1) + (x_2, y_2, w_2) = (x_1 + x_2, y_1 + y_2, w_1 + w_2)$
   - 注意：只有当 $w_1 + w_2 \neq 0$ 时结果才有意义

2. **数乘运算**：
   - $k(x, y, w) = (kx, ky, kw)$
   - 保持点的位置不变（等价类相同）

3. **点的重心坐标**：
   - 三角形顶点 $\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_3$ 的重心坐标
   - $\mathbf{p} = \alpha\mathbf{p}_1 + \beta\mathbf{p}_2 + \gamma\mathbf{p}_3$，其中 $\alpha + \beta + \gamma = 1$
   - 齐次坐标自动保证了 $w$ 分量的正确性

#### 实际应用中的考虑

1. **数值精度**：
   - 避免 $w$ 接近 0 的情况（除法不稳定）
   - 定期规范化齐次坐标（将 $w$ 设为 1）
   - 使用 `glm::perspective` 等库函数时注意精度设置

2. **GPU实现**：
   - 顶点着色器输出 `gl_Position` 是四维齐次坐标
   - 硬件自动执行透视除法（perspective divide）
   - 裁剪在齐次空间进行（效率更高）

3. **调试技巧**：
   - 检查 $w$ 分量是否符合预期（点为1，向量为0）
   - 验证变换后的 $w$ 值（透视投影会改变 $w$）
   - 使用可视化工具观察齐次空间的变换

### 2.1.2 基本二维变换

二维变换构成了计算机图形学的基础，它们可以组合产生复杂的效果。每种变换都有其独特的几何意义和代数性质。

#### 平移变换（Translation）

平移将所有点沿着指定方向移动相同的距离：

$$\mathbf{T}(t_x, t_y) = \begin{bmatrix}
1 & 0 & t_x \\
0 & 1 & t_y \\
0 & 0 & 1
\end{bmatrix}$$

**数学性质**：
- 平移群构成阿贝尔群（交换群）
- 逆变换：$\mathbf{T}^{-1}(t_x, t_y) = \mathbf{T}(-t_x, -t_y)$
- 复合：$\mathbf{T}(t_1)\mathbf{T}(t_2) = \mathbf{T}(t_1 + t_2)$
- 与向量的关系：平移向量 $\mathbf{t} = (t_x, t_y, 0)$

**几何意义**：
- 保持所有几何性质（距离、角度、面积）
- 等距变换（isometry）的一种
- 不改变图形的方向性（orientation）

**实现细节**：
```
点的平移：[x', y', 1]^T = T(tx, ty) × [x, y, 1]^T = [x+tx, y+ty, 1]^T
向量不变：[x', y', 0]^T = T(tx, ty) × [x, y, 0]^T = [x, y, 0]^T
```

#### 旋转变换（Rotation）

绕原点逆时针旋转 $\theta$ 角度的变换矩阵：

$$\mathbf{R}(\theta) = \begin{bmatrix}
\cos\theta & -\sin\theta & 0 \\
\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{bmatrix}$$

**推导过程**：
考虑单位圆上的点 $(r\cos\phi, r\sin\phi)$ 旋转 $\theta$ 后：
$$\begin{align}
x' &= r\cos(\phi + \theta) = r\cos\phi\cos\theta - r\sin\phi\sin\theta \\
y' &= r\sin(\phi + \theta) = r\cos\phi\sin\theta + r\sin\phi\cos\theta
\end{align}$$

这给出了旋转矩阵的形式。

**重要性质**：
1. **正交性**：$\mathbf{R}^T\mathbf{R} = \mathbf{I}$，即 $\mathbf{R}^{-1} = \mathbf{R}^T$
2. **行列式**：$\det(\mathbf{R}) = \cos^2\theta + \sin^2\theta = 1$（保持面积和方向）
3. **特征值**：$\lambda = e^{i\theta}, e^{-i\theta}, 1$（复特征值体现旋转本质）
4. **群性质**：
   - $\mathbf{R}(\alpha)\mathbf{R}(\beta) = \mathbf{R}(\alpha + \beta)$
   - $\mathbf{R}(\theta)^n = \mathbf{R}(n\theta)$
   - $\mathbf{R}(0) = \mathbf{I}$（单位元）

**与复数的联系**：
二维旋转可用复数乘法表示：
$$z' = e^{i\theta}z = (\cos\theta + i\sin\theta)(x + iy)$$

这揭示了旋转的深层数学结构。

**小角度近似**：
当 $\theta$ 很小时：
$$\mathbf{R}(\theta) \approx \begin{bmatrix}
1 & -\theta & 0 \\
\theta & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

这在增量旋转计算中很有用。

#### 缩放变换（Scaling）

缩放变换改变物体的大小：

$$\mathbf{S}(s_x, s_y) = \begin{bmatrix}
s_x & 0 & 0 \\
0 & s_y & 0 \\
0 & 0 & 1
\end{bmatrix}$$

**分类与性质**：
1. **均匀缩放**（$s_x = s_y = s$）：
   - 保持形状和角度
   - 相似变换的一部分
   - 面积变化：$s^2$ 倍

2. **非均匀缩放**（$s_x \neq s_y$）：
   - 改变长宽比
   - 不保持角度（除了与坐标轴平行的角）
   - 椭圆变换的基础

3. **反射缩放**（$s_x < 0$ 或 $s_y < 0$）：
   - 改变方向性（orientation）
   - $\det(\mathbf{S}) < 0$ 时发生镜像
   - 组合使用可产生各种对称

**特征分解**：
缩放矩阵是对角矩阵，其特征向量是坐标轴方向：
- 特征值：$\lambda_1 = s_x, \lambda_2 = s_y, \lambda_3 = 1$
- 特征向量：$(1,0,0)^T, (0,1,0)^T, (0,0,1)^T$

**缩放中心**：
绕任意点 $(c_x, c_y)$ 缩放：
$$\mathbf{S}_c = \mathbf{T}(c_x, c_y) \mathbf{S}(s_x, s_y) \mathbf{T}(-c_x, -c_y)$$

展开后：
$$\mathbf{S}_c = \begin{bmatrix}
s_x & 0 & c_x(1-s_x) \\
0 & s_y & c_y(1-s_y) \\
0 & 0 & 1
\end{bmatrix}$$

#### 错切变换（Shearing）

错切变换使物体发生倾斜：

**水平错切**：
$$\mathbf{H}_x(s) = \begin{bmatrix}
1 & s & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

**垂直错切**：
$$\mathbf{H}_y(s) = \begin{bmatrix}
1 & 0 & 0 \\
s & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

**几何解释**：
- $\mathbf{H}_x(s)$：每个点的 $x$ 坐标增加 $sy$，$y$ 坐标不变
- 垂直线保持垂直，水平线倾斜角度 $\arctan(s)$
- 矩形变成平行四边形

**重要性质**：
1. **面积守恒**：$\det(\mathbf{H}) = 1$
2. **可逆性**：$\mathbf{H}_x^{-1}(s) = \mathbf{H}_x(-s)$
3. **复合错切**：
   $$\mathbf{H}_x(s_1)\mathbf{H}_x(s_2) = \mathbf{H}_x(s_1 + s_2)$$
4. **与旋转的关系**：
   任意旋转可分解为三次错切：
   $$\mathbf{R}(\theta) = \mathbf{H}_x(-\tan\frac{\theta}{2})\mathbf{H}_y(\sin\theta)\mathbf{H}_x(-\tan\frac{\theta}{2})$$

**应用场景**：
- 斜体文字渲染
- 图像倾斜校正
- 仿射纹理映射

#### 反射变换（Reflection）

反射产生镜像效果：

**基本反射**：
- 关于 $x$ 轴：$\mathbf{F}_x = \text{diag}(1, -1, 1)$
- 关于 $y$ 轴：$\mathbf{F}_y = \text{diag}(-1, 1, 1)$
- 关于原点：$\mathbf{F}_o = \text{diag}(-1, -1, 1)$

**关于任意直线的反射**：
对于直线 $\mathbf{n} \cdot \mathbf{p} = d$（$\mathbf{n}$ 是单位法向量）：

$$\mathbf{F} = \mathbf{I} - 2\mathbf{n}\mathbf{n}^T$$

对于直线 $ax + by + c = 0$：
1. 归一化：$\mathbf{n} = \frac{1}{\sqrt{a^2 + b^2}}(a, b)$
2. 构造反射矩阵：
$$\mathbf{F} = \begin{bmatrix}
1-2n_x^2 & -2n_xn_y & -2cn_x/\sqrt{a^2+b^2} \\
-2n_xn_y & 1-2n_y^2 & -2cn_y/\sqrt{a^2+b^2} \\
0 & 0 & 1
\end{bmatrix}$$

**性质**：
1. **对合性**：$\mathbf{F}^2 = \mathbf{I}$（反射两次回到原位）
2. **行列式**：$\det(\mathbf{F}) = -1$（改变方向性）
3. **特征值**：$\lambda = 1, -1, 1$
   - 特征值1对应反射轴上的不动点
   - 特征值-1对应垂直于反射轴的方向

**组合反射**：
- 两个反射的复合是旋转
- 三个反射的复合可能是反射或滑移反射
- 这构成了二维等距变换群的生成元

#### 变换的分类与层次

**保持性质的层次结构**：
1. **欧几里得变换**（刚体变换）：
   - 保持：距离、角度、面积、方向性
   - 包含：平移、旋转
   
2. **相似变换**：
   - 保持：角度、形状比例
   - 包含：欧几里得变换 + 均匀缩放
   
3. **仿射变换**：
   - 保持：平行性、面积比
   - 包含：相似变换 + 错切 + 非均匀缩放
   
4. **投影变换**：
   - 保持：直线性、交比
   - 最一般的线性变换

### 2.1.3 基本三维变换

三维变换将二维概念扩展到三维空间，但引入了更多的复杂性，特别是在旋转表示和组合方面。

#### 三维平移

三维平移扩展了二维情况：

$$\mathbf{T}(t_x, t_y, t_z) = \begin{bmatrix}
1 & 0 & 0 & t_x \\
0 & 1 & 0 & t_y \\
0 & 0 & 1 & t_z \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**数学性质**：
- 构成三维阿贝尔群：$\mathbf{T}(\mathbf{a})\mathbf{T}(\mathbf{b}) = \mathbf{T}(\mathbf{a}+\mathbf{b})$
- 逆变换：$\mathbf{T}^{-1}(t_x, t_y, t_z) = \mathbf{T}(-t_x, -t_y, -t_z)$
- 与旋转的交换性：$\mathbf{T}\mathbf{R} \neq \mathbf{R}\mathbf{T}$（一般情况）

**应用场景**：
- 物体世界定位
- 相机运动（与视图矩阵相关）
- 粒子系统的位置更新
- 骨骼动画的关节偏移

#### 三维缩放

三维缩放影响物体的体积和比例：

$$\mathbf{S}(s_x, s_y, s_z) = \begin{bmatrix}
s_x & 0 & 0 & 0 \\
0 & s_y & 0 & 0 \\
0 & 0 & s_z & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**几何意义**：
1. **体积变化**：$V' = |s_x \cdot s_y \cdot s_z| \cdot V$
2. **表面积变化**：取决于具体形状，如球体：$A' = (s_xs_y + s_ys_z + s_zs_x)/3 \cdot A$
3. **法向量变换**：需要使用逆转置矩阵

**特殊缩放类型**：
1. **均匀缩放**（$s_x = s_y = s_z = s$）：
   - 保持形状和角度
   - 用于LOD（Level of Detail）系统
   - 透视投影中的深度缩放

2. **平面缩放**（如 $s_z = 0$）：
   - 将三维物体投影到平面
   - 阴影生成的基础
   - 注意：产生奇异矩阵

3. **各向异性缩放**：
   - 模拟物理变形（挤压、拉伸）
   - 椭球体生成
   - 纹理空间变换

**缩放中心的选择**：
```
绕点P缩放 = T(P) × S(sx,sy,sz) × T(-P)
常见中心：物体中心、包围盒中心、质心
```

#### 三维旋转

三维旋转是最复杂的基本变换，有多种表示方法。

**基本轴旋转**

绕 x 轴旋转（俯仰 Pitch）：
$$\mathbf{R}_x(\theta) = \begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & \cos\theta & -\sin\theta & 0 \\
0 & \sin\theta & \cos\theta & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

绕 y 轴旋转（偏航 Yaw）：
$$\mathbf{R}_y(\theta) = \begin{bmatrix}
\cos\theta & 0 & \sin\theta & 0 \\
0 & 1 & 0 & 0 \\
-\sin\theta & 0 & \cos\theta & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

绕 z 轴旋转（滚转 Roll）：
$$\mathbf{R}_z(\theta) = \begin{bmatrix}
\cos\theta & -\sin\theta & 0 & 0 \\
\sin\theta & \cos\theta & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**记忆规则**：
1. 旋转轴对应的行列保持单位向量形式
2. 其他2×2子矩阵是二维旋转矩阵
3. $\mathbf{R}_y$ 的符号相反（保持右手系）

**旋转矩阵的性质**：
1. **正交性**：$\mathbf{R}^T\mathbf{R} = \mathbf{I}$
2. **行列式**：$\det(\mathbf{R}) = 1$（保持手性）
3. **特征值**：$\lambda = 1, e^{i\theta}, e^{-i\theta}$
   - 实特征值1对应旋转轴
   - 复特征值对描述旋转角度
4. **群结构**：SO(3)特殊正交群
   - 封闭性：旋转的复合仍是旋转
   - 非交换：$\mathbf{R}_1\mathbf{R}_2 \neq \mathbf{R}_2\mathbf{R}_1$

**欧拉角表示**

任意旋转可分解为三次基本旋转：

1. **ZYX顺序**（偏航-俯仰-滚转）：
   $$\mathbf{R} = \mathbf{R}_z(\psi)\mathbf{R}_y(\theta)\mathbf{R}_x(\phi)$$

2. **XYZ顺序**（滚转-俯仰-偏航）：
   $$\mathbf{R} = \mathbf{R}_x(\phi)\mathbf{R}_y(\theta)\mathbf{R}_z(\psi)$$

**万向锁（Gimbal Lock）问题**：
- 当中间旋转角为 ±90° 时发生
- 第一和第三次旋转轴对齐
- 失去一个旋转自由度
- 导致动画不连续

**解决方案**：
1. 限制旋转范围
2. 使用四元数
3. 轴角表示
4. 旋转矩阵直接插值

#### 三维反射

三维反射扩展了二维情况：

**基本平面反射**：
- XY平面（z=0）：$\mathbf{F}_{xy} = \text{diag}(1, 1, -1, 1)$
- YZ平面（x=0）：$\mathbf{F}_{yz} = \text{diag}(-1, 1, 1, 1)$
- XZ平面（y=0）：$\mathbf{F}_{xz} = \text{diag}(1, -1, 1, 1)$

**任意平面反射**：
对于平面 $\mathbf{n} \cdot \mathbf{p} = d$（$\mathbf{n}$ 是单位法向量）：

$$\mathbf{F} = \mathbf{I} - 2\mathbf{n}\mathbf{n}^T$$

展开为4×4矩阵：
$$\mathbf{F} = \begin{bmatrix}
1-2n_x^2 & -2n_xn_y & -2n_xn_z & -2dn_x \\
-2n_yn_x & 1-2n_y^2 & -2n_yn_z & -2dn_y \\
-2n_zn_x & -2n_zn_y & 1-2n_z^2 & -2dn_z \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**应用**：
- 镜面效果
- 对称建模
- 阴影生成（平面投影）

#### 三维错切

三维错切比二维复杂，有多种形式：

**XY错切**（沿z变化）：
$$\mathbf{H}_{xy}(h_x, h_y) = \begin{bmatrix}
1 & 0 & h_x & 0 \\
0 & 1 & h_y & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**一般错切矩阵**：
$$\mathbf{H} = \begin{bmatrix}
1 & h_{xy} & h_{xz} & 0 \\
h_{yx} & 1 & h_{yz} & 0 \\
h_{zx} & h_{zy} & 1 & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}$$

**性质**：
- 保持体积：$\det(\mathbf{H}) = 1$
- 平行平面保持平行
- 可分解为基本错切的组合

#### 变换的组合与优化

**组合顺序的重要性**：
```
局部到世界：Scale → Rotate → Translate
世界到局部：Translate^(-1) → Rotate^(-1) → Scale^(-1)
```

**常见组合模式**：
1. **TRS变换**：$\mathbf{M} = \mathbf{T}\mathbf{R}\mathbf{S}$
   - 最常用的组合顺序
   - 便于独立控制各分量
   - 适合关键帧动画

2. **枢轴点变换**：
   $$\mathbf{M} = \mathbf{T}_{\text{pos}}\mathbf{T}_{\text{pivot}}\mathbf{R}\mathbf{S}\mathbf{T}_{-\text{pivot}}$$
   - 支持任意点为中心的变换
   - 用于3D软件的变换工具

3. **层次变换**：
   $$\mathbf{M}_{\text{child}} = \mathbf{M}_{\text{parent}} \mathbf{M}_{\text{local}}$$
   - 场景图结构
   - 骨骼动画系统
   - 相对坐标系

**性能优化考虑**：
1. **矩阵缓存**：
   - 预计算静态变换
   - 缓存中间结果
   - 避免重复计算

2. **SIMD优化**：
   - 4×4矩阵适合SSE/AVX
   - 批量顶点变换
   - 并行矩阵乘法

3. **精度维护**：
   - 定期正交化旋转部分
   - 使用双精度关键计算
   - 避免极小缩放值

### 2.1.4 任意轴旋转（Rodrigues' Rotation Formula）

绕任意单位向量 $\mathbf{n} = (n_x, n_y, n_z)$ 旋转 $\theta$ 角度：

$$\mathbf{R}(\mathbf{n}, \theta) = \cos\theta \mathbf{I} + (1-\cos\theta)\mathbf{n}\mathbf{n}^T + \sin\theta \mathbf{N}$$

其中 $\mathbf{N}$ 是 $\mathbf{n}$ 的反对称矩阵（叉积矩阵）：

$$\mathbf{N} = \begin{bmatrix}
0 & -n_z & n_y \\
n_z & 0 & -n_x \\
-n_y & n_x & 0
\end{bmatrix}$$

**几何推导**：
对于向量 $\mathbf{v}$，将其分解为平行和垂直于旋转轴的分量：
- $\mathbf{v}_\parallel = (\mathbf{n} \cdot \mathbf{v})\mathbf{n}$（投影）
- $\mathbf{v}_\perp = \mathbf{v} - \mathbf{v}_\parallel$

旋转后：
- $\mathbf{v}'_\parallel = \mathbf{v}_\parallel$（不变）
- $\mathbf{v}'_\perp = \cos\theta \mathbf{v}_\perp + \sin\theta (\mathbf{n} \times \mathbf{v}_\perp)$

**矩阵形式展开**：
$$\mathbf{R}(\mathbf{n}, \theta) = \begin{bmatrix}
c + n_x^2(1-c) & n_xn_y(1-c) - n_zs & n_xn_z(1-c) + n_ys \\
n_yn_x(1-c) + n_zs & c + n_y^2(1-c) & n_yn_z(1-c) - n_xs \\
n_zn_x(1-c) - n_ys & n_zn_y(1-c) + n_xs & c + n_z^2(1-c)
\end{bmatrix}$$

其中 $c = \cos\theta$，$s = \sin\theta$。

**与四元数的联系**：
四元数 $q = \cos(\theta/2) + \sin(\theta/2)(n_xi + n_yj + n_zk)$ 对应同样的旋转。

**计算优化**：
- 预计算 $\cos\theta$、$\sin\theta$、$(1-\cos\theta)$
- 利用对称性减少乘法次数
- 小角度近似：$\cos\theta \approx 1$，$\sin\theta \approx \theta$

### 2.1.5 变换的分解与组合

**仿射变换的一般形式**：
$$\mathbf{A} = \begin{bmatrix}
\mathbf{M}_{3\times3} & \mathbf{t} \\
\mathbf{0}^T & 1
\end{bmatrix}$$

其中 $\mathbf{M}_{3\times3}$ 可分解为：
$$\mathbf{M} = \mathbf{R} \cdot \mathbf{S} \cdot \mathbf{H}$$

**变换的分类**：
1. **刚体变换（Rigid Transform）**：保持距离和角度
   - 组成：旋转 + 平移
   - 性质：$\mathbf{M} = \mathbf{T} \cdot \mathbf{R}$
   - 自由度：6（3个平移 + 3个旋转）

2. **相似变换（Similarity Transform）**：保持角度和形状比例
   - 组成：均匀缩放 + 旋转 + 平移
   - 性质：$\mathbf{M} = \mathbf{T} \cdot \mathbf{R} \cdot \mathbf{S}(s)$
   - 自由度：7

3. **仿射变换（Affine Transform）**：保持平行性
   - 组成：缩放 + 错切 + 旋转 + 平移
   - 性质：直线映射为直线，平行线保持平行
   - 自由度：12

4. **投影变换（Projective Transform）**：保持直线性
   - 最一般的线性变换
   - 自由度：15（$4\times4$ 矩阵减去齐次缩放）

**变换顺序的重要性**：

先平移后旋转：
$$\mathbf{M}_1 = \mathbf{R} \cdot \mathbf{T} = \begin{bmatrix}
\mathbf{R}_{3\times3} & \mathbf{R}_{3\times3}\mathbf{t} \\
\mathbf{0}^T & 1
\end{bmatrix}$$

先旋转后平移：
$$\mathbf{M}_2 = \mathbf{T} \cdot \mathbf{R} = \begin{bmatrix}
\mathbf{R}_{3\times3} & \mathbf{t} \\
\mathbf{0}^T & 1
\end{bmatrix}$$

两者结果完全不同！

**绕任意点旋转**：
绕点 $\mathbf{c}$ 旋转 $\theta$ 角度：
$$\mathbf{M} = \mathbf{T}(\mathbf{c}) \cdot \mathbf{R}(\theta) \cdot \mathbf{T}(-\mathbf{c})$$

这是变换组合的经典应用：先平移到原点，旋转，再平移回去。

**变换的插值**：
线性插值两个变换时，不能直接插值矩阵元素：
- 平移：可以线性插值
- 旋转：使用四元数球面插值（SLERP）
- 缩放：使用对数空间插值

## 2.2 模型、视图、投影变换

### 2.2.1 坐标系统与变换管线

图形渲染管线中的坐标系统：

1. **模型空间（Model Space）**
   - 定义：物体的局部坐标系，通常以物体中心或关键点为原点
   - 用途：建模、动画骨骼绑定
   - 特点：每个物体有独立的模型空间

2. **世界空间（World Space）**
   - 定义：场景的全局坐标系，所有物体的统一参考系
   - 用途：物体定位、光照计算、物理模拟
   - 坐标系选择：通常 Y-up（OpenGL）或 Z-up（某些CAD系统）

3. **视图空间（View/Camera Space）**
   - 定义：以相机为原点，-Z 为视线方向的坐标系
   - 用途：视锥体裁剪、深度排序
   - 约定：右手系（OpenGL）或左手系（DirectX）

4. **裁剪空间（Clip Space）**
   - 定义：投影变换后的齐次坐标系
   - 范围：未进行透视除法，$w$ 分量包含深度信息
   - 用途：视锥体裁剪、透视插值

5. **NDC空间（Normalized Device Coordinates）**
   - 定义：透视除法后的归一化坐标
   - 范围：$[-1,1]^3$（OpenGL）或 $[-1,1]^2 \times [0,1]$（DirectX）
   - 特点：与设备无关的标准化坐标

6. **屏幕空间（Screen Space）**
   - 定义：最终的像素坐标
   - 范围：$[0,W] \times [0,H] \times [0,1]$
   - 用途：光栅化、片元着色

**完整的变换管线**：
$$\text{顶点} \xrightarrow{\mathbf{M}} \text{世界} \xrightarrow{\mathbf{V}} \text{视图} \xrightarrow{\mathbf{P}} \text{裁剪} \xrightarrow{\text{÷w}} \text{NDC} \xrightarrow{\mathbf{VP}} \text{屏幕}$$

**坐标系的手性（Handedness）**：
- 右手系：拇指=X，食指=Y，中指=Z（OpenGL默认）
- 左手系：Z轴方向相反（DirectX默认）
- 转换：通过缩放 $\mathbf{S}(1,1,-1)$ 实现

### 2.2.2 模型变换（Model Transform）

模型变换将物体从模型空间变换到世界空间：
$$\mathbf{M}_{\text{model}} = \mathbf{T}_{\text{world}} \cdot \mathbf{R}_{\text{world}} \cdot \mathbf{S}_{\text{local}}$$

**模型变换的组成**：
1. **局部缩放**：调整物体大小
2. **局部旋转**：设置物体朝向
3. **世界平移**：放置物体位置

**层次化变换（Hierarchical Transforms）**：
场景图中的父子关系：
$$\mathbf{M}_{\text{child}} = \mathbf{M}_{\text{parent}} \cdot \mathbf{M}_{\text{local}}$$

例如，机械臂的变换链：
- 基座 → 上臂 → 前臂 → 手掌 → 手指

**实例化渲染（Instanced Rendering）**：
对于大量相同模型（如森林中的树木）：
```
基础模型矩阵 + 每实例变换数据 = 最终变换
```

优化技巧：
- GPU实例化：使用实例缓冲区存储每个实例的变换
- LOD（细节层次）：根据距离选择不同精度的模型
- 视锥体剔除：只变换可见物体

**骨骼动画的模型变换**：
顶点最终位置由多个骨骼影响：
$$\mathbf{v}' = \sum_{i} w_i \mathbf{M}_{\text{bone},i} \mathbf{M}_{\text{bind},i}^{-1} \mathbf{v}$$

其中：
- $w_i$：骨骼权重
- $\mathbf{M}_{\text{bone},i}$：骨骼当前变换
- $\mathbf{M}_{\text{bind},i}^{-1}$：绑定姿态的逆变换

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
