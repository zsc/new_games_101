# 第7章：光线追踪基础

光线追踪是计算机图形学中实现真实感渲染的核心技术之一。不同于光栅化的前向渲染方式，光线追踪采用逆向追踪光线的方式，从相机出发追踪光线在场景中的传播路径。本章将深入探讨光线追踪的基本原理、加速结构的设计与实现，以及光线与物体相交算法的优化技术。我们将重点关注算法的数学基础和性能优化策略，为后续的全局光照和高级渲染技术打下坚实基础。

## 7.1 光线追踪基本原理

### 7.1.1 光线的数学表示

光线可以用参数方程表示：
$$\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$$

其中：
- $\mathbf{o}$ 是光线起点（origin）
- $\mathbf{d}$ 是归一化的方向向量（direction）
- $t \geq 0$ 是参数，表示沿光线方向的距离

### 7.1.2 基本光线追踪算法

光线追踪的核心步骤：

1. **主光线生成（Primary Ray Generation）**
   对于图像平面上的每个像素 $(i, j)$，生成从相机位置出发的主光线：
   $$\mathbf{d}_{i,j} = \text{normalize}(\mathbf{p}_{i,j} - \mathbf{e})$$
   
   其中 $\mathbf{e}$ 是相机位置，$\mathbf{p}_{i,j}$ 是像素在世界坐标系中的位置。

2. **光线-场景相交测试**
   对每条光线，找到最近的相交点：
   $$t_{\text{hit}} = \min_{k \in \text{objects}} t_k$$
   
   其中 $t_k$ 是光线与第 $k$ 个物体的相交参数。

3. **着色计算**
   在相交点处计算颜色，考虑：
   - 直接光照（使用阴影光线）
   - 反射（递归追踪反射光线）
   - 折射（对透明物体）

### 7.1.3 递归光线追踪

递归光线追踪通过追踪次级光线来模拟复杂的光照效果：

$$L_o(\mathbf{x}, \omega_o) = L_e(\mathbf{x}, \omega_o) + \int_{\Omega} f_r(\mathbf{x}, \omega_i, \omega_o) L_i(\mathbf{x}, \omega_i) (\omega_i \cdot \mathbf{n}) d\omega_i$$

实际实现中，通常限制递归深度或使用俄罗斯轮盘赌（Russian Roulette）来终止递归。

### 7.1.4 相机模型与光线生成

**透视投影相机**的光线生成：

给定视场角 $\text{fov}$、宽高比 $\text{aspect}$，像素 $(i, j)$ 对应的光线方向：

$$\begin{aligned}
u &= \frac{2i - \text{width}}{\text{width}} \cdot \text{aspect} \cdot \tan(\text{fov}/2) \\
v &= \frac{2j - \text{height}}{\text{height}} \cdot \tan(\text{fov}/2) \\
\mathbf{d} &= \text{normalize}(u\mathbf{right} + v\mathbf{up} - \mathbf{forward})
\end{aligned}$$

## 7.2 加速结构（BVH、KD-Tree）

### 7.2.1 空间数据结构的必要性

朴素的光线追踪需要测试每条光线与场景中所有物体的相交，复杂度为 $O(N)$。对于包含百万级三角形的场景，这是不可接受的。空间加速结构可以将平均复杂度降低到 $O(\log N)$。

### 7.2.2 层次包围盒（BVH）

**BVH构建算法**：

1. **自顶向下构建**
   ```
   选择分割轴和分割位置
   将物体分为两组
   递归构建子树
   ```

2. **SAH（Surface Area Heuristic）**
   最优分割的代价函数：
   $$C = C_{\text{trav}} + \frac{A_L}{A} \cdot N_L \cdot C_{\text{isect}} + \frac{A_R}{A} \cdot N_R \cdot C_{\text{isect}}$$
   
   其中：
   - $C_{\text{trav}}$ 是遍历节点的代价
   - $C_{\text{isect}}$ 是相交测试的代价
   - $A_L, A_R$ 是左右子节点的表面积
   - $N_L, N_R$ 是左右子节点包含的图元数

3. **BVH遍历**
   使用栈式遍历或递归遍历：
   ```
   如果是叶节点：测试所有图元
   否则：
     测试左右子节点的包围盒
     优先遍历较近的子节点
   ```

### 7.2.3 KD-Tree

**KD-Tree特点**：
- 空间分割而非物体分割
- 分割平面与坐标轴对齐
- 可能需要处理跨越分割平面的三角形

**构建策略**：
1. 选择分割轴（循环或基于方差）
2. 选择分割位置（中位数或SAH）
3. 递归构建，直到满足终止条件

**遍历优化**：
- 使用 Mailboxing 避免重复相交测试
- 早期终止（当找到相交点后）

### 7.2.4 加速结构比较

| 特性 | BVH | KD-Tree |
|-----|-----|---------|
| 构建复杂度 | $O(N\log N)$ | $O(N\log N)$ |
| 内存占用 | 较低 | 较高 |
| 遍历效率 | 中等 | 高 |
| 动态更新 | 支持 | 困难 |
| GPU友好性 | 高 | 低 |

## 7.3 光线-物体相交算法优化

### 7.3.1 光线-三角形相交

**Möller-Trumbore算法**：

给定三角形顶点 $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$，光线 $\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$：

$$\begin{aligned}
\mathbf{e}_1 &= \mathbf{v}_1 - \mathbf{v}_0 \\
\mathbf{e}_2 &= \mathbf{v}_2 - \mathbf{v}_0 \\
\mathbf{h} &= \mathbf{d} \times \mathbf{e}_2 \\
a &= \mathbf{e}_1 \cdot \mathbf{h}
\end{aligned}$$

如果 $|a| < \epsilon$，光线与三角形平行。否则：

$$\begin{aligned}
f &= 1/a \\
\mathbf{s} &= \mathbf{o} - \mathbf{v}_0 \\
u &= f(\mathbf{s} \cdot \mathbf{h})
\end{aligned}$$

如果 $u < 0$ 或 $u > 1$，无相交。继续：

$$\begin{aligned}
\mathbf{q} &= \mathbf{s} \times \mathbf{e}_1 \\
v &= f(\mathbf{d} \cdot \mathbf{q})
\end{aligned}$$

如果 $v < 0$ 或 $u + v > 1$，无相交。否则：

$$t = f(\mathbf{e}_2 \cdot \mathbf{q})$$

### 7.3.2 光线-包围盒相交

**Slab方法**：

对于轴对齐包围盒（AABB），计算光线与每对平行平面的相交：

$$\begin{aligned}
t_{\text{min}} &= \max(t_{x,\text{min}}, t_{y,\text{min}}, t_{z,\text{min}}) \\
t_{\text{max}} &= \min(t_{x,\text{max}}, t_{y,\text{max}}, t_{z,\text{max}})
\end{aligned}$$

相交条件：$t_{\text{min}} \leq t_{\text{max}}$ 且 $t_{\text{max}} \geq 0$

**优化技巧**：
- 预计算 $1/\mathbf{d}$ 避免除法
- 使用 SSE/AVX 指令并行计算
- 提前退出机制

### 7.3.3 光线-球体相交

球体方程：$||\mathbf{p} - \mathbf{c}||^2 = r^2$

代入光线方程得到二次方程：
$$at^2 + bt + c = 0$$

其中：
$$\begin{aligned}
a &= \mathbf{d} \cdot \mathbf{d} = 1 \text{ (归一化方向)} \\
b &= 2\mathbf{d} \cdot (\mathbf{o} - \mathbf{c}) \\
c &= ||\mathbf{o} - \mathbf{c}||^2 - r^2
\end{aligned}$$

判别式 $\Delta = b^2 - 4ac$，相交参数：
$$t = \frac{-b \pm \sqrt{\Delta}}{2a}$$

### 7.3.4 相交测试优化策略

1. **早期拒绝（Early Rejection）**
   - 使用简单包围体进行预测试
   - 利用空间连贯性

2. **SIMD并行化**
   - 同时测试多条光线
   - 批量处理三角形

3. **缓存优化**
   - 数据结构对齐
   - 热数据分离

4. **精度考虑**
   - 使用稳定的数值算法
   - 处理自相交问题（shadow acne）

## 本章小结

本章介绍了光线追踪的核心概念和关键技术：

**核心概念**：
- 光线的参数表示：$\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$
- 递归光线追踪实现全局光照效果
- 空间加速结构将复杂度从 $O(N)$ 降至 $O(\log N)$

**关键算法**：
- Möller-Trumbore 三角形相交算法
- SAH启发式用于构建最优BVH
- Slab方法用于包围盒相交测试

**性能优化**：
- 使用BVH或KD-Tree加速结构
- SIMD并行化相交测试
- 缓存友好的数据布局

**重要公式**：
- SAH代价函数：$C = C_{\text{trav}} + \sum_i \frac{A_i}{A} N_i C_{\text{isect}}$
- 透视投影光线生成
- 各种几何体的相交测试公式

光线追踪为后续的全局光照、路径追踪等高级渲染技术奠定了基础。掌握本章内容是理解现代渲染技术的关键。

## 练习题

### 基础题

**练习7.1** 推导正交投影相机的光线生成公式。

<details>
<summary>提示</summary>
正交投影中，所有光线方向相同，只是起点不同。考虑视图体积的定义。
</details>

<details>
<summary>答案</summary>

对于正交投影，光线方向都是 $\mathbf{d} = -\mathbf{forward}$（相机前方）。

光线起点：
$$\mathbf{o}_{i,j} = \mathbf{e} + u\mathbf{right} + v\mathbf{up}$$

其中：
$$\begin{aligned}
u &= \frac{2i - \text{width}}{\text{width}} \cdot \frac{\text{width}_{\text{view}}}{2} \\
v &= \frac{2j - \text{height}}{\text{height}} \cdot \frac{\text{height}_{\text{view}}}{2}
\end{aligned}$$

$\text{width}_{\text{view}}$ 和 $\text{height}_{\text{view}}$ 是视图体积的宽高。
</details>

**练习7.2** 给定光线 $\mathbf{o} = (0, 0, 0)$, $\mathbf{d} = (1, 0, 0)$ 和三角形顶点 $\mathbf{v}_0 = (2, -1, -1)$, $\mathbf{v}_1 = (2, 1, -1)$, $\mathbf{v}_2 = (2, 0, 1)$，计算相交点和重心坐标。

<details>
<summary>提示</summary>
使用 Möller-Trumbore 算法或先计算平面方程。
</details>

<details>
<summary>答案</summary>

三角形在平面 $x = 2$ 上。光线沿 $x$ 轴正方向，所以 $t = 2$。

相交点：$\mathbf{p} = (2, 0, 0)$

验证重心坐标：设 $\mathbf{p} = u\mathbf{v}_0 + v\mathbf{v}_1 + w\mathbf{v}_2$，其中 $u + v + w = 1$。

解得：$u = 0$, $v = 0.5$, $w = 0.5$

因此相交点在边 $\mathbf{v}_1\mathbf{v}_2$ 的中点。
</details>

**练习7.3** 证明BVH中使用SAH构建的树，期望遍历代价是最优的。

<details>
<summary>提示</summary>
考虑光线均匀分布的假设，以及条件概率。
</details>

<details>
<summary>答案</summary>

假设光线均匀分布，击中子节点的概率正比于其表面积。

期望代价：
$$E[C] = P(\text{hit}) \cdot C_{\text{hit}} + P(\text{miss}) \cdot C_{\text{miss}}$$

对于内部节点：
$$E[C] = C_{\text{trav}} + P(L) \cdot E[C_L] + P(R) \cdot E[C_R]$$

其中 $P(L) = A_L/A$，$P(R) = A_R/A$。

SAH正是最小化这个期望代价的贪心策略。通过递归论证，可以证明局部最优导致全局最优（在独立性假设下）。
</details>

**练习7.4** 设计一个算法，检测光线是否穿过由多个三角形组成的封闭网格。

<details>
<summary>提示</summary>
考虑奇偶规则（odd-even rule）或计算有符号体积。
</details>

<details>
<summary>答案</summary>

方法1：计数法
- 从光线起点向任意方向发射测试光线
- 计算与网格的相交次数
- 奇数次相交表示起点在内部

方法2：有符号体积法
- 对每个三角形，计算光线起点与三角形构成的四面体有符号体积
- 累加所有体积
- 非零值表示在内部

方法3：法向一致性
- 检查所有相交点处的法向与光线方向的点积符号
- 一致的符号模式可判断内外
</details>

### 挑战题

**练习7.5** 推导四叉BVH（4个子节点）相对于二叉BVH的理论优势和劣势。在什么情况下四叉BVH更优？

<details>
<summary>提示</summary>
考虑树的深度、节点访问次数、SIMD利用率和缓存性能。
</details>

<details>
<summary>答案</summary>

优势：
1. 树深度减少：$\log_4 N = \frac{1}{2}\log_2 N$
2. SIMD友好：可以同时测试4个包围盒
3. 减少遍历开销：更少的节点访问

劣势：
1. 构建复杂：需要考虑更多分割组合
2. 节点更大：每个节点存储4个子节点信息
3. 包围盒可能更松：难以找到最优4路分割

适用场景：
- GPU实现（SIMD宽度大）
- 场景较为均匀分布
- 内存带宽充足

理论分析：设节点访问代价为 $C_n$，包围盒测试代价为 $C_b$。
- 二叉：$C_{\text{total}} = \log_2 N \cdot C_n + 2\log_2 N \cdot C_b$
- 四叉：$C_{\text{total}} = \frac{1}{2}\log_2 N \cdot C_n + 2\log_2 N \cdot C_b$

当 $C_n \gg C_b$ 时，四叉BVH更优。
</details>

**练习7.6** 设计并分析一种自适应的光线追踪算法，能够根据场景复杂度动态调整采样密度。

<details>
<summary>提示</summary>
考虑图像空间的梯度、几何复杂度和着色复杂度。
</details>

<details>
<summary>答案</summary>

自适应采样策略：

1. **初始稀疏采样**
   - 以低分辨率（如每4×4像素1个样本）进行初始采样
   - 记录每个样本的：深度、法线、材质ID、颜色

2. **复杂度估计**
   ```
   复杂度 = w1·深度方差 + w2·法线差异 + w3·颜色梯度 + w4·材质边界
   ```

3. **自适应细分**
   - 高复杂度区域：增加采样密度
   - 使用四叉树结构管理采样点
   - 最大细分级别限制

4. **插值重建**
   - 平滑区域：双线性插值
   - 边缘区域：最近邻或引导滤波

5. **性能分析**
   - 最坏情况：$O(N)$（全分辨率）
   - 最好情况：$O(N/16)$（4×4块）
   - 实际：通常节省50-80%的光线

关键优化：
- 使用GPU的分层深度缓冲
- 时间连贯性：复用前帧信息
- 并行化：基于tile的处理
</details>

**练习7.7** 分析光线追踪中的数值精度问题，并提出一套完整的解决方案。

<details>
<summary>提示</summary>
考虑浮点误差累积、自相交、薄物体穿透等问题。
</details>

<details>
<summary>答案</summary>

主要精度问题：

1. **自相交（Shadow Acne）**
   - 原因：相交点计算的浮点误差
   - 解决：偏移光线起点
   ```
   ε = 1e-4 * max(|x|, |y|, |z|)
   origin_offset = ε * normal
   ```

2. **光线起点误差传播**
   - 使用误差界限追踪
   - Welzl的误差分析：$\delta = \gamma_n \cdot |t|$
   - 其中 $\gamma_n = \frac{n\epsilon}{1-n\epsilon}$

3. **薄物体穿透**
   - 双面测试
   - 保守包围盒扩展
   - 使用区间算术

4. **大场景的精度损失**
   - 局部坐标系变换
   - 分层精度表示
   - 使用双精度关键计算

5. **完整解决方案**
   ```
   结构：
   - 使用相对坐标系
   - 包围盒适当扩展：box.min -= ε, box.max += ε
   - 相交测试使用稳定算法
   
   算法：
   - Kahan求和处理累积
   - 使用有理数算术验证关键决策
   - 自适应精度：远处物体降低精度要求
   ```

误差分析：
总误差 ≤ 相交误差 + 变换误差 + 着色误差
典型值：~1e-3 到 1e-5 相对误差
</details>

**练习7.8** 探讨如何将光线追踪扩展到非欧几何空间（如球面几何或双曲几何）。

<details>
<summary>提示</summary>
考虑测地线、平行传输和曲率的影响。
</details>

<details>
<summary>答案</summary>

非欧几何光线追踪：

1. **球面几何（正曲率）**
   - 光线沿大圆路径传播
   - 参数化：$\mathbf{r}(t) = \cos(t)\mathbf{o} + \sin(t)\mathbf{d}$
   - 相交测试需要考虑周期性

2. **双曲几何（负曲率）**
   - 使用Poincaré球模型或双曲面模型
   - 测地线方程：
   $$\frac{d^2x^i}{dt^2} + \Gamma^i_{jk}\frac{dx^j}{dt}\frac{dx^k}{dt} = 0$$

3. **算法修改**
   ```
   光线传播：
   - 使用数值积分求解测地线方程
   - 自适应步长控制
   
   相交测试：
   - 将物体变换到局部欧式坐标系
   - 或直接在曲面坐标系求解
   
   加速结构：
   - 需要考虑空间的拓扑结构
   - BVH需要使用测地距离
   ```

4. **应用场景**
   - 相对论可视化
   - 全景渲染
   - 艺术效果

5. **性能考虑**
   - 预计算测地线查找表
   - 使用局部近似
   - GPU上的并行积分器
</details>

## 常见陷阱与错误

### 1. 数值精度问题

**问题**：自相交导致的阴影痤疮（Shadow Acne）
```
错误：shadow_ray.origin = hit_point
正确：shadow_ray.origin = hit_point + ε * normal
```

**问题**：光线方向未归一化导致t值含义错误
```
检查：assert(|ray.direction| ≈ 1.0)
```

### 2. 性能陷阱

**问题**：BVH构建时未考虑空间局部性
- 使用Morton编码改善缓存性能
- 考虑节点大小与缓存行对齐

**问题**：过早优化导致代码复杂
- 先实现正确的算法
- 使用性能分析工具定位瓶颈

### 3. 算法错误

**问题**：BVH遍历时错误的节点访问顺序
```
错误：总是先访问左子节点
正确：先访问距离较近的子节点
```

**问题**：忽略了背面剔除的边界情况
- 透明物体需要双面测试
- 法线翻转的处理

### 4. 边界情况

**问题**：退化三角形（共线顶点）导致除零错误
- 预处理移除退化几何
- 相交测试中添加检查

**问题**：极小或极大的场景规模
- 使用相对误差而非绝对误差
- 考虑分层LOD

## 最佳实践检查清单

### 设计阶段
- [ ] 选择合适的加速结构（BVH vs KD-Tree）
- [ ] 确定精度要求和误差容限
- [ ] 规划内存布局和缓存优化策略
- [ ] 考虑目标平台（CPU/GPU）特性

### 实现阶段
- [ ] 光线方向始终保持归一化
- [ ] 正确处理数值精度问题
- [ ] 实现鲁棒的相交测试算法
- [ ] 添加适当的debug可视化

### 优化阶段
- [ ] 使用性能分析确定瓶颈
- [ ] 实现SIMD加速的相交测试
- [ ] 优化内存访问模式
- [ ] 考虑光线包（Ray Packet）技术

### 验证阶段
- [ ] 测试各种退化情况
- [ ] 验证数值稳定性
- [ ] 比较不同场景下的性能
- [ ] 确保结果的确定性（用于调试）

### 维护阶段
- [ ] 文档化关键算法选择
- [ ] 提供性能调优参数
- [ ] 支持渐进式渲染
- [ ] 保持代码的可扩展性