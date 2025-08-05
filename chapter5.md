# 第5章：高级着色技术

## 章节概述

本章深入探讨计算机图形学中的高级着色技术。我们将从插值技术的数学基础开始，深入理解GPU如何在三角形内部进行属性插值。随后，我们将探索高级纹理映射技术，包括MIP映射、各向异性过滤等现代GPU技术。最后，我们将学习程序化纹理生成和着色器优化技术，为实现高性能渲染打下基础。

## 学习目标

- 掌握透视正确插值的数学原理
- 理解高级纹理映射技术的实现细节
- 学会设计和优化程序化纹理
- 掌握着色器性能优化的关键技术

## 5.1 插值技术详解

### 5.1.1 重心坐标系

在三角形内部进行插值是光栅化的核心操作。对于三角形ABC和内部点P，重心坐标$(u, v, w)$定义为：

$$P = u \cdot A + v \cdot B + w \cdot C$$

其中 $u + v + w = 1$。

#### 几何意义

重心坐标具有直观的几何解释：
- $u$表示点P到边BC的"相对距离"
- 当$u=1$时，P位于顶点A；当$u=0$时，P位于边BC上
- 重心坐标在仿射变换下保持不变

#### 计算方法

**方法1：面积比**

$$u = \frac{\text{Area}(PBC)}{\text{Area}(ABC)}, \quad v = \frac{\text{Area}(APC)}{\text{Area}(ABC)}, \quad w = \frac{\text{Area}(ABP)}{\text{Area}(ABC)}$$

**方法2：叉积形式**

对于2D情况：
$$u = \frac{(B-P) \times (C-P)}{(B-A) \times (C-A)}$$

对于3D情况，需要投影到合适的平面：
$$u = \frac{|(B-P) \times (C-P) \cdot \mathbf{n}|}{|(B-A) \times (C-A) \cdot \mathbf{n}|}$$

其中$\mathbf{n}$是三角形法向量。

**方法3：直接求解线性系统**

设$P = (p_x, p_y)$，可以建立方程：
$$\begin{pmatrix} a_x & b_x & c_x \\ a_y & b_y & c_y \\ 1 & 1 & 1 \end{pmatrix} \begin{pmatrix} u \\ v \\ w \end{pmatrix} = \begin{pmatrix} p_x \\ p_y \\ 1 \end{pmatrix}$$

使用Cramer法则求解：
$$u = \frac{\begin{vmatrix} p_x & b_x & c_x \\ p_y & b_y & c_y \\ 1 & 1 & 1 \end{vmatrix}}{\begin{vmatrix} a_x & b_x & c_x \\ a_y & b_y & c_y \\ 1 & 1 & 1 \end{vmatrix}}$$

#### 数值稳定性

在实际实现中，需要考虑：
1. **退化情况**：当三角形接近共线时，分母接近0
2. **精度问题**：使用double精度计算中间结果
3. **边界处理**：$u, v, w \in [0, 1]$的严格判断

#### 优化实现

对于光栅化中的大量计算，可以使用增量算法：
$$u(x+1, y) = u(x, y) + \Delta u_x$$
$$u(x, y+1) = u(x, y) + \Delta u_y$$

其中$\Delta u_x$和$\Delta u_y$是常数，可预计算。

### 5.1.2 透视正确插值

屏幕空间的线性插值在透视投影下会产生错误。这是因为透视投影是非线性变换，破坏了属性的线性关系。

#### 问题的本质

考虑世界空间中的线段，其属性（如纹理坐标）线性变化。经过透视投影后：
1. 空间位置被非线性压缩
2. 但屏幕空间的光栅化假设线性插值
3. 导致属性插值出现扭曲

#### 数学推导

设顶点属性为$\phi_0, \phi_1, \phi_2$，视空间深度为$z_0, z_1, z_2$。

**关键观察**：在齐次坐标系中，属性$\phi/w$（其中$w=z$）保持线性。

**推导步骤**：

1. 齐次空间的线性插值：
   $$\left(\frac{\phi}{z}\right)_p = u\frac{\phi_0}{z_0} + v\frac{\phi_1}{z_1} + w\frac{\phi_2}{z_2}$$

2. 深度的倒数也线性插值：
   $$\frac{1}{z_p} = u\frac{1}{z_0} + v\frac{1}{z_1} + w\frac{1}{z_2}$$

3. 恢复原始属性：
   $$\phi_p = \frac{\left(\frac{\phi}{z}\right)_p}{\frac{1}{z_p}} = \frac{u\frac{\phi_0}{z_0} + v\frac{\phi_1}{z_1} + w\frac{\phi_2}{z_2}}{u\frac{1}{z_0} + v\frac{1}{z_1} + w\frac{1}{z_2}}$$

#### 实现细节

**顶点着色器输出**：
```glsl
// 输出属性除以w
out_texcoord = texcoord / gl_Position.w;
out_inv_w = 1.0 / gl_Position.w;
```

**片元着色器恢复**：
```glsl
// 硬件自动插值后恢复
texcoord = in_texcoord / in_inv_w;
```

#### 性能优化

1. **预计算优化**：在顶点着色器计算$1/z$，避免片元着色器除法
2. **向量化**：同时处理多个属性的透视校正
3. **精度控制**：关键计算使用highp精度

#### 特殊情况

**屏幕空间属性**：某些属性（如屏幕空间导数）不需要透视校正，使用`noperspective`限定符：
```glsl
noperspective out vec2 screen_coord;
```

**常见错误案例**：
- 法线插值：需要先归一化再插值
- 颜色插值：线性颜色空间vs sRGB空间
- 切线空间：需要保持正交性

### 5.1.3 属性插值的硬件实现

现代GPU使用专门的插值器硬件，这是固定功能管线的关键组件之一。

#### 插值模式详解

1. **flat shading**: 
   - 使用provoking vertex（通常是最后一个顶点）的属性值
   - 整个三角形使用同一值，无插值计算
   - 适用于：图元ID、材质索引等离散属性

2. **smooth shading**: 
   - 执行完整的透视正确插值
   - 默认模式，适用于大多数属性
   - 硬件自动处理$1/w$的计算和恢复

3. **noperspective**: 
   - 屏幕空间线性插值，跳过透视校正
   - 适用于：屏幕空间效果、后处理坐标
   - 性能更高但可能产生透视失真

4. **centroid sampling**:
   - 在MSAA情况下，确保采样点在三角形内部
   - 避免边缘artifacts但可能影响导数计算

#### 硬件架构

**参数缓存**：
- 存储三角形的插值参数（梯度、起始值）
- 典型大小：32-64个三角形
- LRU替换策略

**插值管线**：
```
顶点属性 → 参数计算 → 参数缓存 → 片元插值 → 属性输出
           ↓                      ↑
      梯度计算单元            扫描线生成器
```

#### 增量算法实现

对于属性$\phi$在屏幕坐标$(x,y)$的值：

$$\phi(x,y) = \phi_0 + \frac{\partial \phi}{\partial x}(x-x_0) + \frac{\partial \phi}{\partial y}(y-y_0)$$

梯度计算：
$$\frac{\partial \phi}{\partial x} = \frac{(\phi_1 - \phi_0)(y_2 - y_0) - (\phi_2 - \phi_0)(y_1 - y_0)}{(x_1 - x_0)(y_2 - y_0) - (x_2 - x_0)(y_1 - y_0)}$$

**扫描线遍历优化**：
```
// 沿扫描线的增量
for (x = x_start; x <= x_end; x++) {
    phi += dphi_dx;
    inv_w += dinv_w_dx;
}
// 换行时的调整
phi_row += dphi_dy;
inv_w_row += dinv_w_dy;
```

#### SIMD并行策略

**2×2 Quad并行**：
- GPU以2×2像素块（quad）为单位处理
- 便于计算屏幕空间导数
- 共享插值参数，减少带宽

**导数计算**：
$$\frac{\partial \phi}{\partial x} \approx \phi_{x+1,y} - \phi_{x,y}$$
$$\frac{\partial \phi}{\partial y} \approx \phi_{x,y+1} - \phi_{x,y}$$

#### 精度与性能权衡

**定点数优化**：
- 屏幕坐标使用定点数（如16.8格式）
- 避免浮点累积误差
- 提高遍历精度

**属性压缩**：
- 16位浮点用于颜色、纹理坐标
- 32位浮点用于位置、深度
- 打包多个小属性到单个寄存器

**缓存优化策略**：
1. **顶点缓存**：避免重复变换
2. **参数缓存**：重用三角形设置
3. **纹理缓存**：利用2D局部性

### 5.1.4 高阶插值

线性插值可能不足以表示复杂的属性变化，特别是在曲面渲染和高质量着色中。

#### 二次插值

**Bézier三角形表示**：
对于二次Bézier三角形，使用6个控制点：
$$\phi(u,v,w) = \sum_{i+j+k=2} \binom{2}{i,j,k} u^i v^j w^k \phi_{ijk}$$

其中$\binom{2}{i,j,k} = \frac{2!}{i!j!k!}$是多项式系数。

**展开形式**：
$$\phi(u,v,w) = u^2\phi_{200} + v^2\phi_{020} + w^2\phi_{002} + 2uv\phi_{110} + 2vw\phi_{011} + 2uw\phi_{101}$$

**控制点布局**：
```
      φ200
     /    \
  φ110    φ101
   /        \
φ020--φ011--φ002
```

#### 三次插值

**PN三角形（Point-Normal Triangles）**：
用于几何细分和光滑表面：

1. **几何控制点**：
   $$b_{ijk} = \frac{i\mathbf{P}_1 + j\mathbf{P}_2 + k\mathbf{P}_3}{3} + \frac{ij}{3}\mathbf{w}_{12} + \frac{jk}{3}\mathbf{w}_{23} + \frac{ki}{3}\mathbf{w}_{31}$$
   
   其中$\mathbf{w}_{ij}$是修正向量。

2. **法线插值**：
   使用二次插值保证$C^1$连续性

#### 有理插值

**NURBS三角形**：
$$\phi(u,v,w) = \frac{\sum_{i+j+k=n} w_{ijk} \phi_{ijk} B_{ijk}^n(u,v,w)}{\sum_{i+j+k=n} w_{ijk} B_{ijk}^n(u,v,w)}$$

其中$w_{ijk}$是权重，$B_{ijk}^n$是Bernstein基函数。

#### 自适应插值

**误差驱动细分**：
根据插值误差动态选择插值阶数：

$$E = \max_{\mathbf{p} \in T} |\phi_{linear}(\mathbf{p}) - \phi_{exact}(\mathbf{p})|$$

当$E > \epsilon$时，使用高阶插值或细分三角形。

#### 实现考虑

**存储开销**：
- 线性：3个系数
- 二次：6个系数
- 三次：10个系数

**计算复杂度**：
- 线性：3次乘法，2次加法
- 二次：6次乘法，5次加法，3次平方
- 三次：10次乘法，9次加法，需要立方计算

**硬件支持**：
- 现代GPU主要支持线性插值
- 高阶插值通过着色器实现
- 曲面细分着色器提供硬件加速

**应用场景**：
1. **位移贴图**：需要平滑的几何变化
2. **法线插值**：避免Gouraud着色的棱角
3. **程序化纹理**：平滑的参数变化

## 5.2 高级纹理映射

### 5.2.1 纹理坐标的生成与变换

纹理坐标是将二维图像映射到三维几何体的桥梁。不同的生成方法适用于不同的几何形状和应用场景。

#### 参数化映射

**UV展开算法**：

1. **保角映射（Conformal Mapping）**：
   - 保持局部角度，适合纹理绘制
   - 使用复变函数理论：$f: \mathbb{C} \to \mathbb{C}$
   - 满足Cauchy-Riemann方程：$\frac{\partial u}{\partial x} = \frac{\partial v}{\partial y}$, $\frac{\partial u}{\partial y} = -\frac{\partial v}{\partial x}$

2. **保面积映射（Area-Preserving）**：
   - 保持三角形面积比
   - 优化目标：$\min \sum_T \left(\frac{A_{3D}}{A_{2D}} - 1\right)^2$

3. **LSCM（Least Squares Conformal Maps）**：
   - 最小二乘保角映射
   - 线性系统求解，适合大规模网格

#### 投影映射

**平面投影**：
```
u = dot(position - origin, uAxis) / uScale
v = dot(position - origin, vAxis) / vScale
```

**圆柱投影**：
$$u = \frac{1}{2\pi} \arctan2(y, x) + 0.5$$
$$v = \frac{z - z_{min}}{z_{max} - z_{min}}$$

**球面投影**：
$$u = \frac{1}{2\pi} \arctan2(y, x) + 0.5$$
$$v = \frac{1}{\pi} \arccos\left(\frac{z}{|\mathbf{r}|}\right)$$

注意：在极点处存在奇异性，需要特殊处理。

#### 环境映射

**立方体贴图（Cube Mapping）**：

1. **面选择**：
   ```
   absX = |direction.x|
   absY = |direction.y|
   absZ = |direction.z|
   
   if (absX >= absY && absX >= absZ) {
       face = (direction.x > 0) ? +X : -X
   }
   // 类似处理Y和Z
   ```

2. **坐标计算**：
   对于+X面：
   $$u = 0.5 \times \left(1 + \frac{-direction.z}{direction.x}\right)$$
   $$v = 0.5 \times \left(1 + \frac{-direction.y}{direction.x}\right)$$

**球面环境贴图**：
- Equirectangular：直接使用球面坐标
- Dual-Paraboloid：两个抛物面覆盖全球

#### 动态纹理坐标生成

**三平面映射（Triplanar Mapping）**：
```glsl
vec3 blendWeights = abs(normal);
blendWeights = normalize(max(blendWeights, 0.00001));

vec3 xaxis = texture(tex, position.yz) * blendWeights.x;
vec3 yaxis = texture(tex, position.xz) * blendWeights.y;
vec3 zaxis = texture(tex, position.xy) * blendWeights.z;

return xaxis + yaxis + zaxis;
```

**投影纹理**：
从光源视角投影：
$$\mathbf{uv} = \frac{1}{2} + \frac{1}{2} \times \frac{\mathbf{P}_{light}}{w_{light}}$$

#### 纹理坐标变换

**仿射变换**：
$$\begin{pmatrix} u' \\ v' \end{pmatrix} = \begin{pmatrix} a & b \\ c & d \end{pmatrix} \begin{pmatrix} u \\ v \end{pmatrix} + \begin{pmatrix} t_x \\ t_y \end{pmatrix}$$

**透视变换**：
$$u' = \frac{au + bv + c}{gu + hv + 1}$$
$$v' = \frac{du + ev + f}{gu + hv + 1}$$

### 5.2.2 MIP映射原理

MIP（Multum In Parvo，“多在小中”）映射是解决纹理走样的核心技术。它通过预计算不同分辨率的纹理金字塔，在渲染时选择合适的级别。

#### 走样问题分析

**欠采样（Undersampling）**：
- 当纹理频率超过Nyquist频率时发生
- 产生摩尔纹、闪烁等artifacts
- 根本原因：屏幕像素覆盖多个纹素

**过采样（Oversampling）**：
- 浪费带宽和计算资源
- 可能导致纹理缓存效率低下

#### MIP级别计算

**屏幕空间导数**：
使用雾可比矩阵计算纹理坐标的变化率：

$$\mathbf{J} = \begin{pmatrix} \frac{\partial u}{\partial x} & \frac{\partial u}{\partial y} \\ \frac{\partial v}{\partial x} & \frac{\partial v}{\partial y} \end{pmatrix}$$

**像素覆盖区域估计**：
1. **各向同性估计**（标准MIP）：
   $$\rho = \max\left(||\mathbf{J} \cdot (1, 0)^T||, ||\mathbf{J} \cdot (0, 1)^T||\right)$$

2. **完整LOD计算**：
   $$\text{LOD} = \log_2(\rho \cdot \text{textureSize})$$

3. **LOD偏移与约束**：
   $$\text{LOD}_{final} = \text{clamp}(\text{LOD} + \text{bias}, \text{minLOD}, \text{maxLOD})$$

#### MIP金字塔生成

**滤波核选择**：

1. **简单平均（Box Filter）**：
   $$\text{MIP}_{i+1}(x,y) = \frac{1}{4}\sum_{dx,dy \in \{0,1\}} \text{MIP}_i(2x+dx, 2y+dy)$$

2. **高斯滤波**：
   $$w(x,y) = \frac{1}{2\pi\sigma^2} e^{-\frac{x^2+y^2}{2\sigma^2}}$$

3. **Lanczos滤波**：
   $$L(x) = \begin{cases}
   \text{sinc}(x) \cdot \text{sinc}(x/a) & |x| < a \\
   0 & \text{otherwise}
   \end{cases}$$

**边界处理**：
- Clamp: 防止越界采样
- Wrap: 周期性纹理
- Mirror: 镜像反射

#### 三线性过滤

在两个MIP级别之间进行线性插值：

$$\text{color} = (1-\alpha) \cdot \text{sample}(\lfloor\text{LOD}\rfloor, u, v) + \alpha \cdot \text{sample}(\lceil\text{LOD}\rceil, u, v)$$

其中$\alpha = \text{fract}(\text{LOD})$。

#### 导数计算优化

**显式导数传递**：
```glsl
vec4 textureGrad(sampler2D tex, vec2 uv, vec2 dPdx, vec2 dPdy)
```

**自动导数计算**：
GPU在2×2 quad内使用差分：
$$\frac{\partial f}{\partial x} \approx f(x+1, y) - f(x, y)$$

#### 存储优化

**内存占用**：
完整MIP链增加约33%的存储：
$$\text{Total} = \text{Base} \times \sum_{i=0}^{\infty} \frac{1}{4^i} = \frac{4}{3} \times \text{Base}$$

**紧凑布局**：
```
+--------+----+--+-+
| Level0 | L1 |L2|3|
|        +----+--+-+
|        | L1 |L2|3|
+--------+----+--+-+
```

### 5.2.3 各向异性过滤

各向异性过滤解决了标准MIP映射在倾斜视角下的模糊问题。当纹理在不同方向上拉伸程度不同时，各向同性的MIP滤波会过度模糊。

#### 像素足迹分析

**像素在纹理空间的投影**：
屏幕像素在纹理空间形成一个近似椭圆的区域。

**雅可比矩阵**：
$$\mathbf{J} = \begin{pmatrix} \frac{\partial u}{\partial x} & \frac{\partial u}{\partial y} \\ \frac{\partial v}{\partial x} & \frac{\partial v}{\partial y} \end{pmatrix}$$

**形状矩阵**：
$$\mathbf{M} = \mathbf{J}^T \mathbf{J} = \begin{pmatrix} E & F \\ F & G \end{pmatrix}$$

其中：
- $E = \left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial x}\right)^2$
- $F = \frac{\partial u}{\partial x}\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\frac{\partial v}{\partial y}$
- $G = \left(\frac{\partial u}{\partial y}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2$

#### 椭圆参数计算

**特征值分解**：
$$\lambda_{1,2} = \frac{E + G \pm \sqrt{(E-G)^2 + 4F^2}}{2}$$

**椭圆主轴**：
- 长轴：$a = \sqrt{\lambda_{max}}$
- 短轴：$b = \sqrt{\lambda_{min}}$
- 旋转角：$\theta = \frac{1}{2}\arctan\left(\frac{2F}{E-G}\right)$

**各向异性比**：
$$r = \frac{a}{b} = \sqrt{\frac{\lambda_{max}}{\lambda_{min}}}$$

#### EWA过滤器

**椭圆加权平均（EWA）**：
$$c = \frac{\sum_{(s,t) \in E} w(s,t) \cdot \text{texture}(s,t)}{\sum_{(s,t) \in E} w(s,t)}$$

**高斯权重函数**：
$$w(s,t) = \exp\left(-\frac{1}{2}(s,t)\mathbf{M}^{-1}(s,t)^T\right)$$

#### 实现方法

**1. Ripmap**：
- 存储各种宽高比的预滤波纹理
- 内存开销：原始纹理的4倍
- 索引：$(\log_2(width), \log_2(height))$

**2. 多采样近似**：
```glsl
vec4 anisotropicSample(sampler2D tex, vec2 uv, vec2 ddx, vec2 ddy) {
    float maxAniso = 16.0;
    
    // 计算各向异性方向
    vec2 dx = ddx * textureSize(tex, 0);
    vec2 dy = ddy * textureSize(tex, 0);
    float px = dot(dx, dx);
    float py = dot(dy, dy);
    
    float maxLod = 0.5 * log2(max(px, py));
    float minLod = 0.5 * log2(min(px, py));
    
    float anisoRatio = min(2.0 * (maxLod - minLod), log2(maxAniso));
    float samples = exp2(anisoRatio);
    
    // 沿主轴采样
    vec2 majorAxis = (px > py) ? normalize(ddx) : normalize(ddy);
    vec4 color = vec4(0.0);
    
    for (float i = 0; i < samples; i++) {
        float t = (i + 0.5) / samples - 0.5;
        vec2 offset = t * majorAxis * length(ddx + ddy);
        color += textureLod(tex, uv + offset, minLod);
    }
    
    return color / samples;
}
```

**3. 硬件实现策略**：
- 限制最大采样数（通常8x或8x）
- 使用查找表加速三角函数
- 缓存友好的采样模式

#### 性能与质量权衡

**自适应采样**：
- 根据各向异性比动态调整采样数
- 远处物体使用较少采样
- 近处倾斜表面增加采样

**LOD偏移优化**：
$$\text{LOD}_{aniso} = \text{LOD}_{iso} - \frac{1}{2}\log_2(\text{samples})$$

### 5.2.4 纹理压缩技术

纹理压缩是平衡内存占用和视觉质量的关键技术。现代GPU在硬件级别支持实时解压。

#### 块压缩原理

**基本思想**：
- 将纹理分成4×4的块
- 每块独立压缩，便于随机访问
- 利用块内颜色相关性

#### BC1/DXT1格式

**数据结构**（每块64位）：
```
color0: 16 bits (R5G6B5)
color1: 16 bits (R5G6B5)
indices: 32 bits (2 bits × 16 pixels)
```

**颜色插值规则**：
当color0 > color1：
- $c_0 = color0$
- $c_1 = color1$
- $c_2 = \frac{2}{3}c_0 + \frac{1}{3}c_1$
- $c_3 = \frac{1}{3}c_0 + \frac{2}{3}c_1$

当color0 ≤ color1（透明模式）：
- $c_0 = color0$
- $c_1 = color1$
- $c_2 = \frac{1}{2}c_0 + \frac{1}{2}c_1$
- $c_3 = $ 透明/黑色

**压缩比**：$\frac{16 \times 24}{64} = 6:1$

#### BC7/BPTC格式

**高级特性**：
- 8种压缩模式
- 每块128位
- 支持完整Alpha通道
- 适合高质量贴图

**模式选择策略**：
```
Mode 0: 3 subsets, 4-bit endpoints
Mode 1: 2 subsets, 6-bit endpoints
Mode 2: 3 subsets, 5-bit endpoints
...
Mode 7: 2 subsets, 5-bit RGBA
```

**分区模式**：
将块分为多个子集，每个子集独立插值：
$$c = \sum_{i} w_i \cdot endpoint_i$$

#### ASTC格式

**灵活块大小**：
- 从4×4到8×8、甚至12×12
- 可变压缩率：8 bpp切0.89 bpp
- 支持HDR

**编码原理**：
1. **颜色端点编码**：使用BISE（Bounded Integer Sequence Encoding）
2. **权重编码**：量化到3、5、7位等
3. **分区模式**：支持1-4个分区

#### 纹理压缩选择策略

**根据内容类型**：

1. **漫反射贴图**：
   - BC1/DXT1：一般质量足够
   - BC7：高质量需求

2. **法线贴图**：
   - BC5：两通道，适合存储XY
   - 重建Z：$z = \sqrt{1 - x^2 - y^2}$

3. **HDR纹理**：
   - BC6H：半精度浮点
   - ASTC HDR模式

**性能考虑**：
```
缓存命中率 ∝ 1/压缩后大小
带宽需求 = 采样率 × 压缩后大小
```

#### 实时压缩优化

**快速编码算法**：
1. **PCA颜色选择**：使用主成分分析选择端点
2. **误差度量**：$E = \sum_{i} ||c_i - \hat{c}_i||^2$
3. **快速模式决策**：基于块的统计特性

**GPU加速编码**：
- 并行处理多个块
- 使用计算着色器
- 实时纹理流

### 5.2.5 虚拟纹理

虚拟纹理系统的页表管理：

$$\text{PhysicalAddr} = \text{PageTable}[\text{VirtualPage}] + \text{PageOffset}$$

缺页处理的优先级计算：
$$\text{Priority} = \frac{\text{ScreenCoverage}}{\text{LOD}^2}$$

## 5.3 程序化纹理与着色器优化

### 5.3.1 噪声函数

Perlin噪声的构造：

1. 梯度场生成：在整数格点随机分配梯度向量
2. 插值计算：
$$n(x,y,z) = \sum_{i,j,k \in \{0,1\}} w_i w_j w_k \cdot \mathbf{g}_{ijk} \cdot (x-x_{ijk}, y-y_{ijk}, z-z_{ijk})$$

其中$w_i = \text{fade}(x - \lfloor x \rfloor)$，fade函数：
$$\text{fade}(t) = 6t^5 - 15t^4 + 10t^3$$

### 5.3.2 分形纹理

fBm (Fractional Brownian Motion)：
$$f(x) = \sum_{i=0}^{n} \frac{\text{noise}(2^i x)}{2^i}$$

Turbulence函数：
$$t(x) = \sum_{i=0}^{n} \frac{|\text{noise}(2^i x)|}{2^i}$$

### 5.3.3 着色器优化技术

**算术强度优化**：
- 使用MAD指令：$a \times b + c$
- 向量化操作：同时处理RGBA
- 常量折叠：预计算静态表达式

**纹理带宽优化**：
- 纹理打包：将多个灰度纹理合并到RGBA通道
- 纹理图集：减少绑定切换
- 压缩格式选择：根据内容特性选择

**分支优化**：
动态分支的成本模型：
$$\text{Cost} = p \cdot C_{\text{true}} + (1-p) \cdot C_{\text{false}} + C_{\text{divergence}}$$

优化策略：
- 使用条件赋值代替分支
- 静态分支：使用uniform变量
- 动态分支合并：相似分支路径

### 5.3.4 计算着色器优化

**工作组大小优化**：
占用率计算：
$$\text{Occupancy} = \frac{\text{ActiveWarps}}{\text{MaxWarps}}$$

影响因素：
- 寄存器使用量
- 共享内存使用量
- 工作组大小

**内存访问模式**：
- Coalesced访问：连续线程访问连续地址
- Bank冲突避免：共享内存的步长访问
- 纹理缓存利用：2D局部性

## 本章小结

本章深入探讨了高级着色技术的核心概念：

### 关键公式回顾

1. **透视正确插值**：
   $$\phi = \frac{u\frac{\phi_0}{z_0} + v\frac{\phi_1}{z_1} + w\frac{\phi_2}{z_2}}{u\frac{1}{z_0} + v\frac{1}{z_1} + w\frac{1}{z_2}}$$

2. **MIP级别计算**：
   $$\text{LOD} = \log_2(\max(||\frac{\partial \mathbf{u}}{\partial x}||, ||\frac{\partial \mathbf{u}}{\partial y}||) \cdot \text{textureSize})$$

3. **Perlin噪声fade函数**：
   $$\text{fade}(t) = 6t^5 - 15t^4 + 10t^3$$

4. **fBm分形**：
   $$f(x) = \sum_{i=0}^{n} \frac{\text{noise}(2^i x)}{2^i}$$

### 核心概念

- **插值技术**：重心坐标、透视校正、硬件实现
- **纹理映射**：MIP映射、各向异性过滤、虚拟纹理
- **程序化纹理**：噪声函数、分形、优化技术
- **性能优化**：算术强度、内存访问、分支预测

## 练习题

### 基础题

**练习5.1：重心坐标计算**
给定三角形顶点$A(0,0)$，$B(4,0)$，$C(2,3)$和点$P(2,1)$，计算P的重心坐标。

*提示：使用面积法或叉积法计算。*

<details>
<summary>答案</summary>

使用叉积法：
- 三角形面积：$\text{Area}(ABC) = \frac{1}{2}|(B-A) \times (C-A)| = \frac{1}{2}|4 \times 3 - 0| = 6$
- $u = \frac{\text{Area}(PBC)}{\text{Area}(ABC)} = \frac{|(B-P) \times (C-P)|/2}{6} = \frac{|2 \times 2 - (-1) \times 0|/2}{6} = \frac{1}{3}$
- $v = \frac{\text{Area}(APC)}{\text{Area}(ABC)} = \frac{|(P-A) \times (C-A)|/2}{6} = \frac{|2 \times 3 - 1 \times 2|/2}{6} = \frac{1}{3}$
- $w = 1 - u - v = \frac{1}{3}$

因此重心坐标为$(1/3, 1/3, 1/3)$。
</details>

**练习5.2：透视插值误差**
考虑一条从$z=1$到$z=10$的线段，其纹理坐标从0到1线性变化。在屏幕空间$x=0.5$处，比较线性插值和透视正确插值的结果。

*提示：设置合适的投影矩阵参数。*

<details>
<summary>答案</summary>

假设透视投影后，$z=1$映射到$x'=0$，$z=10$映射到$x'=1$。

线性插值：$u_{\text{linear}} = 0.5$

透视正确插值：
- 在$x'=0.5$处，$\frac{1}{z} = 0.5 \times \frac{1}{1} + 0.5 \times \frac{1}{10} = 0.55$
- 因此$z = \frac{1}{0.55} \approx 1.82$
- $u = \frac{0.5 \times \frac{0}{1} + 0.5 \times \frac{1}{10}}{0.55} = \frac{0.05}{0.55} \approx 0.091$

误差：$|0.5 - 0.091| = 0.409$，说明线性插值会产生显著误差。
</details>

**练习5.3：MIP级别选择**
纹理大小为1024×1024，当前像素的纹理坐标导数为$\frac{\partial u}{\partial x} = 0.01$，$\frac{\partial v}{\partial x} = 0.02$，$\frac{\partial u}{\partial y} = 0.015$，$\frac{\partial v}{\partial y} = 0.01$。计算应选择的MIP级别。

*提示：使用最大各向异性导数。*

<details>
<summary>答案</summary>

计算各方向的导数大小：
- $||\frac{\partial \mathbf{u}}{\partial x}|| = \sqrt{0.01^2 + 0.02^2} = \sqrt{0.0005} \approx 0.0224$
- $||\frac{\partial \mathbf{u}}{\partial y}|| = \sqrt{0.015^2 + 0.01^2} = \sqrt{0.000325} \approx 0.0180$

最大导数：$\max(0.0224, 0.0180) = 0.0224$

LOD级别：
$$\text{LOD} = \log_2(0.0224 \times 1024) = \log_2(22.94) \approx 4.52$$

应选择MIP级别4或5（通常会在两者之间进行三线性插值）。
</details>

### 挑战题

**练习5.4：自定义插值基函数**
设计一个C²连续的插值基函数$h(t)$，满足$h(0)=0$，$h(1)=1$，$h'(0)=h'(1)=0$，$h''(0)=h''(1)=0$，并与Perlin的fade函数比较。

*提示：考虑7次多项式。*

<details>
<summary>答案</summary>

设$h(t) = at^7 + bt^6 + ct^5 + dt^4 + et^3 + ft^2 + gt + h$

边界条件：
- $h(0) = 0 \Rightarrow h = 0$
- $h(1) = 1 \Rightarrow a + b + c + d + e + f + g = 1$
- $h'(0) = 0 \Rightarrow g = 0$
- $h'(1) = 0 \Rightarrow 7a + 6b + 5c + 4d + 3e + 2f = 0$
- $h''(0) = 0 \Rightarrow f = 0$
- $h''(1) = 0 \Rightarrow 42a + 30b + 20c + 12d + 6e = 0$

解得：$h(t) = -20t^7 + 70t^6 - 84t^5 + 35t^4$

与Perlin的fade函数$f(t) = 6t^5 - 15t^4 + 10t^3$相比，新函数具有更高的连续性但计算成本更高。
</details>

**练习5.5：各向异性采样优化**
给定纹理空间的雅可比矩阵$\mathbf{J} = \begin{pmatrix} 2 & 1 \\ 0 & 0.5 \end{pmatrix}$，设计一个采样策略，使用最少的采样点覆盖95%的能量。

*提示：分析椭圆的主轴。*

<details>
<summary>答案</summary>

计算形状矩阵：
$$\mathbf{M} = \mathbf{J}^T\mathbf{J} = \begin{pmatrix} 4 & 2 \\ 2 & 1.25 \end{pmatrix}$$

特征值分解：
- 特征值：$\lambda_1 = 5.13$，$\lambda_2 = 0.12$
- 主轴长度：$a = \sqrt{\lambda_1} = 2.26$，$b = \sqrt{\lambda_2} = 0.35$
- 各向异性比：$\frac{a}{b} = 6.5$

采样策略：
1. 沿主轴方向采样$n = \lceil 2a \rceil = 5$个点
2. 垂直方向采样$m = \lceil 2b \rceil = 1$个点
3. 使用高斯权重：$w(x,y) = \exp(-\frac{x^2}{2a^2} - \frac{y^2}{2b^2})$

总共需要5个采样点即可覆盖95%的能量。
</details>

**练习5.6：虚拟纹理缓存策略**
设计一个虚拟纹理系统的页面替换算法，考虑：页面大小128×128，物理内存可存储256页，每帧预测未来3帧的访问模式。

*提示：结合LRU和预测信息。*

<details>
<summary>答案</summary>

混合替换算法设计：

1. **优先级计算**：
   $$P_i = \alpha \cdot \text{LRU}_i + \beta \cdot \text{Future}_i + \gamma \cdot \text{LOD}_i$$
   
   其中：
   - $\text{LRU}_i$：最近使用时间（归一化）
   - $\text{Future}_i$：未来3帧的预测访问概率
   - $\text{LOD}_i$：MIP级别权重（高分辨率优先）

2. **参数设置**：
   - $\alpha = 0.4$（历史权重）
   - $\beta = 0.5$（预测权重）
   - $\gamma = 0.1$（分辨率权重）

3. **实现细节**：
   - 维护访问历史的循环缓冲区
   - 使用运动向量预测未来访问
   - 预加载边界页面

4. **性能优化**：
   - 批量替换：一次替换多个低优先级页面
   - 异步加载：使用双缓冲避免停顿
   - 压缩传输：使用GPU解压
</details>

**练习5.7：着色器自动优化**
给定着色器代码片段，设计一个优化流程，目标是减少50%的纹理带宽。考虑纹理重用、计算与访存平衡等因素。

*提示：构建数据依赖图。*

<details>
<summary>答案</summary>

优化流程设计：

1. **静态分析**：
   - 构建SSA形式的中间表示
   - 识别纹理访问模式
   - 计算复用距离

2. **优化变换**：
   a) **纹理合并**：
      - 识别相同坐标的多次采样
      - 合并到单个float4采样
   
   b) **循环优化**：
      - 提升循环不变的纹理访问
      - 向量化相邻像素的访问
   
   c) **预计算**：
      - 将view-dependent计算移到顶点着色器
      - 使用查找表替代复杂函数

3. **性能模型**：
   $$\text{Bandwidth} = \sum_i \text{AccessCount}_i \times \text{CacheMiss}_i \times \text{DataSize}_i$$

4. **验证方法**：
   - 使用GPU性能计数器验证
   - A/B测试不同优化组合
   - 确保视觉质量不降低

预期结果：通过纹理打包和访问模式优化，可实现40-60%的带宽减少。
</details>

## 常见陷阱与错误

### 1. 插值精度问题
- **错误**：使用float16进行透视除法
- **正确**：透视除法必须使用float32以上精度
- **原因**：$1/z$的动态范围很大，低精度会导致严重误差

### 2. MIP链生成错误
- **错误**：简单的2×2平均
- **正确**：使用适当的滤波核（如Lanczos）
- **影响**：低质量MIP链导致摩尔纹和模糊

### 3. 纹理采样边界处理
- **错误**：忽略纹理边界的半像素偏移
- **正确**：考虑纹素中心对齐
- **公式**：$u' = \frac{u \times (size - 1) + 0.5}{size}$

### 4. 着色器精度损失
- **错误**：`normalize(normal + tangent)`
- **正确**：先normalize再插值，或使用四元数
- **原因**：向量加法破坏单位长度约束

### 5. 虚拟纹理反馈延迟
- **问题**：页面请求到加载完成有多帧延迟
- **解决**：预测性加载、多分辨率回退
- **权衡**：内存使用vs响应速度

## 最佳实践检查清单

### 插值实现检查
- [ ] 透视正确插值用于所有view-dependent属性
- [ ] 重心坐标计算考虑数值稳定性
- [ ] 正确处理退化三角形（面积接近0）
- [ ] 属性插值的精度匹配渲染需求

### 纹理系统设计
- [ ] MIP链生成使用高质量滤波
- [ ] 各向异性级别根据硬件能力设置
- [ ] 纹理格式选择平衡质量与带宽
- [ ] 虚拟纹理页面大小优化（通常128×128）
- [ ] 纹理图集考虑边界填充（防止渗色）

### 着色器优化验证
- [ ] 使用性能分析工具测量瓶颈
- [ ] 算术强度达到目标硬件的理想值
- [ ] 寄存器使用不超过硬件限制
- [ ] 分支发散度控制在可接受范围
- [ ] 纹理缓存命中率监控
- [ ] 考虑不同硬件的移植性

### 质量保证
- [ ] 各向异性过滤在掠射角下无明显瑕疵
- [ ] 程序化纹理在不同尺度下连续
- [ ] LOD切换无明显跳变
- [ ] 纹理压缩artifacts在可接受范围
- [ ] 性能优化未引入视觉退化
