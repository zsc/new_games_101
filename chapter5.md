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

重心坐标的计算可以通过面积比得到：

$$u = \frac{\text{Area}(PBC)}{\text{Area}(ABC)}, \quad v = \frac{\text{Area}(APC)}{\text{Area}(ABC)}, \quad w = \frac{\text{Area}(ABP)}{\text{Area}(ABC)}$$

使用叉积计算：

$$u = \frac{(B-P) \times (C-P) \cdot \mathbf{n}}{(B-A) \times (C-A) \cdot \mathbf{n}}$$

其中$\mathbf{n}$是三角形法向量。

### 5.1.2 透视正确插值

屏幕空间的线性插值在透视投影下会产生错误。考虑属性$\phi$在三角形顶点的值为$\phi_0, \phi_1, \phi_2$，正确的插值公式为：

$$\phi = \frac{u\frac{\phi_0}{z_0} + v\frac{\phi_1}{z_1} + w\frac{\phi_2}{z_2}}{u\frac{1}{z_0} + v\frac{1}{z_1} + w\frac{1}{z_2}}$$

其中$z_i$是顶点的视空间深度。

推导过程：在投影变换下，属性$\phi/z$是线性的，因此：

$$\frac{\phi}{z} = u\frac{\phi_0}{z_0} + v\frac{\phi_1}{z_1} + w\frac{\phi_2}{z_2}$$

$$\frac{1}{z} = u\frac{1}{z_0} + v\frac{1}{z_1} + w\frac{1}{z_2}$$

### 5.1.3 属性插值的硬件实现

现代GPU使用专门的插值器硬件。常见的插值模式包括：

1. **flat shading**: 使用provoking vertex的属性值
2. **smooth shading**: 透视正确插值
3. **noperspective**: 屏幕空间线性插值

插值器的优化策略：
- 增量计算：利用$\Delta u, \Delta v$沿扫描线递增
- SIMD并行：同时插值多个片元
- 缓存优化：重用顶点属性

### 5.1.4 高阶插值

对于曲面片元，可能需要高阶插值。二次插值示例：

$$\phi(u,v) = \sum_{i+j \leq 2} a_{ij} u^i v^j$$

系数通过解线性系统确定。

## 5.2 高级纹理映射

### 5.2.1 纹理坐标的生成与变换

纹理坐标生成方法：

1. **参数化映射**：通过UV展开
2. **投影映射**：平面、圆柱、球面投影
3. **环境映射**：立方体贴图坐标计算

球面映射的数学表达：
$$u = \frac{1}{2\pi} \arctan\left(\frac{y}{x}\right) + 0.5$$
$$v = \frac{1}{\pi} \arccos(z) $$

### 5.2.2 MIP映射原理

MIP映射解决纹理采样的走样问题。对于屏幕空间导数：

$$\frac{\partial u}{\partial x} = \frac{\partial u}{\partial s} \cdot \frac{\partial s}{\partial x}$$

其中$s$是屏幕坐标。LOD级别计算：

$$\text{LOD} = \log_2 \left( \max \left( \sqrt{\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial x}\right)^2}, \sqrt{\left(\frac{\partial u}{\partial y}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2} \right) \cdot \text{textureSize} \right)$$

### 5.2.3 各向异性过滤

各向异性过滤考虑纹理在不同方向的拉伸。EWA (Elliptical Weighted Average) 过滤器：

纹理空间的椭圆方程：
$$\mathbf{J} = \begin{pmatrix} \frac{\partial u}{\partial x} & \frac{\partial u}{\partial y} \\ \frac{\partial v}{\partial x} & \frac{\partial v}{\partial y} \end{pmatrix}$$

椭圆的形状矩阵：
$$\mathbf{M} = \mathbf{J}^T \mathbf{J}$$

### 5.2.4 纹理压缩技术

常见压缩格式的数学原理：

**BC1 (DXT1)**：4×4块压缩到64位
- 2个16位颜色端点
- 32位索引（每像素2位）

压缩率：$\frac{16 \times 24}{64} = 6:1$

**BC7**：高质量压缩，支持alpha通道
- 多种模式选择
- 更精细的颜色量化

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
