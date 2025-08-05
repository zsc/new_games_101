# 第4章：着色基础

在前面的章节中，我们学习了如何将三维几何对象投影到二维屏幕并进行光栅化。然而，仅有几何形状是不够的——我们需要为这些形状赋予真实的外观。本章将介绍着色的基础知识，包括光照模型、着色频率以及纹理映射的基本概念。这些内容构成了计算机图形学中创建逼真图像的核心基础。

## 4.1 光照与基本着色模型

### 4.1.1 光的物理基础

光是电磁波，在图形学中我们主要关注可见光谱（波长约380-780纳米）。光与物体表面的相互作用决定了我们所看到的颜色和明暗。

**光的基本属性：**
- 强度（Intensity）：光的能量大小
- 方向（Direction）：光传播的方向
- 颜色（Color）：由光谱分布决定

### 4.1.2 Lambertian反射模型

最简单的反射模型是Lambertian（朗伯）反射，描述了完全漫反射表面的行为。

**朗伯余弦定律：**
$$I = I_0 \cos\theta = I_0 (\mathbf{n} \cdot \mathbf{l})$$

其中：
- $I$ 是反射光强度
- $I_0$ 是入射光强度
- $\theta$ 是表面法线与光线方向的夹角
- $\mathbf{n}$ 是归一化的表面法线
- $\mathbf{l}$ 是归一化的光线方向

### 4.1.3 Phong光照模型

Phong光照模型是计算机图形学中最经典的局部光照模型，由三个分量组成：

$$I = I_a + I_d + I_s$$

**1. 环境光（Ambient）：**
$$I_a = k_a I_{a,light}$$

**2. 漫反射（Diffuse）：**
$$I_d = k_d I_{light} \max(0, \mathbf{n} \cdot \mathbf{l})$$

**3. 镜面反射（Specular）：**
$$I_s = k_s I_{light} \max(0, \mathbf{r} \cdot \mathbf{v})^p$$

其中反射向量：
$$\mathbf{r} = 2(\mathbf{n} \cdot \mathbf{l})\mathbf{n} - \mathbf{l}$$

参数说明：
- $k_a, k_d, k_s$：材质的环境光、漫反射、镜面反射系数
- $p$：镜面反射指数（Shininess），控制高光的大小
- $\mathbf{v}$：观察方向

### 4.1.4 Blinn-Phong模型

Blinn-Phong是Phong模型的改进版本，使用半角向量代替反射向量：

$$\mathbf{h} = \frac{\mathbf{l} + \mathbf{v}}{||\mathbf{l} + \mathbf{v}||_2}$$

镜面反射项变为：
$$I_s = k_s I_{light} \max(0, \mathbf{n} \cdot \mathbf{h})^p$$

优点：
- 计算更高效（避免计算反射向量）
- 在某些情况下更符合物理规律
- 处理掠射角时表现更好

### 4.1.5 多光源处理

实际场景中通常有多个光源，总光照是所有光源贡献的叠加：

$$I_{total} = I_a + \sum_{i=1}^{n} (I_{d,i} + I_{s,i})$$

**光源类型：**
1. **点光源（Point Light）**：从一点发出，强度随距离衰减
   $$I(d) = \frac{I_0}{d^2}$$

2. **方向光（Directional Light）**：模拟无限远的光源（如太阳）
   - 所有位置的光线方向相同
   - 无衰减

3. **聚光灯（Spot Light）**：具有方向和角度限制的点光源
   $$I = I_0 \cdot \text{falloff}(\theta) \cdot \text{attenuation}(d)$$

## 4.2 着色频率与图形管线

### 4.2.1 着色频率的概念

着色频率决定了在渲染管线的哪个阶段计算光照：

**1. 平面着色（Flat Shading）**
- 每个三角形使用单一法线
- 计算一次光照，整个三角形使用相同颜色
- 优点：计算快速
- 缺点：产生明显的多边形边界

**2. Gouraud着色（顶点着色）**
- 在顶点处计算光照
- 三角形内部通过插值得到颜色
- 优点：比平面着色平滑
- 缺点：高光可能失真

**3. Phong着色（像素着色）**
- 插值法线而非颜色
- 在每个像素处计算光照
- 优点：高质量的光照效果
- 缺点：计算量大

### 4.2.2 法线插值

在Phong着色中，需要插值法线：

**重心坐标插值：**
$$\mathbf{n} = \alpha \mathbf{n}_1 + \beta \mathbf{n}_2 + \gamma \mathbf{n}_3$$

**注意：** 插值后需要重新归一化：
$$\mathbf{n}_{normalized} = \frac{\mathbf{n}}{||\mathbf{n}||_2}$$

### 4.2.3 透视校正插值

在透视投影下，屏幕空间的线性插值不等于3D空间的线性插值。

**透视校正公式：**
$$\frac{1}{z} = \alpha \frac{1}{z_1} + \beta \frac{1}{z_2} + \gamma \frac{1}{z_3}$$

对于属性$A$的插值：
$$A = \frac{\alpha \frac{A_1}{z_1} + \beta \frac{A_2}{z_2} + \gamma \frac{A_3}{z_3}}{\alpha \frac{1}{z_1} + \beta \frac{1}{z_2} + \gamma \frac{1}{z_3}}$$

### 4.2.4 现代图形管线中的着色

现代GPU使用可编程着色器：

**1. 顶点着色器（Vertex Shader）**
- 处理每个顶点
- 执行模型视图投影变换
- 计算顶点属性（可选的光照计算）

**2. 片段着色器（Fragment Shader）**
- 处理每个像素
- 执行光照计算（通常在这里）
- 纹理采样和最终颜色计算

**3. 延迟着色（Deferred Shading）**
- 第一遍：渲染几何信息到G-Buffer
- 第二遍：使用G-Buffer计算光照
- 优点：减少光照计算次数
- 缺点：内存带宽需求高，不支持透明物体

## 4.3 纹理映射基础

### 4.3.1 纹理映射的概念

纹理映射是将2D图像映射到3D表面的技术，用于增加表面细节而不增加几何复杂度。

**基本流程：**
1. 为3D模型的每个顶点指定纹理坐标$(u,v)$
2. 光栅化时插值纹理坐标
3. 使用插值后的坐标从纹理图像中采样颜色

### 4.3.2 纹理坐标系统

**UV坐标：**
- $u$：横向坐标，通常范围$[0,1]$
- $v$：纵向坐标，通常范围$[0,1]$
- $(0,0)$通常在左下角或左上角（取决于API）

**纹理坐标的获取方式：**
1. 手动指定（美术师在建模软件中设置）
2. 参数化映射（球面、圆柱、平面映射）
3. 自动展开（UV unwrapping算法）

### 4.3.3 纹理采样

**最近邻采样（Nearest Neighbor）：**
```
color = texture[round(u * width), round(v * height)]
```
- 优点：快速，保持像素化效果
- 缺点：产生锯齿

**双线性插值（Bilinear Interpolation）：**
对四个最近的纹素进行加权平均：

$$c = (1-s)(1-t)c_{00} + s(1-t)c_{10} + (1-s)t c_{01} + st c_{11}$$

其中$s = u \cdot width - \lfloor u \cdot width \rfloor$，$t$类似。

### 4.3.4 纹理放大与缩小

**放大（Magnification）问题：**
- 纹理分辨率低于屏幕采样率
- 解决方案：双线性或双三次插值

**缩小（Minification）问题：**
- 纹理分辨率高于屏幕采样率
- 导致走样和摩尔纹
- 解决方案：Mipmap

### 4.3.5 Mipmap技术

Mipmap是预计算的纹理金字塔，每一级分辨率减半：

**Mipmap级别选择：**
$$L = \log_2 \max\left(\sqrt{\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial x}\right)^2}, \sqrt{\left(\frac{\partial u}{\partial y}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2}\right)$$

**三线性插值（Trilinear Interpolation）：**
1. 在相邻两个mipmap级别分别进行双线性插值
2. 在两个结果之间进行线性插值

**各向异性过滤（Anisotropic Filtering）：**
- 处理非均匀缩放的情况
- 使用椭圆形采样区域而非正方形
- 提供更好的质量但计算成本更高

### 4.3.6 纹理应用

**1. 漫反射贴图（Diffuse Map）：**
$$k_d = \text{texture}(u, v)$$

**2. 法线贴图（Normal Map）：**
- 存储扰动后的法线信息
- 切线空间法线：RGB值映射到$[-1,1]$
$$\mathbf{n}_{tangent} = 2 \cdot \text{texture}(u,v) - 1$$

**3. 凹凸贴图（Bump Map）：**
- 存储高度信息
- 通过梯度计算法线扰动

**4. 环境贴图（Environment Map）：**
- 用于反射和环境光照
- 常见格式：立方体贴图、球面贴图

## 本章小结

本章介绍了计算机图形学中着色的基础知识：

**核心概念：**
1. **光照模型**：
   - Lambertian反射：$I = I_0 (\mathbf{n} \cdot \mathbf{l})$
   - Phong模型：$I = I_a + I_d + I_s$
   - Blinn-Phong改进：使用半角向量$\mathbf{h}$

2. **着色频率**：
   - 平面着色：每个三角形一个颜色
   - Gouraud着色：顶点计算，内部插值
   - Phong着色：像素级计算，插值法线

3. **纹理映射**：
   - UV坐标系统：$(u,v) \in [0,1]^2$
   - 采样方法：最近邻、双线性插值
   - Mipmap：解决纹理缩小问题
   - 各向异性过滤：处理非均匀缩放

**关键公式：**
- 反射向量：$\mathbf{r} = 2(\mathbf{n} \cdot \mathbf{l})\mathbf{n} - \mathbf{l}$
- 半角向量：$\mathbf{h} = \frac{\mathbf{l} + \mathbf{v}}{||\mathbf{l} + \mathbf{v}||_2}$
- 透视校正插值：$A = \frac{\sum \alpha_i \frac{A_i}{z_i}}{\sum \alpha_i \frac{1}{z_i}}$
- Mipmap级别：$L = \log_2 \max(||\frac{\partial \mathbf{u}}{\partial x}||, ||\frac{\partial \mathbf{u}}{\partial y}||)$

## 练习题

### 基础题

**练习4.1：Phong光照计算**
给定：
- 表面法线 $\mathbf{n} = (0, 1, 0)$
- 光线方向 $\mathbf{l} = \frac{1}{\sqrt{2}}(1, 1, 0)$
- 观察方向 $\mathbf{v} = (0, 0, 1)$
- 材质参数：$k_a = 0.1$, $k_d = 0.7$, $k_s = 0.2$, $p = 32$
- 光源强度：$I_{light} = 1.0$, $I_{a,light} = 0.2$

计算Phong模型下的最终颜色强度。

*Hint: 分别计算环境光、漫反射和镜面反射分量。*

<details>
<summary>答案</summary>

1. 环境光：$I_a = k_a I_{a,light} = 0.1 \times 0.2 = 0.02$

2. 漫反射：
   - $\mathbf{n} \cdot \mathbf{l} = (0,1,0) \cdot \frac{1}{\sqrt{2}}(1,1,0) = \frac{1}{\sqrt{2}}$
   - $I_d = k_d I_{light} (\mathbf{n} \cdot \mathbf{l}) = 0.7 \times 1.0 \times \frac{1}{\sqrt{2}} \approx 0.495$

3. 镜面反射：
   - $\mathbf{r} = 2(\mathbf{n} \cdot \mathbf{l})\mathbf{n} - \mathbf{l} = 2 \times \frac{1}{\sqrt{2}}(0,1,0) - \frac{1}{\sqrt{2}}(1,1,0) = \frac{1}{\sqrt{2}}(-1,1,0)$
   - $\mathbf{r} \cdot \mathbf{v} = \frac{1}{\sqrt{2}}(-1,1,0) \cdot (0,0,1) = 0$
   - $I_s = k_s I_{light} (\mathbf{r} \cdot \mathbf{v})^p = 0$

总强度：$I = I_a + I_d + I_s = 0.02 + 0.495 + 0 = 0.515$
</details>

**练习4.2：透视校正插值**
三角形顶点的屏幕坐标和深度值为：
- $P_1 = (0, 0)$, $z_1 = 2$, 纹理坐标 $u_1 = 0$
- $P_2 = (1, 0)$, $z_2 = 4$, 纹理坐标 $u_2 = 1$
- $P_3 = (0, 1)$, $z_3 = 4$, 纹理坐标 $u_3 = 0$

计算点$(0.5, 0.25)$处的透视校正纹理坐标。

*Hint: 先计算重心坐标，然后应用透视校正公式。*

<details>
<summary>答案</summary>

1. 重心坐标：$(0.5, 0.25)$对应$\alpha = 0.25$, $\beta = 0.5$, $\gamma = 0.25$

2. 透视校正：
   - $\frac{1}{z} = 0.25 \times \frac{1}{2} + 0.5 \times \frac{1}{4} + 0.25 \times \frac{1}{4} = 0.125 + 0.125 + 0.0625 = 0.3125$
   - $z = \frac{1}{0.3125} = 3.2$

3. 纹理坐标：
   - $u = \frac{0.25 \times \frac{0}{2} + 0.5 \times \frac{1}{4} + 0.25 \times \frac{0}{4}}{0.3125} = \frac{0.125}{0.3125} = 0.4$
</details>

**练习4.3：Mipmap级别计算**
在屏幕空间中，纹理坐标的偏导数为：
- $\frac{\partial u}{\partial x} = 0.01$, $\frac{\partial v}{\partial x} = 0.02$
- $\frac{\partial u}{\partial y} = 0.03$, $\frac{\partial v}{\partial y} = 0.01$

纹理大小为$1024 \times 1024$。计算应该使用的mipmap级别。

*Hint: 使用给定的mipmap级别公式。*

<details>
<summary>答案</summary>

1. 计算纹理空间中的偏导数（乘以纹理大小）：
   - $\frac{\partial s}{\partial x} = 0.01 \times 1024 = 10.24$
   - $\frac{\partial t}{\partial x} = 0.02 \times 1024 = 20.48$
   - $\frac{\partial s}{\partial y} = 0.03 \times 1024 = 30.72$
   - $\frac{\partial t}{\partial y} = 0.01 \times 1024 = 10.24$

2. 计算梯度长度：
   - $||\nabla_x|| = \sqrt{10.24^2 + 20.48^2} = \sqrt{524.29} \approx 22.9$
   - $||\nabla_y|| = \sqrt{30.72^2 + 10.24^2} = \sqrt{1048.58} \approx 32.4$

3. Mipmap级别：
   - $L = \log_2(32.4) \approx 5.02$
   - 使用级别5（或在5和6之间插值）
</details>

### 挑战题

**练习4.4：半球积分**
证明：对于完全漫反射的Lambertian表面，当入射光在半球上均匀分布时，反射的总能量为$\pi$倍的表面反射率。

*Hint: 在半球上对$\cos\theta$进行积分。*

<details>
<summary>答案</summary>

设表面反射率为$\rho$，入射辐射度为常数$L_i$。

反射的总能量：
$$E = \int_{\Omega} \rho L_i \cos\theta \, d\omega$$

在球坐标系中：
$$E = \rho L_i \int_0^{2\pi} \int_0^{\pi/2} \cos\theta \sin\theta \, d\theta \, d\phi$$

计算积分：
$$\int_0^{\pi/2} \cos\theta \sin\theta \, d\theta = \int_0^{\pi/2} \frac{1}{2}\sin(2\theta) \, d\theta = \frac{1}{2}$$

因此：
$$E = \rho L_i \times 2\pi \times \frac{1}{2} = \pi \rho L_i$$

证毕。
</details>

**练习4.5：法线贴图的切线空间**
给定三角形顶点：
- $P_1 = (0,0,0)$, $UV_1 = (0,0)$
- $P_2 = (1,0,0)$, $UV_2 = (1,0)$
- $P_3 = (0,1,1)$, $UV_3 = (0,1)$

计算切线空间的TBN矩阵。

*Hint: 使用边向量和UV差值构建线性方程组。*

<details>
<summary>答案</summary>

1. 计算边向量：
   - $\mathbf{e}_1 = P_2 - P_1 = (1,0,0)$
   - $\mathbf{e}_2 = P_3 - P_1 = (0,1,1)$
   - $\Delta u_1 = 1, \Delta v_1 = 0$
   - $\Delta u_2 = 0, \Delta v_2 = 1$

2. 切线和副切线满足：
   $$\begin{bmatrix} \mathbf{e}_1 \\ \mathbf{e}_2 \end{bmatrix} = \begin{bmatrix} \Delta u_1 & \Delta v_1 \\ \Delta u_2 & \Delta v_2 \end{bmatrix} \begin{bmatrix} \mathbf{T} \\ \mathbf{B} \end{bmatrix}$$

3. 求解得：
   - $\mathbf{T} = (1,0,0)$
   - $\mathbf{B} = (0,1,1)$

4. 法线：
   - $\mathbf{N} = \mathbf{T} \times \mathbf{B} = (0,-1,1)$
   - 归一化：$\mathbf{N} = \frac{1}{\sqrt{2}}(0,-1,1)$

5. TBN矩阵：
   $$TBN = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} \\ 0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \end{bmatrix}$$
</details>

**练习4.6：能量守恒约束**
对于物理真实的BRDF，证明Blinn-Phong模型不满足能量守恒。提出一种归一化方案。

*Hint: 考虑BRDF的半球积分应小于等于1。*

<details>
<summary>答案</summary>

Blinn-Phong的镜面反射项：
$$f_s = \frac{k_s}{4} (\mathbf{n} \cdot \mathbf{h})^p$$

半球积分：
$$\int_{\Omega} f_s \cos\theta_i \, d\omega_i$$

这个积分与$p$有关，当$p$较小时可能大于1，违反能量守恒。

归一化方案：
$$f_s = \frac{p + 2}{2\pi} k_s (\mathbf{n} \cdot \mathbf{h})^p$$

这确保了：
$$\int_{\Omega} f_s \cos\theta_i \, d\omega_i \leq k_s \leq 1$$
</details>

**练习4.7：各向异性纹理过滤**
推导各向异性过滤中椭圆参数的计算方法，给定纹理坐标的雅可比矩阵：
$$J = \begin{bmatrix} \frac{\partial u}{\partial x} & \frac{\partial u}{\partial y} \\ \frac{\partial v}{\partial x} & \frac{\partial v}{\partial y} \end{bmatrix}$$

*Hint: 考虑单位圆在纹理空间中的变形。*

<details>
<summary>答案</summary>

屏幕空间的单位圆映射到纹理空间形成椭圆。椭圆由$J^TJ$的特征值和特征向量决定。

1. 计算$M = J^TJ$：
   $$M = \begin{bmatrix} \left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial u}{\partial y}\right)^2 & \frac{\partial u}{\partial x}\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial v}{\partial y} \\ \frac{\partial u}{\partial x}\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\frac{\partial v}{\partial y} & \left(\frac{\partial v}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 \end{bmatrix}$$

2. 特征值$\lambda_1, \lambda_2$给出椭圆的长短轴：
   - 长轴：$a = \sqrt{\lambda_{\max}}$
   - 短轴：$b = \sqrt{\lambda_{\min}}$

3. 各向异性比：$r = \frac{a}{b} = \sqrt{\frac{\lambda_{\max}}{\lambda_{\min}}}$

4. 采样数量：$n = \min(\lceil r \rceil, \text{max\_anisotropy})$
</details>

**练习4.8：Phong插值的非线性效应**
考虑一个球面上的三角形，其三个顶点的法线分别指向外。分析在Phong着色下，三角形内部插值法线的长度变化，并讨论对光照计算的影响。

*Hint: 考虑法线插值后的长度不再为1。*

<details>
<summary>答案</summary>

设三个顶点的单位法线为$\mathbf{n}_1, \mathbf{n}_2, \mathbf{n}_3$。

插值法线：
$$\mathbf{n} = \alpha \mathbf{n}_1 + \beta \mathbf{n}_2 + \gamma \mathbf{n}_3$$

长度：
$$||\mathbf{n}||^2 = \alpha^2 + \beta^2 + \gamma^2 + 2\alpha\beta(\mathbf{n}_1 \cdot \mathbf{n}_2) + 2\beta\gamma(\mathbf{n}_2 \cdot \mathbf{n}_3) + 2\alpha\gamma(\mathbf{n}_1 \cdot \mathbf{n}_3)$$

由于$\mathbf{n}_i \cdot \mathbf{n}_j < 1$（顶点法线不平行），且$\alpha^2 + \beta^2 + \gamma^2 < 1$（除非在顶点），所以$||\mathbf{n}|| < 1$。

影响：
1. 未归一化的法线导致光照变暗
2. 三角形中心最暗（法线最短）
3. 必须在像素着色器中重新归一化

这解释了为什么Phong着色比Gouraud着色计算量大但效果更好。
</details>

## 常见陷阱与错误

### 1. 法线方向错误
- **问题**：法线指向内部导致光照反转
- **调试**：检查$\mathbf{n} \cdot \mathbf{l}$的符号，使用`max(0, dot)`
- **解决**：确保法线始终指向外部，考虑双面光照

### 2. 未归一化的向量
- **问题**：计算点积时使用未归一化的向量
- **症状**：光照强度异常，高光位置错误
- **解决**：在计算前始终归一化方向向量

### 3. 透视插值错误
- **问题**：直接在屏幕空间插值3D属性
- **症状**：纹理扭曲，特别是在大三角形上
- **解决**：使用透视校正插值公式

### 4. Mipmap边界伪影
- **问题**：mipmap级别突变导致可见边界
- **症状**：纹理出现明显的过渡线
- **解决**：使用三线性插值而非双线性

### 5. 伽马校正遗漏
- **问题**：在线性空间计算但显示时未校正
- **症状**：图像过暗或过亮
- **解决**：渲染到线性空间，显示时应用伽马校正

### 6. 纹理坐标越界
- **问题**：UV坐标超出$[0,1]$范围
- **症状**：纹理重复、镜像或夹紧
- **解决**：正确设置纹理包装模式（wrap/clamp/mirror）

## 最佳实践检查清单

### 光照计算
- [ ] 所有方向向量都已归一化
- [ ] 使用`max(0, dot)`避免负值光照
- [ ] 考虑了光源衰减（点光源）
- [ ] 正确处理背面（双面光照或剔除）
- [ ] 避免在顶点着色器计算高频光照效果

### 着色管线
- [ ] 根据需求选择合适的着色频率
- [ ] 使用透视校正插值
- [ ] 在片段着色器中重新归一化插值后的法线
- [ ] 考虑使用延迟着色优化多光源场景
- [ ] 正确设置和使用深度测试

### 纹理映射
- [ ] UV坐标在合理范围内
- [ ] 选择合适的过滤模式
- [ ] 为缩小情况生成mipmap
- [ ] 考虑各向异性过滤提升质量
- [ ] 纹理格式适合用途（sRGB/线性）
- [ ] 合理设置纹理包装模式

### 性能优化
- [ ] 避免在片段着色器中的分支
- [ ] 预计算可以在顶点着色器完成的内容
- [ ] 使用纹理缓存友好的访问模式
- [ ] 考虑LOD（细节层次）减少远处物体的计算
- [ ] 合理使用早期深度测试（Early-Z）
