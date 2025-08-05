# 第4章：着色基础

在前面的章节中，我们学习了如何将三维几何对象投影到二维屏幕并进行光栅化。然而，仅有几何形状是不够的——我们需要为这些形状赋予真实的外观。本章将介绍着色的基础知识，包括光照模型、着色频率以及纹理映射的基本概念。这些内容构成了计算机图形学中创建逼真图像的核心基础。

## 4.1 光照与基本着色模型

### 4.1.1 光的物理基础

光是电磁波，在图形学中我们主要关注可见光谱（波长约380-780纳米）。光与物体表面的相互作用决定了我们所看到的颜色和明暗。

**光的基本属性：**
- 强度（Intensity）：光的能量大小，通常用辐射度量学中的辐照度（Irradiance）或辐亮度（Radiance）表示
- 方向（Direction）：光传播的方向，在局部光照模型中假设为直线传播
- 颜色（Color）：由光谱分布决定，实践中常用RGB三通道近似
- 偏振（Polarization）：电磁波的振动方向，在高级渲染中用于表现某些材质特性

**辐射度量学基础：**
理解光照计算需要掌握几个关键的辐射度量学概念：

1. **辐射通量（Radiant Flux）**：$\Phi$，单位时间内通过表面的总能量，单位：瓦特(W)
   $$\Phi = \frac{dQ}{dt}$$

2. **辐照度（Irradiance）**：$E$，单位面积接收的辐射通量，单位：$W/m^2$
   $$E = \frac{d\Phi}{dA}$$

3. **辐射强度（Radiant Intensity）**：$I$，点光源在特定方向上的功率，单位：$W/sr$
   $$I = \frac{d\Phi}{d\omega}$$

4. **辐亮度（Radiance）**：$L$，最重要的量，描述光线的"亮度"，单位：$W/(m^2 \cdot sr)$
   $$L = \frac{d^2\Phi}{dA \cos\theta \, d\omega}$$

**光与表面的相互作用：**
当光线到达物体表面时，会发生以下现象：

1. **反射（Reflection）**：光线从表面弹回，分为镜面反射和漫反射
   - 镜面反射：遵循反射定律，入射角等于反射角
   - 漫反射：光线向各个方向均匀散射
   - 光泽反射：介于两者之间，有方向性但不完全镜面

2. **折射（Refraction）**：光线穿透表面进入物体内部，遵循斯涅尔定律
   $$n_1 \sin\theta_1 = n_2 \sin\theta_2$$
   其中$n_1, n_2$为介质的折射率

3. **吸收（Absorption）**：部分光能被材质吸收转化为热能
   - 遵循比尔-朗伯定律：$I(x) = I_0 e^{-\sigma_a x}$
   - $\sigma_a$为吸收系数，$x$为传播距离

4. **散射（Scattering）**：光线在介质内部发生多次反射和折射
   - 瑞利散射：粒子尺寸远小于波长（天空为蓝色的原因）
   - 米氏散射：粒子尺寸与波长相当（云和雾的外观）

**菲涅尔效应（Fresnel Effect）：**
反射率随入射角变化的现象，由菲涅尔方程描述：

对于非偏振光的反射率：
$$F(\theta) = \frac{1}{2}(F_\parallel^2 + F_\perp^2)$$

其中：
$$F_\parallel = \frac{n_2 \cos\theta_i - n_1 \cos\theta_t}{n_2 \cos\theta_i + n_1 \cos\theta_t}$$
$$F_\perp = \frac{n_1 \cos\theta_i - n_2 \cos\theta_t}{n_1 \cos\theta_i + n_2 \cos\theta_t}$$

Schlick近似（常用于实时渲染）：
$$F(\theta) = F_0 + (1 - F_0)(1 - \cos\theta)^5$$
其中$F_0$是垂直入射时的反射率。

在基本着色模型中，我们主要关注反射现象，将其简化为漫反射和镜面反射两个分量。

### 4.1.2 Lambertian反射模型

最简单的反射模型是Lambertian（朗伯）反射，描述了完全漫反射表面的行为。这种表面将入射光均匀地散射到各个方向，看起来同样明亮，与观察角度无关。

**朗伯余弦定律：**
$$I = I_0 \cos\theta = I_0 (\mathbf{n} \cdot \mathbf{l})$$

其中：
- $I$ 是反射光强度
- $I_0$ 是入射光强度  
- $\theta$ 是表面法线与光线方向的夹角
- $\mathbf{n}$ 是归一化的表面法线（指向表面外部）
- $\mathbf{l}$ 是归一化的光线方向（从表面指向光源）

**物理解释：**
余弦项$\cos\theta$反映了有效照射面积的变化。当光线垂直入射时（$\theta = 0$），单位面积接收最多能量；当光线平行于表面时（$\theta = 90°$），没有能量到达表面。

从微观角度理解：
- 粗糙表面由无数微小面元组成，朝向随机分布
- 每个微面元都是完美镜面，但整体表现为漫反射
- 余弦项反映了朝向光源的微面元比例

**BRDF形式：**
Lambertian反射的双向反射分布函数（BRDF）为常数：
$$f_r = \frac{\rho}{\pi}$$

其中$\rho$是反照率（albedo），表示表面反射的能量比例。除以$\pi$是为了满足能量守恒：
$$\int_{\Omega} f_r \cos\theta \, d\omega = \int_0^{2\pi} \int_0^{\pi/2} \frac{\rho}{\pi} \cos\theta \sin\theta \, d\theta \, d\phi = \rho$$

**漫反射的渲染方程：**
对于单个光源，出射辐亮度为：
$$L_o(\mathbf{x}, \omega_o) = f_r L_i(\mathbf{x}, \omega_i) (\mathbf{n} \cdot \omega_i)$$

展开后：
$$L_o = \frac{\rho}{\pi} L_i \max(0, \mathbf{n} \cdot \mathbf{l})$$

**实际应用考虑：**
- 需要clamp负值：$\max(0, \mathbf{n} \cdot \mathbf{l})$，避免背面照明
- 双面材质需要特殊处理：$|\mathbf{n} \cdot \mathbf{l}|$
- 可以乘以材质的漫反射系数$k_d$和颜色$\mathbf{c}_d$得到最终颜色

**高级扩展：**

1. **Oren-Nayar模型**：考虑表面粗糙度的改进模型
   $$f_r = \frac{\rho}{\pi}(A + B \max(0, \cos(\phi_i - \phi_o)) \sin\alpha \tan\beta)$$
   其中$A = 1 - 0.5\frac{\sigma^2}{\sigma^2 + 0.33}$，$B = 0.45\frac{\sigma^2}{\sigma^2 + 0.09}$，$\sigma$是表面粗糙度

2. **次表面散射（Subsurface Scattering）**：
   - 光线进入材质内部散射后射出
   - 常见于皮肤、蜡、大理石等半透明材质
   - 需要考虑光线入射点和出射点的距离

3. **预积分皮肤着色**：使用查找表（LUT）存储不同角度的散射结果

### 4.1.3 Phong光照模型

Phong光照模型是计算机图形学中最经典的局部光照模型，由Bui Tuong Phong于1975年提出。它通过三个分量的叠加来近似真实世界的光照效果：

$$I = I_a + I_d + I_s$$

**1. 环境光（Ambient）：**
$$I_a = k_a I_{a,light}$$

环境光模拟间接光照，是对全局光照的粗糙近似。它确保场景中没有完全黑暗的区域，提供基础亮度。

物理意义：
- 真实世界中，光线会在表面间多次反弹
- 环境光是这种复杂相互作用的简化
- 现代渲染中常用环境贴图或球谐函数代替

**2. 漫反射（Diffuse）：**
$$I_d = k_d I_{light} \max(0, \mathbf{n} \cdot \mathbf{l})$$

基于Lambertian反射模型，表示粗糙表面的反射特性。这是物体的主要颜色来源。

**3. 镜面反射（Specular）：**
$$I_s = k_s I_{light} \max(0, \mathbf{r} \cdot \mathbf{v})^p$$

模拟光滑表面的高光效果。反射向量$\mathbf{r}$通过镜面反射定律计算：
$$\mathbf{r} = 2(\mathbf{n} \cdot \mathbf{l})\mathbf{n} - \mathbf{l}$$

**推导反射向量公式：**
设入射向量为$-\mathbf{l}$，法线为$\mathbf{n}$：
1. 将$-\mathbf{l}$分解为平行和垂直于$\mathbf{n}$的分量
2. 平行分量：$(-\mathbf{l} \cdot \mathbf{n})\mathbf{n}$
3. 垂直分量：$-\mathbf{l} - (-\mathbf{l} \cdot \mathbf{n})\mathbf{n}$
4. 反射时，平行分量反向，垂直分量不变
5. $\mathbf{r} = -(-\mathbf{l} \cdot \mathbf{n})\mathbf{n} + [-\mathbf{l} - (-\mathbf{l} \cdot \mathbf{n})\mathbf{n}]$
6. 简化得：$\mathbf{r} = 2(\mathbf{n} \cdot \mathbf{l})\mathbf{n} - \mathbf{l}$

**参数说明：**
- $k_a, k_d, k_s$：材质系数，满足能量守恒时通常$k_a + k_d + k_s \leq 1$
- $p$：镜面反射指数（Shininess），典型值：
  - 粗糙表面：1-10
  - 塑料：10-100
  - 金属：100-1000
  - 镜面：1000+
- $\mathbf{v}$：归一化的观察方向（从表面指向相机）

**镜面反射的微面元理论解释：**
- 表面由许多微小镜面组成
- 只有法线恰好在$\mathbf{l}$和$\mathbf{v}$中间的微面元才能将光反射到观察者
- 指数$p$控制微面元法线分布的集中程度
- 高$p$值表示分布集中，产生小而亮的高光

**RGB颜色扩展：**
实际应用中，每个颜色通道独立计算：
$$\mathbf{I}_{RGB} = k_a \mathbf{c}_a \odot \mathbf{I}_{a,light} + k_d \mathbf{c}_d \odot \mathbf{I}_{light} (\mathbf{n} \cdot \mathbf{l}) + k_s \mathbf{c}_s \odot \mathbf{I}_{light} (\mathbf{r} \cdot \mathbf{v})^p$$

其中$\odot$表示逐元素乘法，$\mathbf{c}_a, \mathbf{c}_d, \mathbf{c}_s$分别是环境光、漫反射和镜面反射颜色。

**Phong模型的局限性：**
1. **非物理准确**：
   - 不遵循能量守恒
   - 镜面反射项不满足互易性（reciprocity）
   - 无法准确模拟粗糙表面的镜面反射

2. **掠射角问题**：
   - 在掠射角下，镜面反射强度下降过快
   - 真实材质在掠射角下反射率增加（菲涅尔效应）

3. **材质表现力有限**：
   - 无法表现菲涅尔效应
   - 对各向异性材质（如拉丝金属）无能为力
   - 不支持分层材质（如汽车漆）

4. **高光形状单一**：
   - 只能产生圆形高光
   - 无法表现拉长或其他形状的高光

### 4.1.4 Blinn-Phong模型

Blinn-Phong是Phong模型的改进版本，由Jim Blinn于1977年提出。核心改进是使用半角向量（halfway vector）代替反射向量：

$$\mathbf{h} = \frac{\mathbf{l} + \mathbf{v}}{||\mathbf{l} + \mathbf{v}||_2}$$

镜面反射项变为：
$$I_s = k_s I_{light} \max(0, \mathbf{n} \cdot \mathbf{h})^p$$

**物理解释：**
半角向量$\mathbf{h}$代表了微表面法线的统计分布。当$\mathbf{n} \cdot \mathbf{h}$较大时，意味着有更多的微表面朝向能够将光线反射到观察者。

**微面元理论深入：**
- 完美镜面反射发生在微面元法线等于$\mathbf{h}$时
- $(\mathbf{n} \cdot \mathbf{h})^p$近似微面元法线分布函数
- 这是后来Cook-Torrance BRDF的前身

**数学关系推导：**
设$\theta$为$\mathbf{n}$与$\mathbf{h}$的夹角，$\alpha$为$\mathbf{r}$与$\mathbf{v}$的夹角：
$$\cos\alpha = \mathbf{r} \cdot \mathbf{v} = 2\cos^2\theta - 1$$

因此：
$$(\mathbf{r} \cdot \mathbf{v})^{p_{Phong}} = (2\cos^2\theta - 1)^{p_{Phong}} \approx (\cos\theta)^{4p_{Phong}} = (\mathbf{n} \cdot \mathbf{h})^{4p_{Phong}}$$

这解释了为什么$p_{Blinn} \approx 4p_{Phong}$。

**优点：**
- 计算更高效（避免计算反射向量）
- 在某些情况下更符合物理规律
- 处理掠射角时表现更好
- 当光源和视点都在无穷远时，$\mathbf{h}$为常量

**实现细节：**
- 需要检查$\mathbf{l} + \mathbf{v}$是否为零向量（光源在视点后方）
- 半角向量的计算可以在顶点着色器中完成并插值
- 现代实时渲染几乎都使用Blinn-Phong而非原始Phong

**改进的归一化Blinn-Phong：**
为了更好的能量守恒，可以使用归一化版本：
$$I_s = \frac{p + 8}{8\pi} k_s I_{light} \max(0, \mathbf{n} \cdot \mathbf{h})^p$$

归一化因子确保镜面反射瓣的积分为$k_s$。

**各向异性扩展：**
对于各向异性材质，可以将标量指数$p$扩展为两个方向的指数：
$$I_s = k_s I_{light} \max(0, \mathbf{n} \cdot \mathbf{h})^{p_u \cos^2\phi + p_v \sin^2\phi}$$

其中$\phi$是半角向量在切平面上的投影与主切线方向的夹角。

### 4.1.5 多光源处理

实际场景中通常有多个光源，总光照是所有光源贡献的叠加：

$$I_{total} = I_a + \sum_{i=1}^{n} (I_{d,i} + I_{s,i})$$

注意环境光通常只计算一次，而非每个光源都有独立的环境光分量。

**光源类型：**

**1. 点光源（Point Light）**
从一点向所有方向发出光线，强度随距离衰减：
$$I(d) = \frac{I_0}{d^2}$$

**物理正确的衰减：**
- 遵循平方反比定律：能量在球面上均匀分布
- 球面面积$A = 4\pi d^2$，故单位面积能量$\propto 1/d^2$

实践中常用更灵活的衰减公式避免除零和过度衰减：
$$f_{att} = \frac{1}{k_c + k_l d + k_q d^2}$$
其中$k_c$（常数项）、$k_l$（线性项）、$k_q$（二次项）是衰减系数。

**改进的衰减模型：**
$$f_{att} = \frac{1}{\max(d, r_{min})^2}$$
其中$r_{min}$是光源的物理半径，避免近距离的无限亮度。

**2. 方向光（Directional Light）**
模拟无限远的光源（如太阳）：
- 所有位置的光线方向相同：$\mathbf{l} = \text{constant}$
- 无衰减：$f_{att} = 1$
- 适合室外场景的主光源
- 可产生平行阴影

**级联阴影贴图（Cascaded Shadow Maps）：**
- 将视锥体分成多个层级
- 每层使用不同分辨率的阴影贴图
- 近处高分辨率，远处低分辨率

**3. 聚光灯（Spot Light）**
具有位置、方向和照射角度的光源：
$$I = I_0 \cdot f_{spot}(\theta) \cdot f_{att}(d)$$

聚光衰减函数：
$$f_{spot}(\theta) = \begin{cases}
(\cos\theta)^{f} & \text{if } \cos\theta > \cos\theta_{outer} \\
\text{smooth}(\theta) & \text{if } \cos\theta_{outer} \geq \cos\theta > \cos\theta_{inner} \\
0 & \text{otherwise}
\end{cases}$$

**平滑过渡函数：**
使用Hermite插值实现软边缘：
$$\text{smooth}(\theta) = \text{smoothstep}(\cos\theta_{outer}, \cos\theta_{inner}, \cos\theta)$$
$$\text{smoothstep}(a, b, x) = t^2(3 - 2t), \quad t = \frac{x - a}{b - a}$$

**投影纹理聚光灯：**
- 使用纹理调制光强分布
- 可实现复杂图案投影（如舞台灯光）
- 纹理坐标从光源空间投影矩阵计算

**4. 面光源（Area Light）**
真实世界的光源都有一定面积，但精确计算需要积分。

**解析解（仅限简单几何）：**
矩形面光源对点的辐照度：
$$E = \int_A \frac{L \cos\theta_i \cos\theta_o}{r^2} dA$$

对于均匀发光的矩形，存在闭式解（使用边界积分）。

**近似方法：**
1. **多点近似**：
   - 将面光源离散为$N$个点光源
   - 权重根据面积元分配
   - 计算复杂度：$O(N)$

2. **线性变换球面光源（LTCs）**：
   - 将复杂BRDF变换为余弦分布
   - 使用预计算的查找表
   - 实时性能优秀

3. **球谐函数近似**：
   - 将光源投影到SH基函数
   - 适合低频光照
   - 系数预计算

**5. 环境光照（Environment Lighting）**
使用环境贴图表示来自所有方向的光照：

**重要性采样：**
- 根据环境贴图亮度分布采样
- 减少噪声，提高收敛速度
- 预计算CDF用于快速采样

**预滤波环境贴图：**
- 不同粗糙度级别的预滤波
- 存储在mipmap链中
- 实时查询，无需运行时积分

**性能优化策略：**

1. **光源剔除**：
   - 视锥剔除：丢弃视野外光源
   - 距离剔除：忽略超过最大影响距离的光源
   - 遮挡剔除：使用层次Z缓冲（Hi-Z）

2. **延迟着色（Deferred Shading）**：
   - 几何复杂度与光源数量解耦
   - G-Buffer存储几何属性
   - 光照以屏幕空间处理

3. **光源分级（Light Hierarchy）**：
   - 主光源：完整计算
   - 次要光源：简化BRDF
   - 环境光源：预积分近似

4. **光照烘焙（Light Baking）**：
   - 静态光源预计算
   - 存储在光照贴图或探针
   - 运行时仅处理动态光源

5. **Tiled/Clustered Shading**：
   - 将屏幕/视锥分块
   - 每块维护光源列表
   - 减少光源遍历开销

## 4.2 着色频率与图形管线

### 4.2.1 着色频率的概念

着色频率（Shading Frequency）决定了在渲染管线的哪个阶段计算光照。不同的着色频率在质量和性能之间做出不同的权衡。

**1. 平面着色（Flat Shading）**
- 每个三角形使用单一法线（通常为面法线）
- 计算一次光照，整个三角形使用相同颜色
- 优点：计算快速，适合低多边形风格或远距离LOD
- 缺点：产生明显的多边形边界（马赫带）
- 应用：CAD软件、低多边形艺术风格

**2. Gouraud着色（顶点着色）**
Henri Gouraud于1971年提出：
- 在顶点处计算光照（使用顶点法线）
- 三角形内部通过插值颜色得到最终结果
- 顶点法线通常是相邻面法线的平均
- 优点：比平面着色平滑，计算效率较高
- 缺点：
  - 高光可能失真（顶点没有高光时整个三角形都没有）
  - 不能正确处理非线性光照效果

**3. Phong着色（像素着色/片段着色）**
Bui Tuong Phong同时提出：
- 插值法线而非颜色
- 在每个像素处计算光照
- 优点：
  - 高质量的光照效果
  - 正确处理高光和非线性效果
  - 支持逐像素的法线贴图
- 缺点：计算量大（现代GPU已解决）

**选择标准：**
- 低多边形模型 + 高质量着色 = 错误的选择
- 高多边形模型 + Gouraud着色 = 可能足够
- 现代实践：几乎总是使用Phong着色（片段着色器）

### 4.2.2 法线插值

在Phong着色中，需要插值法线而非颜色。这是实现高质量光照的关键。

**重心坐标插值：**
$$\mathbf{n} = \alpha \mathbf{n}_1 + \beta \mathbf{n}_2 + \gamma \mathbf{n}_3$$

**为什么需要重新归一化：**
线性插值不保持向量长度。考虑两个相互垂直的单位法线$\mathbf{n}_1 = (1,0,0)$和$\mathbf{n}_2 = (0,1,0)$，它们的中点插值：
$$\mathbf{n}_{mid} = 0.5\mathbf{n}_1 + 0.5\mathbf{n}_2 = (0.5, 0.5, 0) \Rightarrow ||\mathbf{n}_{mid}|| = \frac{1}{\sqrt{2}} < 1$$

因此必须重新归一化：
$$\mathbf{n}_{normalized} = \frac{\mathbf{n}}{||\mathbf{n}||_2}$$

**法线插值的几何意义：**
- 线性插值在切平面上进行
- 归一化将结果投影回单位球面
- 结果是球面上的最短路径（大圆弧）的近似

**优化技巧：**
1. **快速归一化近似：**
   $$\mathbf{n}_{approx} \approx \mathbf{n} \cdot (1.5 - 0.5 \cdot \mathbf{n} \cdot \mathbf{n})$$
   基于牛顿-拉夫逊迭代，当$||\mathbf{n}||$接近1时非常准确

2. **四元数插值：**
   对于大角度旋转，可以使用四元数表示法线方向，使用SLERP插值

### 4.2.3 透视校正插值

在透视投影下，屏幕空间的线性插值不等于3D空间的线性插值。这是因为透视投影是非线性变换。

**问题的根源：**
考虑一条空间3D线段，其端点在不同深度。透视投影后：
- 3D空间的中点不映射到屏幕空间的中点
- 远处的部分被压缩，近处的部分被拉伸

**透视校正公式推导：**
透视投影中，屏幕坐标与3D坐标的关系：
$$x_s = \frac{x}{z}, \quad y_s = \frac{y}{z}$$

对于属性$A$，正确的插值是：
$$\frac{1}{z} = \alpha \frac{1}{z_1} + \beta \frac{1}{z_2} + \gamma \frac{1}{z_3}$$

$$A = \frac{\alpha \frac{A_1}{z_1} + \beta \frac{A_2}{z_2} + \gamma \frac{A_3}{z_3}}{\alpha \frac{1}{z_1} + \beta \frac{1}{z_2} + \gamma \frac{1}{z_3}}$$

**实现细节：**
1. 现代GPU自动处理透视校正
2. 需要存储$1/z$值（通常在深度缓冲中）
3. 所有顶点属性都需要透视校正：
   - 纹理坐标
   - 法线向量
   - 颜色值
   - 任何自定义属性

**没有透视校正的后果：**
- 纹理扭曲（特别是地面、墙壁）
- 光照不正确
- 动画中出现“游泳”效果

### 4.2.4 现代图形管线中的着色

现代GPU使用可编程着色器，提供了灵活的渲染管线：

**1. 顶点着色器（Vertex Shader）**
- 输入：顶点属性（位置、法线、UV等）
- 主要任务：
  - 模型视图投影变换（MVP）
  - 顶点动画（骨骼、变形）
  - 计算世界空间坐标和法线
  - 传递属性到片段着色器
- 输出：裁剪空间坐标和属性

**2. 片段着色器（Fragment/Pixel Shader）**
- 输入：插值后的顶点属性
- 主要任务：
  - 光照计算（Phong/Blinn-Phong/PBR）
  - 纹理采样和混合
  - 法线贴图应用
  - 后期处理效果
- 输出：像素颜色（RGBA）

**3. 延迟着色（Deferred Shading）**

*几何通道：*
- 渲染到G-Buffer（几何缓冲）：
  - 法线（世界空间或视图空间）
  - 漫反射颜色
  - 镜面反射参数
  - 深度值
  - 材质ID或其他属性

*光照通道：*
- 对每个光源：
  - 渲染光源影响范围（光照体积）
  - 从G-Buffer读取数据
  - 计算光照贡献
  - 累加到最终图像

*优缺点对比：*

| 特性 | 前向着色 | 延迟着色 |
|------|----------|----------|
| 光照计算复杂度 | O(lights × fragments) | O(lights + fragments) |
| 内存带宽 | 低 | 高（G-Buffer） |
| 透明物体 | 支持 | 不支持 |
| 抗锯齿 | MSAA可用 | MSAA成本高 |
| 材质灵活性 | 高 | 受限于G-Buffer |

**4. Tile-Based Deferred Rendering (TBDR)**
- 将屏幕分成小块（tiles）
- 每个块维护影响它的光源列表
- 结合了前向和延迟着色的优点
- 适合移动GPU架构

## 4.3 纹理映射基础

### 4.3.1 纹理映射的概念

纹理映射（Texture Mapping）是Edwin Catmull于1974年首次提出的技术，通过将2D图像映射到3D表面来增加视觉细节，而无需增加几何复杂度。

**核心思想：**
- 几何提供形状，纹理提供细节
- 将复杂的表面细节“贴”到简单几何体上
- 可以映射的不仅是颜色，还有法线、高度、粗糙度等

**基本流程：**
1. **参数化**：为3D模型的每个顶点指定纹理坐标$(u,v)$
2. **插值**：光栅化时插值纹理坐标（需要透视校正）
3. **采样**：使用插值后的坐标从纹理图像中采样
4. **应用**：将采样结果用于着色计算

**纹理的数学定义：**
纹理是一个函数$T: \mathbb{R}^2 \rightarrow \mathbb{R}^n$，其中：
- 输入：2D纹理坐标$(u,v)$
- 输出：$n$维属性（颜色、法线、材质参数等）

**纹理映射的优势：**
1. **内存效率**：2D图像比3D几何简单得多
2. **编辑方便**：可以使用成熟的图像编辑工具
3. **重用性**：同一纹理可应用于多个模型
4. **LOD支持**：Mipmap等技术提供自动细节层次

### 4.3.2 纹理坐标系统

**UV坐标：**
- $u$：横向坐标，通常范围$[0,1]$（也称为s坐标）
- $v$：纵向坐标，通常范围$[0,1]$（也称为t坐标）
- 原点位置：
  - OpenGL：左下角$(0,0)$
  - DirectX：左上角$(0,0)$

**纹理坐标的获取方式：**

**1. 手动指定**
- 美术师在建模软件中精确设置
- 适合复杂模型和艺术效果
- 需要专业技能和经验

**2. 参数化映射**

*平面映射：*
$$u = \frac{x - x_{min}}{x_{max} - x_{min}}, \quad v = \frac{y - y_{min}}{y_{max} - y_{min}}$$

*球面映射：*
$$u = \frac{1}{2\pi}\arctan\left(\frac{y}{x}\right) + 0.5$$
$$v = \frac{1}{\pi}\arccos\left(\frac{z}{r}\right)$$

*圆柱映射：*
$$u = \frac{1}{2\pi}\arctan\left(\frac{y}{x}\right) + 0.5$$
$$v = \frac{z - z_{min}}{z_{max} - z_{min}}$$

*立方体映射：*
根据法线方向选择六个面之一，然后进行平面映射

**3. 自动UV展开**
- **目标**：最小化扭曲，最大化利用率
- **常用算法**：
  - LSCM (Least Squares Conformal Maps)
  - ABF (Angle Based Flattening)
  - ARAP (As-Rigid-As-Possible)
- **质量指标**：
  - 角度扭曲：$\sum(\theta_{3D} - \theta_{2D})^2$
  - 面积扭曲：$\sum(A_{3D} - A_{2D})^2$
  - 拉伸率：$\max(\sigma_1/\sigma_2)$

**UV岛（UV Islands）：**
- 复杂模型通常被切割成多个UV岛
- 切缝位置的选择很重要：
  - 隐藏在不显眼的地方
  - 沿着自然边界（如衣服缝合线）
  - 避免高曲率区域

### 4.3.3 纹理采样

纹理采样是根据纹理坐标从纹理图像中获取颜色值的过程。不同的采样方法在质量和性能之间做出不同权衡。

**最近邻采样（Nearest Neighbor/Point Sampling）：**
```
u_pixel = round(u * (width - 1))
v_pixel = round(v * (height - 1))
color = texture[u_pixel, v_pixel]
```
- 优点：快速，保持像素化效果，适合像素艺术风格
- 缺点：产生锯齿，放大时出现块状伪影
- 应用：Minecraft等像素风格游戏

**双线性插值（Bilinear Interpolation）：**

计算步骤：
1. 找到四个最近纹素：
   ```
   u_coord = u * (width - 1)
   v_coord = v * (height - 1)
   u0 = floor(u_coord), u1 = u0 + 1
   v0 = floor(v_coord), v1 = v0 + 1
   ```

2. 计算分数部分：
   ```
   s = u_coord - u0
   t = v_coord - v0
   ```

3. 双线性插值：
   $$c = (1-s)(1-t)c_{00} + s(1-t)c_{10} + (1-s)t c_{01} + st c_{11}$$

可以分解为两次线性插值：
- 水平方向：$c_0 = (1-s)c_{00} + s c_{10}$，$c_1 = (1-s)c_{01} + s c_{11}$
- 垂直方向：$c = (1-t)c_0 + t c_1$

**双三次插值（Bicubic Interpolation）：**
使用16个最近纹素（4×4网格），通过三次多项式插值：
$$c(s,t) = \sum_{i=-1}^{2}\sum_{j=-1}^{2} c_{ij} \cdot h(s-i) \cdot h(t-j)$$

其中$h(x)$是Catmull-Rom核函数：
$$h(x) = \begin{cases}
1.5|x|^3 - 2.5|x|^2 + 1 & |x| \leq 1 \\
-0.5|x|^3 + 2.5|x|^2 - 4|x| + 2 & 1 < |x| \leq 2 \\
0 & |x| > 2
\end{cases}$$

- 优点：更平滑的结果，适合高质量放大
- 缺点：计算成本高，可能产生振铃（ringing）

### 4.3.4 纹理放大与缩小

纹理采样的核心挑战是处理纹理分辨率与屏幕分辨率不匹配的情况。

**放大（Magnification）问题：**

*现象：*
- 纹理分辨率低于屏幕采样率
- 一个纹素覆盖多个像素
- 出现块状或模糊伪影

*解决方案：*
1. **双线性插值**：基本选择，平滑但模糊
2. **双三次插值**：更高质量，但计算昂贵
3. **超分辨率技术**：
   - ESRGAN等深度学习方法
   - 实时超分（DLSS、FSR）

**缩小（Minification）问题：**

*现象：*
- 纹理分辨率高于屏幕采样率
- 多个纹素映射到一个像素
- 出现走样（aliasing）和摩尔纹（moiré patterns）

*走样的根源：*
根据奈奎斯特采样定理，采样频率必须至少是信号最高频率的两倍：
$$f_s > 2 f_{max}$$

当纹理被缩小时，高频细节超过了屏幕采样率的奈奎斯特限制。

*解决方案：*
1. **超采样**：增加采样率，但成本高
2. **预滤波**：在采样前去除高频分量
3. **Mipmap**：预计算的多分辨率纹理金字塔

**纹理过滤模式设置：**
```
GL_TEXTURE_MIN_FILTER: // 缩小时
  - GL_NEAREST
  - GL_LINEAR
  - GL_NEAREST_MIPMAP_NEAREST
  - GL_LINEAR_MIPMAP_LINEAR
  
GL_TEXTURE_MAG_FILTER: // 放大时
  - GL_NEAREST
  - GL_LINEAR
```

### 4.3.5 Mipmap技术

Mipmap（来自拉丁语"multum in parvo"，意为"小中见大"）是Lance Williams于1983年提出的预滤波技术。

**Mipmap构建：**
1. **金字塔结构**：
   - Level 0：原始纹理（$n \times n$）
   - Level 1：$n/2 \times n/2$
   - Level 2：$n/4 \times n/4$
   - ...
   - Level $\log_2(n)$：$1 \times 1$

2. **内存开销**：
   $$\text{Total} = n^2 + \frac{n^2}{4} + \frac{n^2}{16} + ... = \frac{4}{3}n^2$$
   只增加33%的存储空间

3. **生成方法**：
   - 简单平均：每2×2像素平均为1个
   - 高质量滤波：使用Lanczos或Mitchell滤波器

**Mipmap级别选择：**

*基本思想：*
选择使得纹素与像素1:1对应的级别。

*计算公式：*
$$L = \log_2\left(\max\left(\sqrt{\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial x}\right)^2} \cdot width, \sqrt{\left(\frac{\partial u}{\partial y}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2} \cdot height\right)\right)$$

*实际实现：*
- 使用差分近似偏导数
- GPU硬件自动计算
- 可以加上LOD bias调整

**三线性插值（Trilinear Interpolation）：**

当$L$不是整数时，在相邻两个mipmap级别之间插值：

1. 计算两个级别：
   ```
   L0 = floor(L)
   L1 = L0 + 1
   frac = L - L0
   ```

2. 在每个级别进行双线性插值：
   ```
   color0 = bilinear_sample(mipmap[L0], u, v)
   color1 = bilinear_sample(mipmap[L1], u, v)
   ```

3. 级别间插值：
   ```
   final_color = (1 - frac) * color0 + frac * color1
   ```

**Mipmap的局限性：**
1. **过度模糊**：各向同性假设导致斜视角模糊
2. **方块伪影**：不同级别之间的过渡可能可见
3. **内存带宽**：需要访问多个级别

**各向异性过滤（Anisotropic Filtering）：**

*问题：*
Mipmap假设纹理在两个方向上均匀缩放，但实际上：
- 地面和墙壁在斜视角下非均匀压缩
- 需要椭圆形而非圆形的采样区域

*解决方案：*
1. **Ripmap**：各向异性mipmap，存储不同长宽比
2. **EWA滤波**：椭圆加权平均
3. **硬件AF**：沿各向异性方向多次采样

*各向异性级别：*
- 2x AF：2个样本
- 4x AF：4个样本
- 16x AF：16个样本（现代标准）

### 4.3.6 纹理应用

纹理不仅可以存储颜色，还可以存储各种表面属性。现代渲染中常用多张纹理来定义复杂材质。

**1. 漫反射贴图（Diffuse/Albedo Map）：**
- 存储基础颜色，不包含光照信息
- 通常为sRGB格式
- 应用：$\mathbf{c}_d = \text{texture}(u, v)$
- PBR中称为Albedo，表示反照率

**2. 法线贴图（Normal Map）：**

*切线空间法线贴图：*
- RGB通道存储切线空间法线
- 红色(R)：右方向(+X)
- 绿色(G)：上方向(+Y)
- 蓝色(B)：外方向(+Z)
- 转换：$\mathbf{n}_{tangent} = 2 \cdot \text{texture}(u,v) - 1$

*世界空间法线贴图：*
- 直接存储世界空间法线
- 不需要TBN矩阵转换
- 但不能重用于不同朝向的表面

**3. 凹凸贴图（Bump Map）：**
- 灰度图，存储高度偏移
- 通过有限差分计算法线扰动：
  $$\mathbf{n}_{new} = \mathbf{n} + k \cdot (\frac{\partial h}{\partial u} \mathbf{T} + \frac{\partial h}{\partial v} \mathbf{B})$$
- 优点：单通道，节省内存
- 缺点：需要实时计算法线

**4. 位移贴图（Displacement Map）：**
- 真正改变几何位置
- 需要细分曲面支持
- 应用：$\mathbf{p}_{new} = \mathbf{p} + h(u,v) \cdot \mathbf{n}$

**5. 环境贴图（Environment Map）：**

*立方体贴图（Cube Map）：*
- 6个面，每个90°视野
- 根据方向向量选择面和UV
- 硬件支持好，无扭曲

*球面贴图（Sphere Map）：*
- 单张图像，经纬度映射
- 极点扭曲严重
- 转换：$(u,v) = (\frac{\theta}{2\pi}, \frac{\phi}{\pi})$

*双抛物面贴图（Dual-Paraboloid）：*
- 两张图像，分别表示半球
- 质量介于立方体和球面之间

**6. PBR纹理组：**
- **Base Color**：基础颜色
- **Metallic**：金属度(0-1)
- **Roughness**：粗糙度(0-1)
- **Normal**：法线贴图
- **AO**：环境遮蔽
- **Height/Displacement**：高度信息

**7. 特殊效果纹理：**
- **Light Map**：预计算光照
- **Shadow Map**：阴影信息
- **Gloss Map**：镜面反射强度
- **Emission Map**：自发光
- **Opacity Map**：透明度

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
