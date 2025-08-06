# 第1章：计算机图形学概述与数学基础

计算机图形学是研究如何通过计算生成、操作和显示视觉内容的学科。本章将从宏观视角介绍图形学的核心概念，并建立必要的数学基础。我们将深入探讨向量空间、线性变换以及四元数等高级数学工具，这些内容构成了后续章节的理论基石。通过本章学习，读者将掌握图形学中的基本数学语言，并理解各种数学工具在实际应用中的优势与局限。

## 1.1 计算机图形学概述

### 1.1.1 图形学的本质与应用领域

计算机图形学的核心任务是将三维世界的几何信息转换为二维图像。这个过程涉及几何建模、光照计算、投影变换和图像合成等多个环节。在现代应用中，图形学已经渗透到游戏开发、电影特效、科学可视化、虚拟现实、计算机辅助设计等众多领域。

从数学角度看，图形学本质上是一系列函数映射：
$$f: \mathbb{R}^3 \times \mathcal{L} \times \mathcal{M} \rightarrow \mathbb{R}^2 \times \mathcal{C}$$

其中$\mathbb{R}^3$表示三维空间坐标，$\mathcal{L}$表示光照条件，$\mathcal{M}$表示材质属性，$\mathbb{R}^2$表示图像平面坐标，$\mathcal{C}$表示颜色空间。

**核心应用领域的技术要求**

不同应用领域对图形学技术有不同侧重：

1. **游戏引擎**：
   - 实时性要求：16.67ms/帧（60 FPS）
   - 大规模场景管理：LOD、流式加载、遮挡剔除
   - 物理模拟集成：碰撞检测、刚体动力学
   - 特效系统：粒子系统、后处理管线

2. **电影视觉特效**：
   - 物理准确性：光谱渲染、次表面散射
   - 超高分辨率：4K-8K渲染，每帧可能需要数小时
   - 复杂材质：毛发、皮肤、流体
   - 合成技术：绿幕抠像、动作捕捉集成

3. **科学可视化**：
   - 大数据处理：TB级数据的实时可视化
   - 准确性要求：色彩映射的感知均匀性
   - 交互式探索：切片、等值面、流线
   - 不确定性可视化：误差传播、置信区间

4. **计算机辅助设计（CAD）**：
   - 精确几何表示：NURBS、细分曲面
   - 实时编辑：参数化建模、约束求解
   - 制造集成：STL导出、切片算法
   - 仿真验证：有限元分析可视化

5. **虚拟/增强现实（VR/AR）**：
   - 超低延迟：<20ms运动到光子延迟
   - 立体渲染：双目视差、会聚调节冲突
   - 注视点渲染：中心凹渲染优化
   - 空间定位：SLAM、深度重建

**图形学的层次化视角**

从不同抽象层次理解图形学：

1. **信号处理层**：图像是二维信号，渲染是重建过程
   $$I(x,y) = \int\int h(x-x', y-y') S(x', y') dx' dy'$$
   其中$h$是重建核，$S$是理想信号

2. **几何变换层**：坐标系变换链
   $$\mathbf{p}_{screen} = \mathbf{M}_{viewport} \cdot \mathbf{M}_{proj} \cdot \mathbf{M}_{view} \cdot \mathbf{M}_{model} \cdot \mathbf{p}_{local}$$

3. **光传输层**：辐射度量学描述光能传播
   $$L(\mathbf{x}, \omega) = L_e(\mathbf{x}, \omega) + \int_\Omega f_r(\mathbf{x}, \omega', \omega) L(\mathbf{x}', -\omega') G(\mathbf{x}, \mathbf{x}') V(\mathbf{x}, \mathbf{x}') d\omega'$$

**渲染的信息论视角**

从信息论角度，渲染是信息压缩过程：
- **输入信息**：3D场景描述，信息熵$H(S)$
- **输出信息**：2D图像，信息熵$H(I)$
- **信息损失**：$I(S;I) = H(S) - H(S|I)$

理想渲染器最大化互信息$I(S;I)$，同时满足人类视觉系统的感知约束。

**渲染方程的统一框架**

更严格地说，渲染的本质是求解渲染方程（Kajiya, 1986）：
$$L_o(\mathbf{x}, \omega_o) = L_e(\mathbf{x}, \omega_o) + \int_{\Omega} f_r(\mathbf{x}, \omega_i, \omega_o) L_i(\mathbf{x}, \omega_i) (\omega_i \cdot \mathbf{n}) d\omega_i$$

其中：
- $L_o$：出射辐射度（radiance）
- $L_e$：自发光辐射度
- $f_r$：双向反射分布函数（BRDF）
- $L_i$：入射辐射度
- $\Omega$：半球积分域

这个积分方程统一了局部光照和全局光照，是现代渲染算法的理论基础。

**图形学的计算复杂度挑战**

渲染的计算复杂度主要来自：
1. **几何复杂度**：$O(n)$个多边形，每个需要投影和裁剪
2. **着色复杂度**：$O(p)$个像素，每个需要光照计算
3. **光照复杂度**：$O(l)$个光源，可能需要阴影计算
4. **全局光照**：递归的光线传播，复杂度可达$O(n^2p)$

因此，总复杂度在最坏情况下为$O(nlp)$，这解释了为什么需要各种加速算法。

**深入理解复杂度来源**

让我们分析一个典型场景的计算量：
- 场景几何：$n = 10^6$ 个三角形
- 屏幕分辨率：$p = 1920 \times 1080 \approx 2 \times 10^6$ 像素
- 光源数量：$l = 100$ 个光源
- 光线弹射次数：$d = 5$ 次（全局光照）

基础光栅化：$O(n \log n + p)$ （深度排序 + Z-buffer）
带阴影的光栅化：$O(nl + p \cdot l)$ （每个光源一个shadow map）
路径追踪：$O(p \cdot d \cdot \log n)$ （每像素多条光线）

实际计算量：
- 光栅化：约$10^7$次操作
- 带阴影：约$10^8$次操作  
- 路径追踪：约$10^8$次光线求交

这解释了为什么需要：
- 硬件加速（GPU并行化）
- 算法优化（空间数据结构）
- 近似技术（屏幕空间效果）

**复杂度优化策略**

现代图形学通过多种策略降低复杂度：

1. **空间数据结构**：
   - BVH（包围体层次）：构建$O(n\log n)$，查询$O(\log n)$
   - KD-Tree：空间划分，平均查询$O(\log n)$
   - 八叉树：稀疏体素化，内存效率高

2. **时间相干性**：
   - 帧间复用：$I_t = \alpha I_{t-1} + (1-\alpha)I_{new}$
   - 运动向量：$\mathbf{p}_t = \mathbf{p}_{t-1} + \mathbf{v}\Delta t$
   - 时序抗锯齿（TAA）：累积多帧样本

3. **层次化加速**：
   - LOD（细节层次）：$L(d) = L_{max} \cdot \max(0, 1 - \frac{d}{d_{max}})$
   - Mipmap：预过滤纹理，$O(1)$查询
   - 视锥体裁剪：$\mathbf{n} \cdot (\mathbf{p} - \mathbf{p}_0) \geq 0$

4. **近似算法**：
   - 屏幕空间技术：SSAO、SSR
   - 预计算辐照度：球谐函数展开
   - 神经网络降噪：$\tilde{I} = D_\theta(I_{noisy})$

### 1.1.2 实时渲染与离线渲染

图形学应用可分为两大类：

**实时渲染**：追求在有限时间内（通常16-33毫秒）生成图像，常用于游戏和交互应用。主要技术包括：
- 光栅化管线（Z-buffer、着色器）
- 层次细节（LOD）优化
- 时间相关性利用（时序抗锯齿、运动模糊）

**离线渲染**：追求物理准确性和视觉真实感，计算时间可达数小时。核心技术包括：
- 光线追踪与路径追踪
- 全局光照算法
- 基于物理的渲染（PBR）

**实时渲染的关键技术栈**

1. **可见性判断**：
   - 视锥体裁剪：平面方程 $\mathbf{n} \cdot \mathbf{p} + d = 0$
   - 遮挡剔除：硬件遮挡查询、Hi-Z buffer
   - 背面剔除：$(\mathbf{v}_2 - \mathbf{v}_0) \times (\mathbf{v}_1 - \mathbf{v}_0) \cdot \mathbf{view} < 0$
   - 早期深度测试：避免无用的片段着色

2. **着色优化**：
   - 延迟渲染：几何pass + 光照pass分离
   - 前向+渲染：tiled光照剔除
   - 可变速率着色：根据内容动态调整采样率
   - 着色器LOD：距离相关的着色复杂度

3. **抗锯齿技术**：
   - MSAA：硬件多重采样，$4\times$或$8\times$
   - FXAA/SMAA：后处理边缘检测
   - TAA：时间累积 + 运动向量
   - DLSS：AI驱动的超采样

**离线渲染的质量追求**

1. **光传输模拟**：
   - 双向路径追踪（BDPT）：从相机和光源同时追踪
   - Metropolis光传输（MLT）：马尔可夫链采样
   - 光子映射：因果光子 + 最终聚集
   - 辐照度缓存：稀疏采样 + 插值

2. **材质复杂度**：
   - 分层材质：清漆、基底多层混合
   - 体积散射：次表面散射、云雾
   - 光谱渲染：避免同色异谱问题
   - 偏振光：晶体、水面真实感

3. **几何细节**：
   - 置换贴图：真实几何细节
   - 自适应细分：屏幕空间误差控制
   - 毛发渲染：Marschner模型
   - 流体模拟：SPH、FLIP求解器

**性能指标的数学分析**

实时渲染的帧率要求可以表述为：
$$t_{frame} = t_{CPU} + t_{GPU} \leq \frac{1}{FPS_{target}}$$

其中：
- $t_{CPU}$：CPU计算时间（场景管理、物理模拟）
- $t_{GPU}$：GPU渲染时间（顶点处理、片段着色）
- $FPS_{target}$：目标帧率（如60 FPS → 16.67ms）

**算法复杂度对比**

| 渲染方法 | 几何复杂度 | 像素复杂度 | 内存需求 | 典型应用 |
|---------|-----------|-----------|----------|---------|
| 光栅化 | $O(n\log n)$ | $O(p)$ | $O(p)$ (Z-buffer) | 游戏引擎 |
| 光线追踪 | $O(p\log n)$ | $O(pd)$ | $O(n)$ (BVH) | 电影特效 |
| 辐射度 | $O(n^2)$ | $O(p)$ | $O(n^2)$ | 建筑可视化 |

这里$n$是场景中的图元数，$p$是像素数，$d$是递归深度。

**实时与离线的界限模糊化**

现代GPU的发展使得这两者的界限越来越模糊：
- **混合渲染**：实时光线追踪用于反射和阴影，光栅化用于主要可见性
- **时间积累**：利用多帧累积提高质量（如DLSS）
- **预计算技术**：光照贴图、环境光遮蔽、辐照度探针

### 1.1.3 图形管线概述

现代图形管线可以抽象为以下阶段：

1. **应用阶段**：场景管理、碰撞检测、动画更新
2. **几何处理**：顶点变换、裁剪、图元装配
3. **光栅化**：将几何图元转换为像素片段
4. **片段处理**：着色计算、深度测试、混合

每个阶段都涉及特定的数学变换和算法优化，我们将在后续章节详细探讨。

**管线的演化历史**

1. **固定功能管线时代**（1990s）：
   - 硬编码的变换和光照模型
   - Gouraud着色、Phong着色
   - 有限的纹理单元和混合模式
   - OpenGL 1.x、DirectX 7

2. **可编程着色器时代**（2000s）：
   - 顶点着色器：自定义顶点变换
   - 片段着色器：自定义像素着色
   - 着色器模型演进：SM 2.0 → 5.0
   - GLSL、HLSL、Cg语言

3. **统一着色器架构**（2010s）：
   - 几何着色器：图元级操作
   - 曲面细分着色器：动态细分
   - 计算着色器：通用GPU计算
   - Vulkan、DirectX 12低级API

4. **现代GPU创新**（2020s）：
   - 光线追踪核心：硬件BVH遍历
   - 张量核心：AI加速
   - 网格着色器：GPU驱动的渲染
   - 可变速率着色：自适应质量

**管线阶段的详细分析**

各阶段的核心算法和数据流：

1. **应用阶段（CPU）**：
   - 场景图遍历：$O(n)$复杂度
   - 视锥体裁剪：6个平面测试
   - 物理模拟：$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \mathbf{v}_t\Delta t + \frac{1}{2}\mathbf{a}_t\Delta t^2$
   - 动画更新：骨骼矩阵插值

2. **顶点处理（GPU）**：
   - 模型变换：$\mathbf{v}_{world} = \mathbf{M}_{model} \cdot \mathbf{v}_{local}$
   - 视图变换：$\mathbf{v}_{view} = \mathbf{M}_{view} \cdot \mathbf{v}_{world}$
   - 投影变换：$\mathbf{v}_{clip} = \mathbf{M}_{proj} \cdot \mathbf{v}_{view}$
   - 透视除法：$\mathbf{v}_{ndc} = \mathbf{v}_{clip} / w_{clip}$

3. **图元装配与裁剪**：
   - Sutherland-Hodgman算法
   - 齐次裁剪空间：$-w \leq x,y,z \leq w$
   - 背面剔除：$(\mathbf{v}_1 - \mathbf{v}_0) \times (\mathbf{v}_2 - \mathbf{v}_0) \cdot \mathbf{view} < 0$

4. **光栅化**：
   - 扫描线算法：增量式边缘计算
   - 重心坐标：$\lambda_i = \frac{A_i}{A_{total}}$
   - 深度插值：$z = \lambda_0 z_0 + \lambda_1 z_1 + \lambda_2 z_2$

5. **片段处理**：
   - Early-Z测试：$z_{frag} < z_{buffer}$
   - 着色计算：BRDF评估
   - Alpha混合：$\mathbf{c}_{out} = \alpha_{src}\mathbf{c}_{src} + (1-\alpha_{src})\mathbf{c}_{dst}$

**管线的数学表示**

图形管线可以表示为一系列函数的复合：
$$\mathcal{P} = \mathcal{F} \circ \mathcal{R} \circ \mathcal{G} \circ \mathcal{A}$$

其中：
- $\mathcal{A}$：应用阶段，$\mathcal{A}: \mathcal{S} \rightarrow \mathcal{V}$（场景→顶点）
- $\mathcal{G}$：几何阶段，$\mathcal{G}: \mathcal{V} \rightarrow \mathcal{T}$（顶点→变换后的图元）
- $\mathcal{R}$：光栅化，$\mathcal{R}: \mathcal{T} \rightarrow \mathcal{F}$（图元→片段）
- $\mathcal{F}$：片段阶段，$\mathcal{F}: \mathcal{F} \rightarrow \mathcal{I}$（片段→图像）

**可编程管线的革命**

现代GPU提供了可编程阶段：
- **顶点着色器**：$\mathbf{v}_{clip} = \mathbf{M}_{proj} \cdot \mathbf{M}_{view} \cdot \mathbf{M}_{model} \cdot \mathbf{v}_{local}$
- **几何着色器**：可以创建或销毁图元
- **片段着色器**：计算最终颜色 $\mathbf{c} = f_{shader}(\mathbf{n}, \mathbf{l}, \mathbf{v}, \mathbf{material})$
- **计算着色器**：通用并行计算

**并行化的数学基础**

GPU的高效来自于大规模并行：
- **数据并行**：$\forall i \in [0, N): y_i = f(x_i)$
- **SIMD执行**：相同指令作用于多个数据
- **内存访问模式**：合并访问、纹理缓存

并行效率可以用Amdahl定律估算：
$$S = \frac{1}{(1-p) + \frac{p}{n}}$$

其中$p$是可并行化部分的比例，$n$是处理器数量。

### 1.1.4 现代图形学的挑战与机遇

随着AI技术的发展，图形学面临新的机遇：

- **神经渲染**：使用深度学习加速或替代传统渲染算法
- **可微渲染**：将渲染过程设计为可微函数，支持基于梯度的优化
- **实时光线追踪**：硬件加速使得实时光线追踪成为可能
- **虚拟现实与增强现实**：对延迟和分辨率提出更高要求

**新兴技术的深度分析**

1. **神经场景表示**：
   - NeRF类方法：连续体素表示
   - 3D Gaussian Splatting：离散高斯表示
   - Hash编码：$\gamma(p) = [\sin(2^0\pi p), \cos(2^0\pi p), ..., \sin(2^L\pi p), \cos(2^L\pi p)]$
   - 优势：压缩率高，新视角合成质量好

2. **混合渲染管线**：
   - 光栅化主路径 + 光线追踪辅助
   - 降噪器：$\mathcal{D}: I_{noisy} \rightarrow I_{clean}$
   - 重要性采样：$p(\omega) \propto L_i(\omega)f_r(\omega)|\cos\theta|$
   - 自适应采样：基于方差估计动态分配样本

3. **实时全局光照近似**：
   - Voxel Cone Tracing：体素化场景锥形采样
   - Radiance Probes：稀疏辐照度采样
   - Screen Space GI：屏幕空间光线步进
   - 混合方案：远场用探针，近场用屏幕空间

4. **AI驱动的内容生成**：
   - 程序化纹理：$T = G_\theta(z, \mathbf{uv})$
   - 风格迁移：$I_{styled} = \arg\min_I \mathcal{L}_{content}(I, I_c) + \alpha\mathcal{L}_{style}(I, I_s)$
   - 超分辨率：端到端学习上采样
   - 时序稳定性：循环一致性约束

**神经渲染的数学框架**

神经渲染将传统渲染函数$f$替换为神经网络$f_\theta$：
$$\mathbf{I} = f_\theta(\mathbf{G}, \mathbf{V}, \mathbf{L})$$

其中$\theta$是网络参数，通过最小化损失函数学习：
$$\mathcal{L} = \mathbb{E}_{(\mathbf{I}_{gt}, \mathbf{x})} \left[ \|\mathbf{I}_{gt} - f_\theta(\mathbf{x})\|^2 + \lambda R(\theta) \right]$$

典型应用包括：
- **神经辐射场（NeRF）**：$F_\theta: (\mathbf{x}, \mathbf{d}) \rightarrow (\mathbf{c}, \sigma)$
- **神经纹理**：学习的特征图代替传统纹理
- **超分辨率**：$\mathbf{I}_{high} = G_\theta(\mathbf{I}_{low})$

**可微渲染的梯度计算**

可微渲染需要计算$\frac{\partial \mathbf{I}}{\partial \mathbf{p}}$，其中$\mathbf{p}$可以是几何、材质或光照参数。主要挑战包括：

1. **可见性的不连续性**：使用边缘采样或软光栅化
2. **蒙特卡洛积分的梯度**：
   $$\nabla_\theta \mathbb{E}_{p_\theta}[f(x)] = \mathbb{E}_{p_\theta}[f(x) \nabla_\theta \log p_\theta(x)]$$
3. **反向模式自动微分**：高效计算复杂管线的梯度

**硬件加速的新时代**

现代GPU引入了专用硬件：
- **RT Core**：加速光线-三角形相交测试
- **Tensor Core**：加速矩阵运算和AI推理
- **可变速率着色（VRS）**：根据内容自适应调整着色率

性能模型：
$$T_{total} = T_{traverse} + T_{intersect} + T_{shade}$$

其中硬件加速主要优化$T_{intersect}$，可达10-100倍加速。

## 1.2 向量与线性代数

### 1.2.1 向量空间与基本运算

向量是图形学中最基本的数学对象。在$n$维欧几里得空间$\mathbb{R}^n$中，向量$\mathbf{v} = (v_1, v_2, ..., v_n)^T$表示一个有向线段。

**向量加法与标量乘法**：
$$\mathbf{u} + \mathbf{v} = (u_1 + v_1, u_2 + v_2, ..., u_n + v_n)^T$$
$$k\mathbf{v} = (kv_1, kv_2, ..., kv_n)^T$$

这些运算满足向量空间的八条公理，使$\mathbb{R}^n$成为一个线性空间。

**向量的几何解释**

向量可以有多种解释，在图形学中常见的包括：

1. **位移向量**：表示从一点到另一点的移动
   - 平移变换：$\mathbf{p}' = \mathbf{p} + \mathbf{t}$
   - 相对位置：$\mathbf{v}_{AB} = \mathbf{p}_B - \mathbf{p}_A$

2. **方向向量**：表示朝向，通常归一化
   - 表面法线：$\mathbf{n} = \frac{(\mathbf{v}_1 - \mathbf{v}_0) \times (\mathbf{v}_2 - \mathbf{v}_0)}{|(\mathbf{v}_1 - \mathbf{v}_0) \times (\mathbf{v}_2 - \mathbf{v}_0)|}$
   - 光线方向：$\mathbf{d} = \frac{\mathbf{p}_{target} - \mathbf{p}_{origin}}{|\mathbf{p}_{target} - \mathbf{p}_{origin}|}$

3. **速度向量**：表示运动状态
   - 线速度：$\mathbf{v} = \frac{d\mathbf{p}}{dt}$
   - 角速度：$\boldsymbol{\omega} \times \mathbf{r}$（叉积给出切向速度）

4. **力向量**：物理模拟中的作用力
   - 重力：$\mathbf{F}_g = m\mathbf{g}$
   - 弹簧力：$\mathbf{F}_s = -k(\mathbf{x} - \mathbf{x}_0)$

**向量空间的公理系统**

向量空间$(V, +, \cdot)$必须满足：
1. **加法封闭性**：$\forall \mathbf{u}, \mathbf{v} \in V: \mathbf{u} + \mathbf{v} \in V$
2. **加法交换律**：$\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$
3. **加法结合律**：$(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$
4. **零向量存在**：$\exists \mathbf{0} \in V: \mathbf{v} + \mathbf{0} = \mathbf{v}$
5. **加法逆元**：$\forall \mathbf{v} \in V, \exists -\mathbf{v}: \mathbf{v} + (-\mathbf{v}) = \mathbf{0}$
6. **标量乘法封闭**：$\forall k \in \mathbb{R}, \mathbf{v} \in V: k\mathbf{v} \in V$
7. **标量分配律**：$k(\mathbf{u} + \mathbf{v}) = k\mathbf{u} + k\mathbf{v}$
8. **向量分配律**：$(k + l)\mathbf{v} = k\mathbf{v} + l\mathbf{v}$

**仿射空间与向量空间的区别**

图形学中常混淆点和向量：
- **点**：位置，属于仿射空间$\mathcal{A}$
- **向量**：位移，属于向量空间$V$
- **关系**：$\mathbf{p} - \mathbf{q} = \mathbf{v}$（两点之差是向量）
- **仿射组合**：$\mathbf{p} = \sum_i \lambda_i \mathbf{p}_i$，其中$\sum_i \lambda_i = 1$

**内积（点积）**：
$$\mathbf{u} \cdot \mathbf{v} = \sum_{i=1}^{n} u_i v_i = |\mathbf{u}||\mathbf{v}|\cos\theta$$

内积的几何意义：
- 计算投影长度：$\text{proj}_{\mathbf{v}}\mathbf{u} = \frac{\mathbf{u} \cdot \mathbf{v}}{|\mathbf{v}|^2}\mathbf{v}$
- 判断方向关系：$\mathbf{u} \cdot \mathbf{v} > 0$表示夹角小于90°
- 计算功和能量：在物理模拟中广泛应用

**内积的推广形式**

更一般地，内积可以定义为：
$$\langle \mathbf{u}, \mathbf{v} \rangle_{\mathbf{M}} = \mathbf{u}^T \mathbf{M} \mathbf{v}$$

其中$\mathbf{M}$是正定矩阵。这在各向异性材质和椭球体距离计算中很有用。

**外积（叉积）**：仅在三维空间定义
$$\mathbf{u} \times \mathbf{v} = \begin{pmatrix}
u_2v_3 - u_3v_2 \\
u_3v_1 - u_1v_3 \\
u_1v_2 - u_2v_1
\end{pmatrix}$$

外积的性质：
- $|\mathbf{u} \times \mathbf{v}| = |\mathbf{u}||\mathbf{v}|\sin\theta$
- 方向遵循右手定则
- 反交换律：$\mathbf{u} \times \mathbf{v} = -\mathbf{v} \times \mathbf{u}$

**叉积的矩阵形式**

叉积可以表示为反对称矩阵乘法：
$$\mathbf{u} \times \mathbf{v} = [\mathbf{u}]_\times \mathbf{v}$$

其中：
$$[\mathbf{u}]_\times = \begin{pmatrix}
0 & -u_3 & u_2 \\
u_3 & 0 & -u_1 \\
-u_2 & u_1 & 0
\end{pmatrix}$$

这个表示在刚体动力学和相机标定中非常有用。

**三重积公式**

标量三重积：
$$(\mathbf{a} \times \mathbf{b}) \cdot \mathbf{c} = \det[\mathbf{a} | \mathbf{b} | \mathbf{c}]$$

几何意义：三个向量张成的平行六面体的有向体积。

向量三重积（BAC-CAB规则）：
$$\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = \mathbf{b}(\mathbf{a} \cdot \mathbf{c}) - \mathbf{c}(\mathbf{a} \cdot \mathbf{b})$$

### 1.2.2 矩阵与线性变换

线性变换$T: \mathbb{R}^n \rightarrow \mathbb{R}^m$可以用$m \times n$矩阵$\mathbf{A}$表示：
$$T(\mathbf{v}) = \mathbf{A}\mathbf{v}$$

线性变换必须满足：
- 可加性：$T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v})$
- 齐次性：$T(k\mathbf{v}) = kT(\mathbf{v})$

**基本变换矩阵**：

缩放矩阵：
$$\mathbf{S} = \begin{pmatrix}
s_x & 0 & 0 \\
0 & s_y & 0 \\
0 & 0 & s_z
\end{pmatrix}$$

绕$z$轴旋转$\theta$角：
$$\mathbf{R}_z(\theta) = \begin{pmatrix}
\cos\theta & -\sin\theta & 0 \\
\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{pmatrix}$$

**任意轴旋转（Rodrigues公式）**

绕单位轴$\mathbf{k}$旋转$\theta$角的旋转矩阵：
$$\mathbf{R}(\mathbf{k}, \theta) = \mathbf{I} + \sin\theta[\mathbf{k}]_\times + (1-\cos\theta)[\mathbf{k}]_\times^2$$

展开形式：
$$\mathbf{R} = \cos\theta\mathbf{I} + (1-\cos\theta)\mathbf{k}\mathbf{k}^T + \sin\theta[\mathbf{k}]_\times$$

这个公式在骨骼动画和物理模拟中广泛使用。

**矩阵运算的几何意义**：
- 矩阵乘法对应变换的复合
- 逆矩阵对应逆变换
- 行列式表示体积缩放因子
- 特征值和特征向量揭示变换的不变方向

**奇异值分解（SVD）的几何解释**

任意矩阵$\mathbf{A}$可分解为：
$$\mathbf{A} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T$$

几何意义：
1. $\mathbf{V}^T$：第一次旋转（输入空间）
2. $\mathbf{\Sigma}$：沿主轴缩放
3. $\mathbf{U}$：第二次旋转（输出空间）

在图形学中的应用：
- 提取旋转和缩放分量
- 矩阵近似和压缩
- 最小二乘问题求解
- 主成分分析（PCA）

### 1.2.3 齐次坐标与仿射变换

为了统一表示线性变换和平移，引入齐次坐标。三维点$(x, y, z)$的齐次坐标为$(x, y, z, 1)$，向量的齐次坐标为$(x, y, z, 0)$。

仿射变换的一般形式：
$$\mathbf{T} = \begin{pmatrix}
\mathbf{A}_{3×3} & \mathbf{t}_{3×1} \\
\mathbf{0}_{1×3} & 1
\end{pmatrix}$$

其中$\mathbf{A}$是线性变换部分，$\mathbf{t}$是平移向量。

**变换的组合与分解**：
- 任意刚体变换可分解为旋转和平移
- 任意仿射变换可通过SVD分解为旋转、缩放和剪切
- 变换顺序很重要：$\mathbf{T}_1\mathbf{T}_2 \neq \mathbf{T}_2\mathbf{T}_1$

### 1.2.4 正交性与投影

**正交矩阵**满足$\mathbf{Q}^T\mathbf{Q} = \mathbf{I}$，具有以下性质：
- 保持长度和角度
- 行列式为±1
- 逆矩阵等于转置矩阵

**投影矩阵**：
正交投影到子空间$\text{span}(\mathbf{v}_1, ..., \mathbf{v}_k)$：
$$\mathbf{P} = \mathbf{V}(\mathbf{V}^T\mathbf{V})^{-1}\mathbf{V}^T$$

其中$\mathbf{V} = [\mathbf{v}_1 | ... | \mathbf{v}_k]$。

透视投影矩阵（后续章节详述）将三维点映射到规范化设备坐标。

## 1.3 四元数与复数在图形学中的应用

### 1.3.1 复数与二维旋转

复数提供了优雅的二维旋转表示。复数$z = a + bi$可以写成极坐标形式：
$$z = r(\cos\theta + i\sin\theta) = re^{i\theta}$$

二维点$(x, y)$对应复数$x + yi$，绕原点旋转$\phi$角等价于乘以$e^{i\phi}$：
$$(x + yi)e^{i\phi} = (x + yi)(\cos\phi + i\sin\phi)$$
$$= (x\cos\phi - y\sin\phi) + i(x\sin\phi + y\cos\phi)$$

这恰好对应旋转矩阵的作用结果。复数乘法的几何意义：
- 模相乘：$|z_1z_2| = |z_1||z_2|$
- 幅角相加：$\arg(z_1z_2) = \arg(z_1) + \arg(z_2)$

### 1.3.2 四元数基础

四元数是复数的推广，形式为：
$$q = w + xi + yj + zk$$

其中$i^2 = j^2 = k^2 = ijk = -1$，且：
- $ij = k, jk = i, ki = j$
- $ji = -k, kj = -i, ik = -j$

四元数可以写成标量-向量形式：
$$q = (s, \mathbf{v}) = s + v_xi + v_yj + v_zk$$

**四元数运算**：

加法：$(s_1, \mathbf{v}_1) + (s_2, \mathbf{v}_2) = (s_1 + s_2, \mathbf{v}_1 + \mathbf{v}_2)$

乘法（Hamilton积）：
$$q_1q_2 = (s_1s_2 - \mathbf{v}_1 \cdot \mathbf{v}_2, s_1\mathbf{v}_2 + s_2\mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2)$$

共轭：$q^* = (s, -\mathbf{v})$

模：$|q| = \sqrt{qq^*} = \sqrt{s^2 + |\mathbf{v}|^2}$

### 1.3.3 四元数与三维旋转

单位四元数可以表示三维旋转。绕单位轴$\mathbf{n}$旋转$\theta$角对应的四元数为：
$$q = \cos\frac{\theta}{2} + \sin\frac{\theta}{2}(n_xi + n_yj + n_zk)$$

**旋转公式**：点$\mathbf{p}$经四元数$q$旋转后得到：
$$\mathbf{p}' = q\mathbf{p}q^*$$

其中$\mathbf{p}$被视为纯四元数$(0, \mathbf{p})$。

**四元数转旋转矩阵**：
设$q = w + xi + yj + zk$，对应的旋转矩阵为：
$$\mathbf{R} = \begin{pmatrix}
1-2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\
2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\
2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2)
\end{pmatrix}$$

### 1.3.4 四元数的优势与应用

**相比欧拉角的优势**：
1. 无万向锁问题
2. 插值更平滑（球面线性插值SLERP）
3. 数值稳定性更好
4. 存储更紧凑（4个数vs9个数）

**球面线性插值（SLERP）**：
$$\text{slerp}(q_1, q_2, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}q_1 + \frac{\sin(t\Omega)}{\sin\Omega}q_2$$

其中$\cos\Omega = q_1 \cdot q_2$。

**在图形学中的应用**：
- 骨骼动画的关节旋转
- 相机控制（避免万向锁）
- 物理模拟中的刚体方向
- 法线贴图的切线空间变换

### 1.3.5 对偶四元数与刚体变换

对偶四元数将旋转和平移统一表示：
$$\hat{q} = q_r + \epsilon q_d$$

其中$\epsilon^2 = 0$（对偶单位），$q_r$表示旋转，$q_d = \frac{1}{2}tq_r$（$t$是平移四元数）。

对偶四元数的优势：
- 统一表示刚体变换
- 支持平滑插值（包括旋转和平移）
- 避免矩阵乘法的数值误差累积

## 本章小结

本章建立了计算机图形学的数学基础：

**核心概念**：
1. 图形学管线：应用→几何→光栅化→片段处理
2. 向量运算：点积（投影、角度）、叉积（法向、面积）
3. 矩阵变换：线性变换、仿射变换、齐次坐标
4. 四元数：紧凑的旋转表示，避免万向锁

**关键公式**：
- 点积：$\mathbf{u} \cdot \mathbf{v} = |\mathbf{u}||\mathbf{v}|\cos\theta$
- 叉积：$|\mathbf{u} \times \mathbf{v}| = |\mathbf{u}||\mathbf{v}|\sin\theta$
- 四元数旋转：$\mathbf{p}' = q\mathbf{p}q^*$
- SLERP：$\text{slerp}(q_1, q_2, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}q_1 + \frac{\sin(t\Omega)}{\sin\Omega}q_2$

**重要性质**：
- 正交矩阵保持长度和角度
- 单位四元数形成3-球面$S^3$
- 旋转的复合对应四元数乘法
- 仿射变换保持直线和平行关系

## 练习题

### 基础题

**1.1 向量运算基础**
给定向量$\mathbf{a} = (1, 2, 3)$，$\mathbf{b} = (4, -1, 2)$，计算：
- (a) $\mathbf{a} \cdot \mathbf{b}$
- (b) $\mathbf{a} \times \mathbf{b}$
- (c) $\mathbf{a}$在$\mathbf{b}$上的投影向量

*Hint*: 投影公式为$\text{proj}_{\mathbf{b}}\mathbf{a} = \frac{\mathbf{a} \cdot \mathbf{b}}{|\mathbf{b}|^2}\mathbf{b}$

<details>
<summary>答案</summary>

(a) $\mathbf{a} \cdot \mathbf{b} = 1×4 + 2×(-1) + 3×2 = 8$

(b) $\mathbf{a} \times \mathbf{b} = \begin{pmatrix} 2×2-3×(-1) \\ 3×4-1×2 \\ 1×(-1)-2×4 \end{pmatrix} = \begin{pmatrix} 7 \\ 10 \\ -9 \end{pmatrix}$

(c) $|\mathbf{b}|^2 = 16 + 1 + 4 = 21$，所以投影向量为$\frac{8}{21}(4, -1, 2) = (\frac{32}{21}, -\frac{8}{21}, \frac{16}{21})$

</details>

**1.2 矩阵变换组合**
一个物体先绕$z$轴旋转45°，再沿$x$轴方向缩放2倍，最后平移$(3, 1, 0)$。写出复合变换矩阵。

*Hint*: 注意变换顺序，使用齐次坐标

<details>
<summary>答案</summary>

变换顺序（从右到左）：旋转→缩放→平移

$$\mathbf{T} = \begin{pmatrix}
1 & 0 & 0 & 3 \\
0 & 1 & 0 & 1 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
2 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
\cos45° & -\sin45° & 0 & 0 \\
\sin45° & \cos45° & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}$$

$$= \begin{pmatrix}
\sqrt{2} & -\frac{\sqrt{2}}{2} & 0 & 3 \\
\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2} & 0 & 1 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}$$

</details>

**1.3 四元数旋转**
将四元数$q = \frac{1}{2} + \frac{1}{2}i + \frac{1}{2}j + \frac{1}{2}k$转换为旋转矩阵。

*Hint*: 首先验证这是单位四元数

<details>
<summary>答案</summary>

验证：$|q| = \sqrt{(\frac{1}{2})^2 + (\frac{1}{2})^2 + (\frac{1}{2})^2 + (\frac{1}{2})^2} = 1$ ✓

设$w = \frac{1}{2}, x = y = z = \frac{1}{2}$，代入公式：

$$\mathbf{R} = \begin{pmatrix}
1-2(\frac{1}{4}+\frac{1}{4}) & 2(\frac{1}{4}-\frac{1}{4}) & 2(\frac{1}{4}+\frac{1}{4}) \\
2(\frac{1}{4}+\frac{1}{4}) & 1-2(\frac{1}{4}+\frac{1}{4}) & 2(\frac{1}{4}-\frac{1}{4}) \\
2(\frac{1}{4}-\frac{1}{4}) & 2(\frac{1}{4}+\frac{1}{4}) & 1-2(\frac{1}{4}+\frac{1}{4})
\end{pmatrix}$$

$$= \begin{pmatrix}
0 & 0 & 1 \\
1 & 0 & 0 \\
0 & 1 & 0
\end{pmatrix}$$

这是一个120°旋转，轴为$(1,1,1)$方向。

</details>

### 挑战题

**1.4 最小二乘投影**
推导将向量$\mathbf{b}$投影到由列向量$\mathbf{a}_1, \mathbf{a}_2$张成的平面上的投影矩阵$\mathbf{P}$，并证明$\mathbf{P}^2 = \mathbf{P}$。

*Hint*: 使用正规方程$\mathbf{A}^T\mathbf{A}\mathbf{x} = \mathbf{A}^T\mathbf{b}$

<details>
<summary>答案</summary>

设$\mathbf{A} = [\mathbf{a}_1 | \mathbf{a}_2]$，投影为$\mathbf{p} = \mathbf{A}\mathbf{x}$。

最小化$|\mathbf{b} - \mathbf{A}\mathbf{x}|^2$，得正规方程：
$$\mathbf{A}^T\mathbf{A}\mathbf{x} = \mathbf{A}^T\mathbf{b}$$

解得：$\mathbf{x} = (\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\mathbf{b}$

投影：$\mathbf{p} = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\mathbf{b}$

因此投影矩阵：$\mathbf{P} = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T$

证明幂等性：
$$\mathbf{P}^2 = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T\mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T = \mathbf{A}(\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T = \mathbf{P}$$

</details>

**1.5 四元数插值优化**
给定两个单位四元数$q_1$和$q_2$，当$q_1 \cdot q_2 < 0$时，SLERP会选择"长路径"。如何修正这个问题？证明你的方法保持插值的连续性。

*Hint*: 考虑四元数的双覆盖性质

<details>
<summary>答案</summary>

四元数$q$和$-q$表示相同的旋转（双覆盖）。当$q_1 \cdot q_2 < 0$时，$q_1$和$q_2$在4D球面上的夹角大于90°。

修正方法：若$q_1 \cdot q_2 < 0$，将$q_2$替换为$-q_2$。

证明连续性：
设原始插值路径为$\gamma(t) = \text{slerp}(q_1, q_2, t)$，修正后为$\gamma'(t) = \text{slerp}(q_1, -q_2, t)$。

在$t=0$：$\gamma'(0) = q_1 = \gamma(0)$
在$t=1$：$\gamma'(1) = -q_2$，但$-q_2$和$q_2$表示相同旋转

修正后的夹角$\Omega' = \arccos(-q_1 \cdot q_2) = \pi - \Omega < \pi/2$，确保选择短路径。

</details>

**1.6 仿射变换的不动点**
证明：任意二维仿射变换$T(\mathbf{x}) = \mathbf{A}\mathbf{x} + \mathbf{b}$（其中$\det(\mathbf{A}-\mathbf{I}) \neq 0$）有唯一不动点。找出该不动点的表达式。

*Hint*: 不动点满足$T(\mathbf{x}) = \mathbf{x}$

<details>
<summary>答案</summary>

设不动点为$\mathbf{x}_0$，则：
$$\mathbf{A}\mathbf{x}_0 + \mathbf{b} = \mathbf{x}_0$$
$$(\mathbf{A} - \mathbf{I})\mathbf{x}_0 = -\mathbf{b}$$

由于$\det(\mathbf{A}-\mathbf{I}) \neq 0$，矩阵$\mathbf{A}-\mathbf{I}$可逆，因此：
$$\mathbf{x}_0 = -(\mathbf{A} - \mathbf{I})^{-1}\mathbf{b}$$

这是唯一解。

几何意义：不动点是变换的"中心"，所有点围绕它进行变换。

</details>

**1.7 复数与四元数的联系**
将二维旋转的复数表示推广到三维。具体地，如何用两个复数表示一个三维旋转？这与四元数有什么关系？

*Hint*: 考虑Cayley-Dickson构造

<details>
<summary>答案</summary>

四元数可以看作"复数对"：$q = z_1 + z_2j$，其中$z_1 = a + bi$，$z_2 = c + di$。

这给出：$q = a + bi + cj + dk$

关系：
- 当$c = d = 0$时，退化为复数（二维旋转）
- 四元数乘法规则可从复数乘法推广得出

另一种表示：使用两个复数$(z_1, z_2)$表示三维旋转，满足$|z_1|^2 + |z_2|^2 = 1$。这实际上是$SU(2)$群的表示，与单位四元数同构。

对应关系：
$$\begin{pmatrix} z_1 & -\bar{z}_2 \\ z_2 & \bar{z}_1 \end{pmatrix} \leftrightarrow a + bi + cj + dk$$

其中$z_1 = a + bi$，$z_2 = c + di$。

</details>

**1.8 图形管线的能量守恒**
在物理真实的渲染中，为什么需要确保BRDF满足能量守恒？给出数学表述并解释其在实时渲染中的简化方法。

*Hint*: 考虑反射率的积分

<details>
<summary>答案</summary>

能量守恒要求出射能量不超过入射能量。对于BRDF $f_r(\omega_i, \omega_o)$：

$$\int_{\Omega} f_r(\omega_i, \omega_o) \cos\theta_i \, d\omega_i \leq 1$$

对所有出射方向$\omega_o$成立。

物理意义：反射的总能量不能超过入射能量。

实时渲染简化：
1. 预计算积分表（如UE4的环境BRDF查找表）
2. 使用归一化的解析BRDF（如归一化Phong）
3. 分离漫反射和镜面反射，分别确保守恒

例如，对于简化的BRDF：
$$f_r = \frac{k_d}{\pi} + k_s \cdot f_{specular}$$

确保$k_d + k_s \leq 1$（菲涅尔项会进一步调整）。

</details>

## 常见陷阱与错误

### 向量与点的混淆
- **错误**：将向量当作点进行平移变换
- **正确**：向量的齐次坐标第四分量为0，不受平移影响
- **调试技巧**：检查齐次坐标的w分量

### 矩阵乘法顺序
- **错误**：`translate * rotate * scale`（错误顺序）
- **正确**：`translate * rotate * scale * vertex`（从右向左）
- **记忆方法**：变换从右向左应用，最右边的先执行

### 四元数归一化
- **错误**：忘记归一化导致缩放
- **正确**：每次运算后重新归一化
- **优化**：使用平方根的快速近似算法

### 叉积的手性
- **错误**：混淆左手系和右手系
- **正确**：OpenGL用右手系，DirectX传统上用左手系
- **验证**：$\mathbf{x} \times \mathbf{y} = \mathbf{z}$（右手系）

### 浮点精度问题
- **错误**：直接比较浮点数相等
- **正确**：使用epsilon容差：`abs(a - b) < epsilon`
- **选择epsilon**：通常用`1e-6`，但需根据场景调整

### Gimbal Lock（万向锁）
- **症状**：欧拉角在某些角度失去一个自由度
- **解决**：使用四元数或轴角表示
- **转换时机**：仅在最终输出时转换为欧拉角

## 最佳实践检查清单

### 数学库设计
- [ ] 向量和点使用不同的类型或明确区分
- [ ] 提供行主序和列主序的明确说明
- [ ] 实现SIMD优化的向量运算
- [ ] 包含数值稳定的归一化函数
- [ ] 提供便捷的调试输出功能

### 变换系统
- [ ] 使用齐次坐标统一处理变换
- [ ] 缓存常用的变换矩阵
- [ ] 提供矩阵分解功能（提取旋转、缩放等）
- [ ] 实现视锥体裁剪的快速路径
- [ ] 支持不同坐标系之间的转换

### 四元数使用
- [ ] 始终保持四元数归一化
- [ ] 实现SLERP和NLERP两种插值
- [ ] 处理四元数点积为负的情况
- [ ] 提供与欧拉角、轴角的转换函数
- [ ] 考虑使用对偶四元数处理刚体变换

### 性能优化
- [ ] 避免不必要的矩阵求逆
- [ ] 使用快速平方根倒数近似
- [ ] 预计算能复用的中间结果
- [ ] 考虑内存对齐和缓存友好性
- [ ] 使用查表法加速三角函数

### 数值稳定性
- [ ] 使用稳定的正交化算法（如改进的Gram-Schmidt）
- [ ] 避免在接近奇异的情况下求逆
- [ ] 使用适当的epsilon进行浮点比较
- [ ] 定期重新正交化旋转矩阵
- [ ] 考虑使用双精度进行关键计算