# moonbit-seismic: MoonBit 地震波信号处理库

面向地震物理、勘探数据和教学实验的处理工具库。本库完全采用 MoonBit 开发，提供地震信号读取、预处理、数字滤波、互相关匹配、到时拾取、震级估算和走时曲线计算等多维度功能。

A MoonBit library for seismic wave signal processing, tailored for seismic physics, exploration data, and pedagogical experiments.

---

## 核心功能 (Key Features)

- **SAC 二进制文件读取与写入 (SAC I/O)**
  - 支持 SAC (Seismic Analysis Code) 格式的双端字节序（Little-Endian/Big-Endian）二进制数据解析与序列化。
  - 自动识别并解析浮点、整型及字符头字段（包括采样率、起始时间、站名及通道名等），支持直接还原为内存结构体。
- **miniSEED 数据记录简化解析 (miniSEED Parser)**
  - 解析 Fixed Section of Data Header (FSDH) 字节。
  - 支持未压缩 32 位整型和浮点型数据段的解包与重构，按规范计算采样率。
- **地震信号数字滤波 (Digital Filtering)**
  - 基于双线性变换法实现 2 阶 IIR 巴特沃斯滤波器 (Butterworth Lowpass/Highpass/Bandpass/Bandstop)。
  - 提供**零相位滤波模式 (Zero-phase Filtering / filtfilt)**，通过前向和后向双向滤波彻底消除相位偏移，完美还原地震相到时。
  - 提供滑动窗口移动平均平滑滤波器 (Moving Average Smoothing)。
- **信号预处理与分析 (Preprocessing)**
  - 包含去均值 (Demean) 与去线性趋势 (Detrend，基于最小二乘线性回归)。
- **波形互相关与模板匹配 (Cross-Correlation)**
  - 极低时滞滑动窗口归一化互相关匹配 (Normalized Cross-Correlation, NCC)。
  - 支持多通道波形相似度匹配及最佳延迟时间 (Lag) 检测。
- **地震波到时拾取 (Arrival-Time Picking)**
  - 实现经典的 $O(N)$ 线性时间复杂度滑动窗口 STA/LTA（短动/长动均值比）触发算法。
  - 支持基于多级阈值的事件触发Onset（开端）和Offset（结束）提取。
- **震级自动估算 (Magnitude Estimation)**
  - 地方性震级 $M_L$ (Richter Local Magnitude) 计算（带有随距离校正校准）。
  - 体波震级 $M_b$ (Body Wave Magnitude) 与面波震级 $M_s$ (Surface Wave Magnitude) 振幅/周期比自动算式。
- **1D 走时曲线计算 (Travel-Time Curves)**
  - 构建多层 1D 速度模型。
  - 基于直射线与折射平均慢度路径积分（Path-Integral Slowness）高精度计算任意震源深度与震中距下的 P 波与 S 波旅行时。
  - 内置全球标准 `IASPEI91` 地球速度模型。
- **速度模型数值插值 (Velocity Interpolation)**
  - 一维线性插值与**自然三次样条插值 (Natural Cubic Spline Interpolation)**。
  - 保证层速度模型在任意深度的光滑连续过渡。

---

## 项目结构 (Project Layout)

```text
├── moon.mod                # 模块配置 (Module definition)
├── moon.pkg                # 根包导入配置 (Root package setup)
├── types.mbt               # 核心 Trace 结构体与预处理 helpers
├── sac.mbt                 # SAC 二进制数据 I/O 引擎
├── miniseed.mbt            # miniSEED 解包器与打包器
├── filter.mbt              # Butterworth 与平滑滤波器
├── correlation.mbt         # 归一化互相关与模板匹配
├── picking.mbt             # 滑动窗口 STA/LTA 与事件检测
├── magnitude.mbt           # ML / Mb / Ms 震级解算器
├── traveltime.mbt          # 1D 地壳速度模型与走时计算
├── interpolation.mbt       # 线性与三次样条插值器
├── seismic_test.mbt        # 全功能模块黑盒测试套件
└── cmd/
    └── main/
        ├── moon.pkg        # 可执行程序导入包配置
        └── main.mbt        # 可运行 CLI 演示程序 (Demo code)
```

---

## 快速上手 (Quick Start)

### 编译与测试 (Build & Test)

保证已安装最新版 MoonBit 工具链 (0.10.3 或更高版本)。

1. **格式化自查 (Code Formatting Check)**:
   ```bash
   moon fmt
   ```
2. **静态检查与分析 (Static Analysis Check)**:
   ```bash
   moon check --deny-warn
   ```
3. **执行测试套件 (Run Tests)**:
   ```bash
   moon test --deny-warn
   ```
4. **运行演示 CLI 程序 (Run Demo Program)**:
   ```bash
   moon run cmd/main
   ```

---

## 演示输出示例 (Demo Output)

运行 `moon run cmd/main` 后，您将看到如下地震学 data 处理流程的演示输出：

```text
==================================================
   MoonBit Seismic Wave Signal Processing Demo    
==================================================
1. Created synthetic Trace:
   Station: ST01, Channel: BHZ
   Samples: 200, Sampling Rate: 50 Hz
   Mean: 0.09323698847629626, Std Dev: 0.48135320964707284

2. Applied 2nd-order Lowpass Filter (Cutoff: 5.0 Hz, Zero-phase):
   Original Std Dev: 0.48135320964707284
   Filtered Std Dev: 0.4763126889788772

3. Ran STA/LTA Arrival Picker:
   Detected Event Arrival!
   Onset Time: 1.8 s
   Offset Time: 2.44 s

4. Earthquake Magnitude Estimation:
   Richter Local Magnitude (M_L): 1.1850538726071076

5. 1D Travel-Time Calculations (IASPEI91 model, source depth 15km, distance 80km):
   P-wave travel time: 14.033466031120444 s
   S-wave travel time: 24.22443541086265 s

6. Velocity Model Interpolation:
   Cubic Spline interpolated velocity at 15.0 km depth: 6.220394736842105 km/s

7. SAC File Serialization:
   Successfully serialized trace to SAC binary format (1432 bytes).
==================================================
```

---

## 许可证 (License)

本项目采用 OSI 认可的 **Apache-2.0** 许可证开源。详情见根目录 [LICENSE](file:///c:/Users/33046/Desktop/黑客松_moonbit/嫂子/LICENSE) 文件。
