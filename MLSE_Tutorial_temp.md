# 最大似然序列估计 (MLSE) 完整教程

---

## 目录

1. [核心思想](#核心思想)
2. [数学原理](#数学原理)
3. [Viterbi 算法详解](#viterbi-算法详解)
4. [Python 完整实现](#python-完整实现)
5. [应用场景](#应用场景)
6. [与 MMSE 的区别](#与-mmse-的区别)
7. [扩展阅读](#扩展阅读)

---

## 核心思想

### 问题定义

在数字通信系统中，发射机发送一个符号序列 **s** = (s₀, s₁, ..., sₙ₋₁)，经过信道后，接收机观测到的信号 **r** 受到噪声干扰：

$$\mathbf{r} = H \cdot \mathbf{s} + \mathbf{n}$$

其中：
- **s** 是发送的符号序列（离散取值，如 BPSK 的 ±1）
- **H** 是信道冲激响应矩阵
- **n** 是加性高斯白噪声 (AWGN)

### MLSE 的目标

> **在所有可能的符号序列中，找到使接收信号似然概率最大的那个序列。**

数学表达：

$$\hat{\mathbf{s}} = \arg\max_{\mathbf{s} \in \mathcal{S}^N} P(\mathbf{r} | \mathbf{s})$$

在 AWGN 信道下，这等价于 **最小欧几里得距离** 问题：

$$\hat{\mathbf{s}} = \arg\min_{\mathbf{s} \in \mathcal{S}^N} \|\mathbf{r} - H\mathbf{s}\|^2$$

---

## 数学原理

### 1. 似然函数推导

假设噪声 n 是零均值复高斯噪声，协方差矩阵 σ²I，则条件概率密度：

$$P(\mathbf{r} | \mathbf{s}) = \frac{1}{(\pi\sigma^2)^N} \exp\left(-\frac{\|\mathbf{r} - H\mathbf{s}\|^2}{\sigma^2}\right)$$

取对数（对数似然），最大化似然等价于最小化指数部分：

$$\hat{\mathbf{s}} = \arg\min_{\mathbf{s}} \sum_{k=0}^{N-1} |r_k - \sum_{l=0}^{L-1} h_l s_{k-l}|^2$$

### 2. 路径度量（Branch Metric）

定义每个时刻 k 的分支度量：

$$\gamma_k(s_k, s_{k+1}) = |r_k - \sum_{l=0}^{L-1} h_l s_{k-l}|^2$$

其中 L 是信道记忆长度。

### 3. 累计路径度量

$$\Gamma(\mathbf{s}) = \sum_{k=0}^{N-1} \gamma_k(s_k, s_{k+1})$$

MLSE 就是寻找使累计路径度量最小的完整序列。

---

## Viterbi 算法详解

MLSE 的**精确最优解**可以通过 **Viterbi 算法** 高效计算，避免了穷举所有可能的序列。

### 关键概念：Trellis（网格图）

| 参数 | 含义 |
|------|------|
| M | 调制阶数（如 BPSK: M=2, QPSK: M=4） |
| L | 信道记忆长度（ISI 跨度） |
| 状态数 | M^(L-1) |

**BPSK + 3径信道 (L=3)** → 状态数 = 2² = 4 个状态

### Viterbi 算法步骤

```
1. 初始化：设置 t=0 时刻各状态的累计度量 = 0
2. 对每个时刻 t = 1, 2, ..., N:
   a. 对每个状态，计算所有进入该状态的分支度量
   b. 对每个状态，选择累计度量最小的路径（幸存者）
   c. 记录幸存路径和新的累计度量
3. 终止：选择最终累计度量最小的状态
4. 回溯：沿幸存路径反向追溯，得到最优序列
```

### 复杂度对比

| 方法 | 计算复杂度 | 说明 |
|------|----------|------|
| 穷举搜索 | O(M^N) | 指数级，不可行 |
| **Viterbi (MLSE)** | **O(N · M^L)** | 线性于 N，指数于 L |
| 线性均衡 (ZF/MMSE) | O(N) | 次优，但计算简单 |

> ⚠️ **注意**：Viterbi 的复杂度随信道长度 L 指数增长，因此 L 通常限制在 ≤ 5。

---

## Python 完整实现

见同目录下的 `mlse_viterbi_demo.py` 文件，包含：

1. **Viterbi 解调器核心类**
2. **BPSK 调制的 MLSE 检测**
3. **与线性均衡器 (MMSE) 的 BER 对比**
4. **可视化网格图与误码率曲线**

---

## 应用场景

### 1. GSM 系统
GSM 使用 GMSK 调制，在时延扩展信道中，接收机使用 Viterbi 均衡器实现 MLSE，这是 GSM 标准定义的核心技术。

### 2. 磁盘存储读取
硬盘读通道中的 PRML (Partial Response Maximum Likelihood) 检测，本质上就是 MLSE。

### 3. 卷积码译码
卷积码的软判决 Viterbi 译码，与 MLSE 在数学上是同一个问题。

### 4. Wi-Fi / 5G 中的短包检测
在某些短块传输和已知信道条件下，MLSE 仍是性能标杆。

---

## 与 MMSE 的区别

| 特性 | MLSE | MMSE 线性均衡 |
|------|------|--------------|
| **最优性** | 符号序列意义下的最优 | 符号逐点最优（非全局） |
| **噪声放大** | 无 | 存在（零 forcing 更严重） |
| **计算复杂度** | 高 (M^L) | 低 (O(N)) |
| **需要信道信息** | 是 | 是 |
| **对深衰落的鲁棒性** | 强 | 弱 |
| **实现难度** | 中等 | 简单 |

> **通俗比喻**：
> - MMSE 就像 "逐字纠正错别字"（每次只看一个）
> - MLSE 就像 "通读全文后再判断哪个词更合理"（考虑上下文依赖）

---

## 扩展阅读

1. **G. D. Forney Jr.** — "The Viterbi Algorithm" (1973) — 经典论文
2. **Proakis & Salehi** — *Digital Communications* — 第 10 章 自适应均衡
3. **J. R. Barry, E. A. Lee, D. G. Messerschmitt** — *Digital Communication* — 第 6.4 节 MLSE
4. **S. M. Kay** — *Fundamentals of Statistical Signal Processing*

---

*文档生成时间：见代码运行输出*
