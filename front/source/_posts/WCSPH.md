---
title: 《WCSPH》
date: 2022-12-10 21:14:08

categories:
  - R-图形学
  - 流体仿真
tags:
  - WCSPH

img: 'http://43.143.246.2/imgs/wcsph/featureImg.jpg'
mathjax: true
summary: 本章介绍弱可压缩流体仿真技术。
---

# 写在前面🕶️
WCSPH方法出自2007的文章<font color=blue>**《Weakly compressible SPH for free surface flows》**</font>，我将会在文末给出论文原文供大家参考学习。本章所介绍的内容主要源于论文阅读后进行的总结、思考和实践。欢迎大家交流讨论。

# 开始学习😊

## 1、SPH模型
WCSPH使用的N-S方程形式如下：
$$
  \frac{d\pmb{v}}{dt} = -\frac{1}{ρ}\triangledown P + f_{vis} + f_{ext}
$$
其中第一项为压强力贡献的加速度，第二项为粘滞力贡献的加速度，第三项为外力贡献的加速度。接下来分别探究每一项的计算。
### 1.1 外力部分
在模拟中，考虑 $f_{ext} = \pmb{g}$，$\pmb{g}$ 为重力加速度。
### 1.2 压力部分
压力计算基于泰特方程（Tait equation）：
$$
  P = B\big(\Big(\frac{ρ}{ρ_0}\Big)^γ - 1\big)
$$
其中**P**表示压力；**B**表示压力常数，用于控制密度波动: $ \frac{|\triangle ρ|}{ρ_0} $，其中 $ |\triangle ρ| = ρ-ρ_0 $，**B**的基本形式为 $ \frac{kρ_0}{γ} $，经过一些考虑，**B**取<font color=red> $ \frac{ρ_0c_s^2}{γ} $</font>，其中 $c_s$ 为流体中的声速；γ为常数，取7。  
  
压力计算依赖与局部粒子密度，因此通过场量近似计算方程：
<img src='http://43.143.246.2/imgs/wcsph/field_approximate.png' width = "310" height = "210" align=center />
得到密度估计： $ ρ_a = \sum_b m_bW_{ab} $。
进一步的，我们使用对称形式的方程来计算压力贡献的加速度：
$$
  \frac{d\pmb{v}}{dt} = -\frac{1}{ρ}\triangledown P = - \sum_b m_b\Big( \frac{P_a}{ρ_a^2} + \frac{P_b}{ρ_b^2} \Big)\bigtriangledown W_{ab}
$$
在得到密度和压强参数后，通过压强密度比的梯度可得上式，其推导如下：
<img src='http://43.143.246.2/imgs/wcsph/derivation of pressure_acc.png' width = "400" height = "500" align=center />
相比于在不可压缩条件下求解压力的泊松方程，使用理想气体方程求解压力虽然会导致可压缩性，但是使用Tait方程能够在高声速条件下，保持可接受的小密度偏差。

### 1.3 粘滞力部分
在论文中，粘滞力是通过人工方式添加的，其目的是为了保证模拟中的数值稳定性和冲击现象。其贡献的加速度计算如下：
<img src='http://43.143.246.2/imgs/wcsph/vis_acc.png' width = "350" height = "210" align=center />
其中 $v^T_{ab} x_{ab} > 0$ 等价于速度的散度大于0（$\triangledown·v > 0$)；$\prod_{ab}$计算如下：
$$
  \prod_{ab} = -μ\Big(\frac{v^T_{ab} x_{ab}}{|x_ab|^2 + εh^2}\Big)
$$
其中μ是粘性项，用 $ \frac{2αhc_s}{ρ_a + ρ_b} $ 计算，α 为粘度常数，论文中使用的范围在[0.08, 0.5]。${εh}^2$ 项是为了保证分母大于0，论文中 ε 取0.01。<font color=red>另外，粘性项（也叫扩散项）一般在N-S方程表示为 $ \eta\triangledown ^2 v $</font>，$\eta$ 表示动力粘度，$v$ 是速度场，实验中通常也可以直接进行人工拉普拉斯算子计算来得到粘滞力加速度。  
  
## 2、表面张力（surface tension）

