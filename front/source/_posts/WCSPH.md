---
title: 《WCSPH》
date: 2022-12-10 21:14:08

categories:
  - R-图形学
  - 流体仿真
tags:
  - WCSPH
  - research

img: 'http://43.143.246.2/imgs/wcsph/featureImg.jpg'
mathjax: true
summary: 本章介绍弱可压缩流体仿真技术。
---

# 写在前面🕶️

WCSPH方法出自2007年的论文<font color=blue>**《Weakly compressible SPH for free surface flows》**</font>。本章所介绍的内容主要源于论文阅读后进行的总结、思考和实践。欢迎大家交流讨论。

# 开始学习😊

## 1、WCSPH模型

---
WCSPH使用的N-S方程形式如下：
$$
  \frac{d\pmb{v}}{dt} = -\frac{1}{\rho}\bigtriangledown P + f_{vis} + f_{ext}
$$
其中第一项为压强力贡献的加速度，第二项为粘滞力贡献的加速度，第三项为外力贡献的加速度。接下来分别探究每一项的计算。  

### 1.1 压力部分

压力计算基于泰特方程（Tait equation）：
$$
  P = B\big(\Big(\frac{\rho}{\rho_0}\Big)^\gamma - 1\big)
$$
其中**P**表示压力；**B**表示压力常数，用于控制密度波动: $ \frac{|\triangle \rho|}{\rho_0} $，其中 $ |\triangle \rho| = \rho-\rho_0 $，**B**的基本形式为 $ \frac{k\rho_ 0}{\gamma} $，其中 $k$ 是刚度参数(stiffness)，表示流体抵抗形变的能力，经过一些考虑，**B**取<font color=red> $ \frac{\rho_0\ c_s^2}{\gamma} $</font>，其中 $c_s$ 为流体中的声速；$\gamma$ 为常数，论文中取7。  
压力计算依赖与局部粒子密度，因此通过场量近似计算方程：
$$
  A(x) = \sum_{j\in N(x)} m_j\frac{A_j}{\rho_j}W(x-x_j)
$$
得到密度估计： $ \rho_a = \sum_b m_bW_{ab} $。
进一步的，我们使用对称SPH公式来计算压力贡献的加速度：
$$
  \frac{d\pmb{v}}{dt} = -\frac{1}{\rho}\bigtriangledown P = - \sum_b m_b\Big( \frac{P_a}{\rho_a^2} + \frac{P_b}{\rho_b^2} \Big)\bigtriangledown W_{ab}
$$
在得到密度和压强参数后，通过压强密度比的梯度可得上式，其推导如下：
<img src='http://43.143.246.2/imgs/wcsph/derivation of pressure_acc.png' width = "400" height = "500" align=center />  
相比于在不可压缩条件下求解压力的泊松方程，使用理想气体方程求解压力虽然会导致可压缩性，但是使用Tait方程能够在高声速条件下，保持可接受的小密度偏差。

### 1.2 粘滞力部分

论文中，粘滞力是通过人工方式添加的，其目的是为了保证模拟中的数值稳定性和冲击现象。其贡献的加速度计算如下：
<img src='http://43.143.246.2/imgs/wcsph/vis_acc.png' width = "350" height = "210" align=center />
其中 $v^T_{ab} x_{ab} > 0$ 等价于速度的散度大于0（$\bigtriangledown \cdot v > 0$)； $\prod_{ab}$计算如下：
$$
  \prod_{ab} = -\mu\Big(\frac{v^T_{ab} x_{ab}}{|x_ab|^2 + \epsilon h^2}\Big)
$$
其中 $\mu$ 是粘性项，用 $ \frac{2\alpha hc_s}{\rho_a + \rho_b} $ 计算，$\alpha$ 为粘度常数，论文中使用的范围在[0.08, 0.5]；$h$ 是粒子半径。${\epsilon h}^2$ 项是为了保证分母恒正，论文中 $\epsilon$ 取0.01。另外，<font color=red>粘性项（也叫扩散项）一般在N-S方程表示为 $ \eta\bigtriangledown ^2 v $</font>，$\eta$ 表示动力粘度，$v$ 是速度场，实验时也可以直接进行人工拉普拉斯算子计算来得到粘滞力加速度。  

### 1.3 外力部分

在模拟中，考虑 $f_{ext} = \pmb{g}$，$\pmb{g}$ 为重力加速度。  

## 2、表面张力（Surface Tension）

---
对于流体分子，因为内部电荷不均匀分布而产生了分子间作用力。流体内部的分子受到的各个方向的力相互平衡，但是表面的分子受到的分子间作用力并不平衡。例如，在水和空气的交界面上，表面水分子受到的来自空气分子的作用力相比于来自内部水分子的力要弱得多，结果就产生了一个倾向于把表面水分子拉向内部的净作用力，现象如水珠。  
  
论文提出了一种基于内聚力的表面张力模型。该方法使用核函数来对吸引力加权：
$$
  \frac{d\pmb v_a}{dt} = -\frac{\kappa}{m_a}\sum_b m_bW_{ab}
$$
其中，$\kappa$ 为表面张力，是材料的特定常数，0.01接近于水的参数。  
  
在零重力情境下，初始化的流体立方体速度为0并且没有加速度，那么由于表面张力，其应内聚为球形：
![(左)初始化流体立方体; (中)之前的张力计算模型效果;  (右)本文的模型效果](http://43.143.246.2/imgs/wcsph/surface_tension_compare.png)
<!-- <img src='http://43.143.246.2/imgs/wcsph/surface_tension_compare.png' width = "350" height = "210" align=center /> -->

## 3、时间步长（Time Step）

---
时间步长的选择对数值稳定性至关重要，论文中使用的时步计算由CFL条件、粘性项和外力决定：
$$
  \triangle t = \min\bigg(0.25\cdot\min\Big(\frac{h}{|f_a|}\Big), 0.4\cdot\frac{h}{c_s\cdot(1+0.6\alpha)}\bigg)
$$
其中 $f_a$ 表示外力；$\alpha$ 为常数。

## 4、边界条件（Boundary Condition）
---
设流体粒子a，边界粒子k，二者碰撞时施加的力 $f_{ak}$ 计算如下：
$$
  f_{ak} = \frac{m_k}{m_a + m_k} \Gamma(x_a, x_b) \frac{x_a - x_k}{|x_a - x_k|}
$$
其中 $\Gamma$ 定义如下：
<img src='http://43.143.246.2/imgs/wcsph/function_Gamma.png' width = "350" height = "210" align=center />
$\quad\quad$其中 $q = \frac{|x_a-x_b|}{h}$。

## 5、算法流程

``` whatever
Loop_main:

    update_dt();
    for each particle:
      update_neighbors();
  
    for each particle:
        compute_density();
        compute_pressure();
  
    for each particle:
        update_velocity_from_NonePressure();  // viscous force, external force...

    for each particle:
        update_velocity_from_Pressure();
        update_postion();

```

## 6、总结

（1）提出了一种基于Tait方程的弱可压缩SPH方法，模拟对象为低粘牛顿流体
（2）提出了一种改进的基于内聚力的表面张力模型，可处理单相自由表面流动中的高曲率问题
（3）采用显式时间积分方案($x_{t+1}=x_t+v_t\cdot dt$)，除了要求高音速导致的小时步外，这也是需要小时步的原因  
  
**以上便是对wcsph的解读，本文未涉及的部分还有邻居搜索等模块，后续可能会出相关内容。如果这篇文章对你有帮助，请在篇头点个赞赞呦～**
