---
title: 《PCISPH》
date: 2022-12-15 11:05:42

categories:
    - R-图形学
    - 流体仿真
tags:
    - PCISPH
    - research

img: 'http://43.143.246.2/imgs/pcisph/featureImg.png'
mathjax: true
summary: 本章介绍基于预测-矫正的不可压缩SPH技术。
---

# 写在前面🕶️
PCISPH最早出自2009年的论文<font color=blue>**《Predictive-Corrective Incompressible SPH》**</font>。该技术最显著的特点是采用预估-矫正的迭代方式来满足不可压缩性而非求解压力泊松方程。  
WCSPH技术有两个显著的缺点，一是在利用Tait方程计算压强时，为了满足密度偏差阈值，需要不断调整刚度参数$k$，该参数在未进行仿真以前，无法知其合适与否，因此需要动画师不断调参测试；二是其时间步长计算要求高，计算量大，会导致一部分的性能损耗。  
$$
  Tait:\quad P = \frac{k\rho_ 0}{\gamma}\big(\Big(\frac{\rho}{\rho_0}\Big)^\gamma - 1\big)
$$
PCISPH技术能够有效的解决相关问题，并且性能对比WCSPH有很大提升。  

# 开始学习😊
不可压缩的N-S方程如下：
$$
  \frac{d\pmb{v}}{dt} = -\frac{1}{\rho}\bigtriangledown P + \eta\bigtriangledown ^2 v + f_{ext}
$$
PCISPH实际上是在不断计算后一时刻的密度偏差的过程中，矫正上一次计算的压力，在满足条件时，也就是密度偏差小于等于阈值时，跳出迭代、导出当前所得的压力并应用到N-S方程中。

### 1. PCISPH模型
对于外力和粘滞力部分的计算可以沿用WCSPH的设计。
  
对于压力计算，这里还采用对称SPH公式：
$$
  \frac{d\pmb{v}}{dt} = -\frac{1}{\rho}\bigtriangledown P = - \sum_b m_b\Big( \frac{P_a}{\rho_a^2} + \frac{P_b}{\rho_b^2} \Big)\bigtriangledown W_{ab}
$$
  
**技术难点在于通过密度偏差来矫正压力。**首先给出密度估计方程（设粒子质量相同）：
$$
    \rho_a = \sum_b m_bW_{ab}
$$
$$
    =>\quad \rho_a = m\sum_b W_{ab}
$$
展开并带入时间参数得：
<img src='http://43.143.246.2/imgs/pcisph/estimate_density.png' width='450' />
其中 $\pmb d_{ab}(t) = x_a(t)-x_b(t)$，$\triangle\pmb d_{ab}(t) = \triangle x_a(t)-\triangle x_b(t)$。  
假设 $\triangle\pmb d_{ab}(t)$ 相对较小，那么可以对式中的核 $W$ 进行一阶泰勒展开：
<img src='http://43.143.246.2/imgs/pcisph/estimate_density_taylor_series.png' width='550' />
其中第二项其实由非压力部分和压力部分的贡献组成，展开第二项后，上式得：
<img src='http://43.143.246.2/imgs/pcisph/decomposing_equation.png' width='550' />
将非压力贡献的密度计算视为预测密度 $\rho ^*_a$，并且矫正的目的是令 $\rho_a(t+1)=\rho _0$，因此得：
<img src='http://43.143.246.2/imgs/pcisph/predict_density.png' width='400' />
在进行到这一步推导后，式中的加速度项已经可以通过对称SPH公式来计算，但是PCISPH引入了近似值和简化，使得最终只有一个未知压力值 $p_a$：
<img src='http://43.143.246.2/imgs/pcisph/acc_approx.png' width='300' />
将上式带入预测密度计算得：
<img src='http://43.143.246.2/imgs/pcisph/density_approx.png' width='500' />
用 $p_a \approx p_b$，并且近似计算 $\sum_k\bigtriangledown W(x_b-x_k) \approx \bigtriangledown W(x_b-x_a) $，就能继续化简得：
<img src='http://43.143.246.2/imgs/pcisph/k_PCI.png' width='750' />
其中 $\rho _0 - \rho ^*_a$ 可以用密度偏差 $\rho _{err}^{*}$ 表示。至此，我们已经推导出用密度偏差矫正压力的公式。

### 2. 算法流程

注：rho_err $=\rho _{err}^{*}$ ；eta 为阈值；p_error 为上式中的 $p_a$。

``` whatever
Loop_main():

    update_dt();
    for each particle:
        update_neighbors();
        init pressure p(t)=0.0;

    // PCI 循环
    Loop_PCI(rho_err > eta || iter < minIterations):
        for each particle:
            update_velocity_from_NonePressure();
            update_position();

        for each particle:
            predict_density();          // 通过密度估计方程预测密度
            compute_density_error();    // 计算 k_PCI, rho_error
            correct_pressure();         // p(t) += p_error

    for each particle:
        update_velocity_from_Pressure();
        update_postion();
        
```

### 3. 总结

（1）通过PC（predict-correct）技术实现I（Incompressible）条件
（2）可采用大时步进行模拟并减少了每个时步的计算量
（3）对比WCSPH技术，性能提升在一个数量级，并能保持与WCSPH模拟相当的结果
  
一个显而易见的问题在于，本技术使用了很多领域近似，这会导致预测的各项直接或间接物理量有所偏差。若进行这些问题的矫正，则需要额外的工作，如调整时步或重新计算邻居等。  
  
**以上便是对PCISPH的解读，希望能对你有所帮助呦～赞赞！赞赞！**
