clc
clear
load matrix

%%%%%% 求解线性方程
x = K\B;
ux_k_1 = x(1:Nz);
vy_k_1 = x(1 + Nz: 2 * Nz);
p4 = x(1 + 2 * Nz: 2 * Nz + Nd);
p_k_1= [Pding2Pzong(p4, JMV)]';   % 压力插值计算
%%%%%% 求解线性方程

%%%%%% 输出Tecplot后处理结果
E = E * 4;
Nz = Nz;
data = [JXYV, ux_k_1, vy_k_1, sqrt(ux_k_1.^2 + vy_k_1.^2), p_k_1];
JMV4 = JMV_9to4(JMV);    % 二次九结点单元拆分为四个线性单元, 并建立IEN
%%%%%% 输出Tecplot后处理结果

clear B E II JBV JMV JMV4 JXYZ K N_matrix Nd Nz   % 清除多余变量
clear data viscosity p4 u v x                     % 清除多余变量
save result_of_n1                                 % 存储结果