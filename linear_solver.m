clc
clear
load matrix

%%%%%% 求解线性方程
x = K\B;   % 高斯消元
ux_k_1 = x(1:Nz);
vy_k_1 = x(1 + Nz: 2 * Nz);
p4 = x(1 + 2 * Nz: 2 * Nz + Nd);
p_k_1= [Pding2Pzong(p4, JMV)]';   % 压力插值计算
%%%%%% 求解线性方程

%%%%%% 输出Tecplot后处理结果
E4 = E * 4;
Nz = Nz;
data = [JXYV, ux_k_1, vy_k_1, sqrt(ux_k_1.^2 + vy_k_1.^2), p_k_1];
JMV4 = JMV_9to4(JMV);    % 二次九结点单元拆分为四个线性单元, 并建立IEN


k = 1; 
for i = 1:length(JMV4(:, 1))
    for j = 1: 4
        %%%%%% 生成Tecplot后处理数据
        JXY_post(k, 1) = JMV4(i, j);
        JXY_post(k, 2) = JXYV(JMV4(i, j), 1);
        JXY_post(k, 3) = JXYV(JMV4(i, j), 2);

        ux_post(k, 1) = ux_k_1(JMV4(i, j));
        vy_post(k, 1) = vy_k_1(JMV4(i, j));
        p_post(k, 1) = p_k_1(JMV4(i, j));
        %%%%%% 生成Tecplot后处理数据

        k = k + 1;
    end
end
clear k;

post_data = [JXY_post, ux_post, vy_post, sqrt(ux_post.^2 + vy_post.^2), p_post];
post_node = [1:1:size(post_data,1)];
post_node2 = reshape(post_node, [E4/16, 64])';   % 当post_node行太长，对其换行
%%%%%% 输出Tecplot后处理结果

clear B E II JBV JMV JXYZ K N_matrix Nd Nz   % 清除多余变量
clear data viscosity p4 u v x i j k               % 清除多余变量
save result_of_n1                                 % 存储结果