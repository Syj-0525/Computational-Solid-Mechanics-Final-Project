function [Be1, Be2] = function_of_Be(e_JXY)
%%%%%%% 初始化Be1和Be2
Be1 = zeros(4, 9);
Be2 = zeros(4, 9);
%%%%%%% 初始化Be1和Be2

%%%%%%% 高斯积分点和权重
gp = [0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.9324695142031521, -0.6612093864662645, -0.2386191860831969];
gw = [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704, 0.3607615730481386, 0.4679139345726910];
kesi = gp;
ita = gp;
%%%%%%% 高斯积分点和权重

%%%%%%% Be1和Be2的数值积分
for i = 1:6
    for j = 1:6
        %%%%%%% 压力插值函数
        Fyp = [1/4 * (1-kesi(i)) * (1 - ita(j))
            1/4 * (1 + kesi(i)) * (1 - ita(j))
            1/4 * (1 + kesi(i)) * (1 + ita(j))
            1/4 * (1 - kesi(i)) * (1 + ita(j))];
        %%%%%%% 压力插值函数

        %%%%%%% 速度插值函数对kesi和ita导数
        fy_kesi = [1/4 * ita(j) * (kesi(i) - 1) * (ita(j) - 1) + 1/4 * kesi(i) * ita(j) * (ita(j) - 1)
            -ita(j) * kesi(i) * (ita(j) - 1 )
            1/4 * ita(j) * (kesi(i) + 1) * (ita(j) - 1) + 1/4 * kesi(i) * ita(j) * (ita(j) - 1)
            1/2 * (kesi(i) - 1) * (1 - ita(j)^2) + 1/2 * kesi(i) * (1 - ita(j)^2)
            -2 * kesi(i) * (1 - ita(j)^2)
            1/2 * (kesi(i) + 1) * (1 - ita(j)^2) + 1/2 * kesi(i) * (1 - ita(j)^2)
            1/4 * ita(j) * (kesi(i) - 1) * (ita(j) + 1) + 1/4 * kesi(i) * ita(j) * (ita(j) + 1)
            -ita(j) * kesi(i) * (ita(j) + 1 )
            1/4 * ita(j) * (kesi(i) + 1) * (ita(j) + 1) + 1/4 * kesi(i) * ita(j) * (ita(j) + 1)];
        fy_ita = [1/4 * kesi(i) * (kesi(i) - 1) * (ita(j) - 1) + 1/4 * kesi(i) * ita(j) * (kesi(i) - 1)
            1/2 * (1 - kesi(i)^2) * (ita(j) - 1) + 1/2 * ita(j) * (1 - kesi(i)^2)
            1/4 * kesi(i) * (kesi(i) + 1) * (ita(j) - 1) + 1/4 * kesi(i) * ita(j) * (kesi(i) + 1)
            - kesi(i) * (kesi(i) - 1) * ita(j)
            -2 * (1 - kesi(i)^2) * (ita(j))
            - kesi(i) * (kesi(i) + 1) * ita(j)
            1/4 * kesi(i) * (kesi(i) - 1) * (ita(j) + 1) + 1/4 * kesi(i) * ita(j) * (kesi(i) - 1)
            1/2 * (1 - kesi(i)^2) * (ita(j) + 1) + 1/2 * (1 - kesi(i)^2) * ita(j)
            1/4 * kesi(i) * (kesi(i) + 1) * (ita(j) + 1) + 1/4 * kesi(i) * ita(j) * (kesi(i) + 1)];
        %%%%%%% 速度插值函数对kesi和ita导数

        %%%%%% Jacobi相关计算
        dx_dkesi = fy_kesi' * e_JXY(:,1);
        dx_dita = fy_ita' * e_JXY(:,1);
        dy_dkesi = fy_kesi' * e_JXY(:,2);
        dy_dita = fy_ita' * e_JXY(:,2);
        Jacobi = [ dx_dkesi dy_dkesi
            dx_dita dy_dita];
        AAAA = inv(Jacobi) * [fy_kesi'; fy_ita'];
        fy_x = AAAA(1,:)';
        fy_y = AAAA(2,:)';
        det_Jacobi = det(Jacobi);
        %%%%%% Jacobi相关计算

        %%%%%% Be1和Be2单元子块矩阵计算
        Be1 = Be1 + gw(i) * gw(j)  * Fyp * fy_x' * det_Jacobi;
        Be2 = Be2 + gw(i) * gw(j)  * Fyp * fy_y' * det_Jacobi;
        %%%%%% Be1和Be2单元子块矩阵计算
    end
end
%%%%%%% Be1和Be2的数值积分