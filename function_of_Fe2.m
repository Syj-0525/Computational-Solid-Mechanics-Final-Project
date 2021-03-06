function [Fe1, Fe2] = My_function_of_Fe(JXYe,JBPe)
%%%%%%% 初始化Fe1和Fe2
Fe1 = zeros(9, 1);
Fe2 = zeros(9, 1);
%%%%%%% 初始化Fe1和Fe2

%%%%%%% 高斯积分点和权重
gp = [0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.9324695142031521, -0.6612093864662645, -0.2386191860831969];
gw = [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704, 0.3607615730481386, 0.4679139345726910];
kesi = gp;
ita = gp;
%%%%%%% 高斯积分点和权重

%%%%%%% (压强)自然边界条件数据
P_side = JBPe(1, 2);    % 边界单元的边界编号
l = JBPe(1, 3);         % 边界边外法向与x轴余弦
m = JBPe(1, 4);         % 边界边外法向与y轴余弦
p_value1 = JBPe(1, 5);  % 边界边第一点处压力值
P_value2 = JBPe(1, 6);  % 边界边第二点处压力值
%%%%%%% (压强)自然边界条件数据

%%%%%%% Fe1, Fe2的数值积分
%%%%%%% 当压力施加于单元第1边时，Fe1和Fe2的计算
if P_side == 1
    for i = 1:6
        fy = [1/4 * kesi(i) * (-1) * (kesi(i) - 1) * ((-1) - 1);
            1/2 * (-1) * (1 - kesi(i)^2) * ((-1) - 1);
            1/4 * kesi(i) * (-1) * (kesi(i) + 1) * ((-1) - 1);
            1/2 * kesi(i) * ( kesi(i) - 1) * (1 - (-1)^2);
            (1 - kesi(i)^2) * (1 - (-1)^2);
            1/2 * kesi(i) * (kesi(i) + 1) * (1 - (-1)^2);
            1/4 * kesi(i) * (-1) * (kesi(i) - 1) * ((-1) + 1);
            1/2 * (-1) * (1 - kesi(i)^2) * ((-1) + 1);
            1/4 * kesi(i) * (-1) * (kesi(i) + 1) * ((-1) + 1);];
        Fyp = [1/4 * (1 - kesi(i)) * (1 - (-1))
            1/4 * (1 + kesi(i)) * (1 - (-1))
            1/4 * (1 + kesi(i)) * (1 + (-1))
            1/4 * (1 - kesi(i)) * (1 + (-1))];
        P = [p_value1, P_value2, 0, 0]';

        %%%%%%% ita = -1时，速度插值函数对kesi的导数
        fy_kesi = [1/4 * (-1) * (kesi(i) - 1) * ((-1) - 1) + 1/4 * kesi(i) * (-1) * ((-1) - 1)
            -(-1) * kesi(i) * ((-1) - 1 )
            1/4 * (-1) * (kesi(i) + 1) * ((-1)- 1) + 1/4 * kesi(i) * (-1) * ((-1) - 1)
            1/2 * (kesi(i) - 1) * (1 - (-1)^2) + 1/2 * kesi(i) * (1 - (-1)^2)
            -2 * kesi(i) * (1 - (-1)^2)
            1/2 * (kesi(i) + 1) * (1 - (-1)^2) + 1/2 * kesi(i) * (1 - (-1)^2)
            1/4 * (-1) * (kesi(i) - 1) * ((-1) + 1) + 1/4 * kesi(i) * (-1) * ((-1) + 1)
            -(-1) * kesi(i) * ((-1) + 1 )
            1/4 * (-1) * (kesi(i) + 1) * ((-1) + 1) + 1/4 * kesi(i) * (-1) * ((-1) + 1)];  
        %%%%%%% ita = -1时，速度插值函数对kesi的导数

        %%%%%% Jacobi相关计算
        dx_dkesi = fy_kesi' * JXYe(:,1);
        dy_dkesi = fy_kesi' * JXYe(:,2);
        det_Jacobi = sqrt(dx_dkesi^2 + dy_dkesi^2);
        %%%%%% Jacobi相关计算

        Fe1 = Fe1 + gw(i) * Fyp' * P * l * fy * det_Jacobi;
        Fe2 = Fe2 + gw(i) * Fyp' * P * m * fy * det_Jacobi;
    end
end
%%%%%%% 当压力施加于单元第1边时，Fe1和Fe2的计算

%%%%%%% 当压力施加于单元第2边时，Fe1和Fe2的计算
if P_side == 2
    for j = 1:6
        fy = [1/4 * 1 * ita(j) * (1 - 1) * (ita(j) - 1);
            1/2 * ita(j) * (1 - 1^2) * (ita(j) - 1);
            1/4 * 1 * ita(j) * (1 + 1) * (ita(j) - 1);
            1/2 * 1 * (1 - 1) * (1 - ita(j)^2);
            (1 - 1^2) * (1 - ita(j)^2);
            1/2 * 1 * (1 + 1) * (1 - ita(j)^2);
            1/4 * 1 * ita(j) * (1 - 1) * (ita(j) + 1);
            1/2 * ita(j) * (1 - 1^2) * (ita(j) + 1);
            1/4 * 1 * ita(j) * (1 + 1) * (ita(j) + 1);];
        Fyp = [1/4 * (1 - 1) * (1 - ita(j))
            1/4 * (1 + 1) * (1 - ita(j))
            1/4 * (1 + 1) * (1 + ita(j))
            1/4 * (1 - 1) * (1 + ita(j))];
        P = [0, p_value1, P_value2, 0]';

        %%%%%%% kesi = 1时，速度插值函数对ita的导数
        fy_ita = [1/4 * 1 * (1 - 1) * (ita(j) - 1) + 1/4 * 1 * ita(j) * (1 - 1)
            1/2 * (1 - 1^2) * (ita(j) - 1) + 1/2 * ita(j) * (1 - 1^2)
            1/4 * 1 * (1 + 1) * (ita(j) - 1) + 1/4 * 1 * ita(j) * (1 + 1)
            - 1 * (1 - 1) * ita(j)
            -2 * (1 - 1^2) * (ita(j))
            - 1 * (1 + 1) * ita(j)
            1/4 * 1 * (1 - 1) * (ita(j) + 1) + 1/4 * 1 * ita(j) * (1 - 1)
            1/2 * (1 - 1^2) * (ita(j) + 1) + 1/2 * (1 - 1^2) * ita(j)
            1/4 * 1 * (1 + 1) * (ita(j) + 1) + 1/4 * 1 * ita(j) * (1 + 1)]; 
        %%%%%%% kesi = 1时，速度插值函数对ita的导数

        %%%%%% Jacobi相关计算
        dx_dita = fy_ita' * JXYe(:,1);
        dy_dita = fy_ita' * JXYe(:,2);
        det_Jacobi = sqrt(dx_dita^2 + dy_dita^2);
        %%%%%% Jacobi相关计算

        Fe1 = Fe1 + gw(j) * Fyp' * P * l * fy * det_Jacobi;
        Fe2 = Fe2 + gw(j) * Fyp' * P * m * fy * det_Jacobi;
    end
end
%%%%%%% 当压力施加于单元第2边时，Fe1和Fe2的计算

%%%%%%% 当压力施加于单元第3边时，Fe1和Fe2的计算
if P_side == 3
    for i = 1:6
        fy = [1/4 * kesi(i) * 1 * (kesi(i) - 1) * (1 - 1);
            1/2 * 1 * (1 - kesi(i)^2) * (1 - 1);
            1/4 * kesi(i) * 1 * (kesi(i) + 1) * (1 - 1);
            1/2 * kesi(i) * ( kesi(i) - 1) * (1 - 1^2);
            (1 - kesi(i)^2) * (1 - 1^2);
            1/2 * kesi(i) * (kesi(i) + 1) * (1 - 1^2);
            1/4 * kesi(i) * 1 * (kesi(i) - 1) * (1 + 1);
            1/2 * 1 * (1 - kesi(i)^2) * (1 + 1);
            1/4 * kesi(i) * 1 * (kesi(i) + 1) * (1 + 1);];
        Fyp = [1/4 * (1 - kesi(i)) * (1 - 1)
            1/4 * (1 + kesi(i)) * (1 - 1)
            1/4 * (1 + kesi(i)) * (1 + 1)
            1/4 * (1 - kesi(i)) * (1 + 1)];
        P = [0, 0, p_value1, P_value2]';

        %%%%%%% ita = 1时，速度插值函数对kesi的导数
        fy_kesi = [1/4 * 1 * (kesi(i) - 1) * (1 - 1) + 1/4 * kesi(i) * 1 * (1 - 1)
            -1 * kesi(i) * (1 - 1 )
            1/4 * 1 * (kesi(i) + 1) * (1 - 1) + 1/4 * kesi(i) * 1 * (1 - 1)
            1/2 * (kesi(i) - 1) * (1 - 1^2) + 1/2 * kesi(i) * (1 - 1^2)
            -2 * kesi(i) * (1 - 1^2)
            1/2 * (kesi(i) + 1) * (1 - 1^2) + 1/2 * kesi(i) * (1 - 1^2)
            1/4 * 1 * (kesi(i) - 1) * (1 + 1) + 1/4 * kesi(i) * 1 * (1 + 1)
            -1 * kesi(i) * (1 + 1 )
            1/4 * 1 * (kesi(i) + 1) * (1 + 1) + 1/4 * kesi(i) * 1 * (1 + 1)]; 
        %%%%%%% ita = 1时，速度插值函数对kesi的导数

        %%%%%% Jacobi相关计算
        dx_dkesi = fy_kesi' * JXYe(:,1);
        dy_dkesi = fy_kesi' * JXYe(:,2);
        det_Jacobi = sqrt(dx_dkesi^2 + dy_dkesi^2);
        %%%%%% Jacobi相关计算

        Fe1 = Fe1 + gw(i) * Fyp' * P * l * fy * det_Jacobi;
        Fe2 = Fe2 + gw(i) * Fyp' * P * m * fy * det_Jacobi;
    end
end
%%%%%%% 当压力施加于单元第3边时，Fe1和Fe2的计算

%%%%%%% 当压力施加于单元第4边时，Fe1和Fe2的计算
if P_side == 4
    for j = 1:6
        fy = [1/4 * (-1) * ita(j) * ((-1) - 1) * (ita(j) - 1);
            1/2 * ita(j) * (1 - (-1)^2) * (ita(j) - 1);
            1/4 * (-1) * ita(j) * ((-1) + 1) * (ita(j) - 1);
            1/2 * (-1) * ((-1) - 1) * (1 - ita(j)^2);
            (1 - (-1)^2) * (1 - ita(j)^2);
            1/2 * (-1) * ((-1) + 1) * (1 - ita(j)^2);
            1/4 * (-1) * ita(j) * ((-1) - 1) * (ita(j) + 1);
            1/2 * ita(j) * (1 - (-1)^2) * (ita(j) + 1);
            1/4 * (-1) * ita(j) * ((-1) + 1) * (ita(j) + 1);];
        Fyp = [1/4 * (1 - (-1)) * (1 - ita(j))
            1/4 * (1 + (-1)) * (1 - ita(j))
            1/4 * (1 + (-1)) * (1 + ita(j))
            1/4 * (1 - (-1)) * (1 + ita(j))];
        P = [p_value1, 0, 0, P_value2]';

        %%%%%%% kesi = -1时，速度插值函数对kesi的导数
        fy_ita = [1/4 * (-1) * ((-1) - 1) * (ita(j) - 1) + 1/4 * (-1) * ita(j) * ((-1) - 1)
            1/2 * (1 - (-1)^2) * (ita(j) - 1) + 1/2 * ita(j) * (1 - (-1)^2)
            1/4 * (-1) * ((-1) + 1) * (ita(j) - 1) + 1/4 * (-1) * ita(j) * ((-1) + 1)
            - (-1) * ((-1) - 1) * ita(j)
            -2 * (1 - (-1)^2) * (ita(j))
            - (-1) * ((-1) + 1) * ita(j)
            1/4 * (-1) * ((-1) - 1) * (ita(j) + 1) + 1/4 * (-1) * ita(j) * ((-1) - 1)
            1/2 * (1 - (-1)^2) * (ita(j) + 1) + 1/2 * (1 - (-1)^2) * ita(j)
            1/4 * (-1) * ((-1) + 1) * (ita(j) + 1) + 1/4 * (-1) * ita(j) * ((-1) + 1)];   
        %%%%%%% kesi = -1时，速度插值函数对kesi的导数

        %%%%%% Jacobi相关计算
        dx_dita = fy_ita' * JXYe(:,1);
        dy_dita = fy_ita' * JXYe(:,2);
        det_Jacobi = sqrt(dx_dita^2 + dy_dita^2);
        %%%%%% Jacobi相关计算

        Fe1 = Fe1 + gw(j) * Fyp' * P * l * fy * det_Jacobi;
        Fe2 = Fe2 + gw(j) * Fyp' * P * m * fy * det_Jacobi;
    end
end
%%%%%%% 当压力施加于单元第4边时，Fe1和Fe2的计算
%%%%%%% Fe1, Fe2的数值积分