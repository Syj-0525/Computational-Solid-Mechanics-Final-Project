clc
clear
%%%%%% 读取网格数据
load msh
%%%%%% 读取网格数据

%%%%%% 设定流体黏度
viscosity = 1000;
%%%%%% 设定流体黏度

%%%%%% 设定边界条件
u1 = 0; v1 = 0;     % 设定1边界速度
u3 = 0.01; v3 = 0;  % 设定3边界速度

%%% Case1 槽道流
%%%%% 编写边界条件数据JBP与JBV
% JBV1 =[ BP1', u1 * ones(size(BP1))', v1 * ones(size(BP1))'];
% JBV3 =[ BP3', u3 * ones(size(BP3))', v3 * ones(size(BP3))'];
% JBV = [JBV1; JBV3];
% P2 = 0;             % 设定边界压强
% P4 = 1000;          % 设定边界压强
% JBP2 = [BE2 ,ones(size(BE2(:,1))) * P2, ones(size(BE2(:,1))) * P2];
% JBP4 = [BE4 ,ones(size(BE4(:,1))) * P4, ones(size(BE4(:,1))) * P4];
% JBP = [JBP2; JBP4];
%%%%%% 编写边界条件数据JBP与JBV

%%% Case2
%%% 方腔边界
JBV1 =[ BP1', u1 * ones(size(BP1))', v1 * ones(size(BP1))'];
JBV3 =[ BP3', u3 * ones(size(BP3))', v3 * ones(size(BP3))'];
u2 = 0; v2 = 0;     % 设定2边界速度
u4 = 0; v4 = 0;     % 设定4边界速度
JBV2 =[ BP2, u2 * ones(size(BP2)), v2 * ones(size(BP2))];
JBV4 =[ BP4, u2 * ones(size(BP4)), v4 * ones(size(BP4))];
JBV = [JBV1; JBV3; JBV2; JBV4];
P2 = 0;             % 设定方腔边界压强
P4 = 0;             % 设定方腔边界压强
JBP2 = [BE2 ,ones(size(BE2(:,1))) * P2, ones(size(BE2(:,1))) * P2];
JBP4 = [BE4 ,ones(size(BE4(:,1))) * P4, ones(size(BE4(:,1))) * P4];
JBP = [JBP2; JBP4];
%%% 方腔边界

clear JBV1 JBV3 BP1 BP3 BP4
clear JBP2 JBP4 P2 P4
clear BE1 BE2 BE3 BE4
clear u1 v1 u3 v3
%%%%%% 设定边界条件

%%%%%% 初始化总体方程的各个子块矩阵
B1 = zeros(Nd, Nz);
B2 = zeros(Nd, Nz);
% D11 = zeros(Nz, Nz);
% D12 = zeros(Nz, Nz);
% D21 = zeros(Nz, Nz);
% D22 = zeros(Nz, Nz);
C1 = zeros(Nz ,Nd);
C2 = zeros(Nz ,Nd);
F1 = zeros(Nz ,1);
F2 = zeros(Nz ,1);

% B1 = sparse(zeros(Nd, Nz));   %采用稀疏阵（稀疏阵赋值特别慢）
% B2 = sparse(zeros(Nd, Nz));
D11 = sparse(zeros(Nz, Nz));
D12 = sparse(zeros(Nz, Nz));
D21 = sparse(zeros(Nz, Nz));
D22 = sparse(zeros(Nz, Nz));
% C1 = sparse(zeros(Nz ,Nd));
% C2 = sparse(zeros(Nz ,Nd));
% F1 = sparse(zeros(Nz ,1));
% F2 = sparse(zeros(Nz ,1));
%%%%%% 初始化总体方程的各个子块矩阵

%%%%%% 计算各单元刚度(子块)矩阵并装配为整体刚度的子块矩阵
for i = 1:E   % 遍历所有单元
    for ie = 1:9
        JXYe(ie,:) = JXYV(JMV(i, ie), :);   % 第i个速度单元的结点坐标数据
    end
    [Be1, Be2] = function_of_Be(JXYe);   % 调用function_of_Be
    [De11, De12, De21, De22] = function_of_De(JXYe, viscosity);   % 调用function_of_De
    [Ce1, Ce2] = function_of_Ce(JXYe);   % 调用function_of_Ce

    %%%%%% 装配B1, B2
    for m = 1:4
        for n = 1:9
            B1(JMP(i, m), JMV(i, n)) = B1(JMP(i, m), JMV(i, n)) + Be1(m, n);
            B2(JMP(i, m), JMV(i, n)) = B2(JMP(i, m), JMV(i, n)) + Be2(m, n);            
        end
    end
    %%%%%% 装配B1, B2

    %%%%%% 装配D11, D12, D21, D22
    for m = 1:9
        for n = 1:9
            D11(JMV(i, m), JMV(i, n)) = D11(JMV(i, m), JMV(i, n)) + De11(m, n);
            D12(JMV(i, m), JMV(i, n)) = D12(JMV(i, m), JMV(i, n)) + De12(m, n);
            D21(JMV(i, m), JMV(i, n)) = D21(JMV(i, m), JMV(i, n)) + De21(m, n);
            D22(JMV(i, m), JMV(i, n)) = D22(JMV(i, m), JMV(i, n)) + De22(m, n);  
        end
    end
    %%%%%% 装配D11, D12, D21, D22

    %%%%%% 装配C1, C2
    for m = 1:9
        for n = 1:4
            C1(JMV(i, m), JMP(i, n)) = C1(JMV(i, m), JMP(i, n)) + Ce1(m, n);
            C2(JMV(i, m), JMP(i, n)) = C2(JMV(i, m), JMP(i, n)) + Ce2(m, n);            
        end
    end
    %%%%%% 装配C1, C2

end
%%%%%% 计算各单元刚度(子块)矩阵并装配为整体刚度的子块矩阵

clear JXYe;
%%%%%% 计算各单元荷载(子块)矩阵并装配为整体荷载的子块矩阵
for i = 1:length(JBP(:,1))   % 遍历所有压强边界单元
    for ie = 1:9
        JXYe(ie,:) = JXYV(JMV(JBP(i, 1), ie), :);   % 第i个压强边界单元的速度结点坐标数据
    end
    [Fe1, Fe2] = function_of_Fe2(JXYe, JBP(i, :));   % 调用My_function_of_Fe

    %%%%%% 装配F1, F2
    for m = 1:9
        F1(JMV(JBP(i, 1), m), 1) = F1(JMV(JBP(i, 1), m), 1) + Fe1(m, 1);
        F2(JMV(JBP(i, 1), m), 1) = F2(JMV(JBP(i, 1), m), 1) + Fe2(m, 1);
    end
    %%%%%% 装配F1, F2
end
%%%%%% 计算各单元荷载(子块)矩阵并装配为整体荷载的子块矩阵


%%%%%% 组合为整体矩阵方程
% K = [D11 D12 -C1
%     D21 D22 -C2
%     B1 B2 zeros(Nd, Nd)];
K = sparse([D11 D12 -C1
    D21 D22 -C2
    B1 B2 zeros(Nd, Nd)]);   % 直接装配为稀疏阵
B = [ -F1; -F2; zeros(Nd, 1)];
%%%%%% 组合为整体矩阵方程

%%%%%% (强制)速度边界条件实现
N_matrix = 2 * Nz + Nd;      % B矩阵行数
for i = 1:length(JBV(:,1))   % 遍历所有速度边界结点
    II = JBV(i, 1);
    u = JBV(i, 2);
    for J = 1:N_matrix       % 遍历B矩阵所有行
        B(J) = B(J) - K(J, II) * u;
    end
    K(II, :) = zeros(1, N_matrix);
    K(:, II) = zeros(N_matrix, 1);
    K(II, II) = 1;
    B(II) = u;
end

for i = 1:length(JBV(:,1))   % 遍历所有速度边界结点
    II = Nz + JBV(i, 1);
    v = JBV(i, 3);
    for J = 1:N_matrix       % 遍历B矩阵所有行
        B(J) = B(J) - K(J, II) * v;
    end
    K(II, :) = zeros(1, N_matrix);
    K(:, II) = zeros(N_matrix, 1);
    K(II, II) = 1;
    B(II) = v;
end
%%%%%% (强制)速度边界条件实现


%%%%%% way2:(强制)速度边界条件实现
% N_matrix = 2 * Nz + Nd;      % B矩阵行数
% for i = 1:length(JBV(:,1))   % 遍历所有速度边界结点
%     II = JBV(i, 1);
%     u = JBV(i, 2);
%     K(II, :) = 0;
%     K(II, II) = 1;
%     B(II) = u;
% end
% 
% for i = 1:length(JBV(:,1))   % 遍历所有速度边界结点
%     II = Nz + JBV(i, 1);
%     v = JBV(i, 3);
%     K(II, :) = 0;
%     K(II, II) = 1;
%     B(II) = v;
% end
%%%%%% way2:(强制)速度边界条件实现

%%% 方腔压力约束点
K(2 * Nz + 1, :) = 0;
K(2 * Nz + 1, 2 * Nz + 1) = 1;
B(2 * Nz + 1) = 0;
%%% 方腔压力约束点

%%%%%% 清理内存, 存储线性矩阵
clear D11 D12 D21 D22 C1 C2 B1 B2 
clear F1 F2 De11 De12 De21 De22 JXYe
clear Be1 Be2 Ce1 Ce2 JBP JMP JXYP
clear Fe1 Fe2 P_element P_side P_value 
clear i ie u v
clear AAA1 AAA2 BBB1 BBB2 CCC1 CCC2 DDD1 DDD2
clear II J JBV JBV2 JBV4 m n thetal N_matrix
save matrix
%%%%%% 清理内存, 存储线性矩阵


