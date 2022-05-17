clc
clear
clf
%%%%%%% 区域几何尺寸及网格划分参数
H = 0.06; % 区域总高（m）
L = 0.08; % 区域总长（m）
Nx = 4; % 水平方向的网格数量
Ny = 4; % 竖直方向网格数量
theta = 0; % 网格平面旋转角度
%%%%%%% 区域几何尺寸及网格划分参数

%%%%%%% 总单元及总结点数
E = Nx * Ny; %总单元数
Nz = (2 * Nx + 1) * (2 * Ny + 1); % 二次单元结点总数
Nd = (Nx + 1) * (Ny + 1); % 线性单元结点总数
%%%%%%% 总单元及总结点数

%%%%%%% 单元间距
Dx = L/Nx/2; % 水平方向网格间距
Dy = H/Ny/2; % 数值方向网格间距
%%%%%%% 单元间距

%%%%%%% 节点分布拓扑
AAA = zeros(Ny * 2 + 1, Nx * 2 + 1);
for i = 1:2:2 * Nx + 1
    AAA(1, i) = (i + 1)/2;
end
for i = 1:Nx
    AAA(1, 2 * i) = (Nx + 1) * (Ny + 1) + i;
end
for i = 1:2 * Nx + 1
    AAA(2, i) = (Nx + 1) * (Ny + 1) + Nx + i;
end
for j = 3:2:2 * Ny + 1
    for i = 1:2:2 * Nx + 1
        AAA(j,i) = (i + 1)/2 + (Nx + 1) * (j - 1)/2;
    end
end
for j = 3:2:2 * Ny + 1
    for i = 1:Nx
        AAA(j, 2 * i) = (Nx + 1)*(Ny + 1) + (Nx + 2 * Nx + 1) * (j - 1)/2 + i;
    end
end
for j = 4:2:2 * Ny 
    for i = 1:2 * Nx + 1
        AAA(j, i) = (Nx + 1)*(Ny + 1) + (j/2) * Nx + (2 * Nx + 1) * (j/2 - 1) + i;
    end
end
%%%%%%% 节点分布拓扑

%%%%%%% JXYV(速度单元结点坐标数据)
%%%%%%% 四边形二次单元JXYV生成
for i = 1:2 * Ny + 1
    for j = 1:2 * Nx + 1
        JXYV(AAA(i,j), 1) = Dx * (j - 1);
        JXYV(AAA(i,j), 2) = Dy * (i - 1);
    end
end
%%%%%%% 四边形二次单元JXYV生成
%%%%%%% 网格平面旋转
for i = 1:length(JXYV(:,1))
    R = sqrt( (JXYV(i, 1) + 1)^2 + JXYV(i, 2)^2 );
    thetal = atan( JXYV(i,2)/(JXYV(i,1) + 1) );  % 节点位置向右平移1
    JXYV(i,1) = R * cos(theta/180 * pi + thetal);
    JXYV(i,2) = R * sin(theta/180 * pi + thetal);
end
%%%%%%% 网格平面旋转

%%%%%%% JMV_IEN(速度单元IEN)
%%%%%%% 四边形二次单元JMV_IEN生成
k = 0;
for i = 1:Ny
    for j = 1:Nx
        k = k + 1;
        JMV(k,:) = [AAA(2 * i - 1, 2 * j - 1), AAA(2 * i - 1, 2 * j),...
            AAA(2 * i - 1, 2 * j + 1), AAA(2 * i, 2 * j - 1), AAA(2 * i, 2 * j),...
            AAA(2 * i, 2 * j + 1), AAA(2 * i + 1, 2 * j -1), AAA(2 * i + 1, 2 * j),...
            AAA(2 * i + 1, 2 * j + 1),];
    end
end
%%%%%%% 四边形二次单元JMV_IEN生成

%%%%%%% JMP_IEN(压力单元IEN)，JXYP(压力单元结点坐标数据)
%%%%%%% 四边形线性单元JMP_IEN和JXYP生成
JMP = [JMV(:,1), JMV(:,3), JMV(:,9), JMV(:,7)];
JXYP = JXYV([1:Nd],:);
%%%%%%% 四边形线性单元JMP_IEN和JXYP生成

%%%%%%% BP(边界节点数据)
%%%%%%% BP数据生成
BP1 = AAA(1,:);
BP2 = AAA(:, 2 * Nx + 1);
BP3 = AAA(2 * Ny + 1, :);
BP4 = AAA(:,1);
%%%%%%% BP数据生成

%%%%%%% BE边界单元数据
%%%%%%% BE数据生成
thetax1 = pi/2 - theta/180 * pi;      % 1号边界外法线方向与x轴夹角
thetay1 = pi - theta/180 * pi;        % 1号边界外法线方向与y轴夹角
thetax2 = theta/180 * pi;             % 2号边界外法线方向与x轴夹角
thetay2 = pi/2 - thetax2;             % 2号边界外法线方向与y轴夹角
thetax3 = pi - pi/2 + theta/180 * pi; % 3号边界外法线方向与x轴夹角
thetay3 = pi - theta/180 * pi + pi;   % 3号边界外法线方向与y轴夹角
thetax4 = (180 + theta)/180 * pi;     % 4号边界外法线方向与x轴夹角
thetay4 = pi/2 + theta/180 * pi;      % 4号边界外法线方向与y轴夹角
AAA1 = ones(Nx, 1) * cos(thetax1);    % 底边方向余弦
AAA2 = ones(Nx, 1) * cos(thetay1);    % 底边方向余弦
BBB1 = ones(Ny, 1) * cos(thetax2);    % 右侧方向余弦
BBB2 = ones(Ny, 1) * cos(thetay2);    % 右侧方向余弦
CCC1 = ones(Nx, 1) * cos(thetax3);    % 上边方向余弦
CCC2 = ones(Nx, 1) * cos(thetay3);    % 上边方向余弦
DDD1 = ones(Ny, 1) * cos(thetax4);    % 左侧方向余弦
DDD2 = ones(Ny, 1) * cos(thetay4);    % 左侧方向余弦
BE1 = [ [1:Nx]', ones(size([1:Nx]')), AAA1, AAA2 ];
BE2 = [ [Nx:Nx:Ny * Nx]', 2 * ones(size([1:Ny]')), BBB1, BBB2 ];
BE3 = [ [Nx * (Ny - 1) + 1: Nx * Ny]', 3 * ones(size([1:Nx]')), CCC1, CCC2 ];
BE4 = [ [1:Nx:(Ny - 1) * Nx + 1]', 4 * ones(size([1:Ny]')), DDD1, DDD2 ];
%%%%%%% BE数据生成

%%%%%%% 调用四边形网格绘制程序
rectangle_grid(JMP,JXYV);
%%%%%%% 调用四边形网格绘制程序

%%%%%%% 清楚多余变量
clear Dx Dy H L Nx Ny i j k
clear theta theta1 R AAA
clear thetax1 thetax2 thetax3 thetax4
clear thetay1 thetay2 thetay3 thetay4
%%%%%%% 清楚多余变量

%%%%%%% 存储网格数据
save msh
save ('../source/matlab_repos/Computational mechanics/Project_Final/msh')
%%%%%%% 存储网格数据