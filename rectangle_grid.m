function rectangle_grid(JMP, JXYV)
E = length(JMP);  % 读取结点信息矩阵的行数（单元总数）
N = length(JXYV); % 读取节点坐标矩阵的行数（结点总数）
x = JXYV(: ,1);   % 建立所有结点x坐标向量
y = JXYV(:, 2);   % 建立所有结点y左边向量
ele_num_prnt = 1; % 是否显示单元号码
pnt_num_prnt = 1; % 是否显示结点号码
axis equal;
hold on;
del_x = 0.001; del_y = 0.001; % 单元号放置位置的调整变量
for k = 1:E
    for l = 1:4
        p = JMP(k, l);
        xx(l) = x(p); yy(l) = y(p); % 读取一个单元中的对应点的坐标
    end
    xx(5) = xx(1); yy(5) = yy(1);
    plot(xx, yy);
    x_cen = sum(xx(1:4))/4; y_cen = sum(yy(1:4))/4;
    if(ele_num_prnt == 1)  % 标记每个单元的单元号码
        text(x_cen-del_x, y_cen-del_y, int2str(k));
    end
end
if(pnt_num_prnt == 1) % 打印结点号
    for n = 1:N
        text(x(n), y(n), ['(',int2str(n),')']);
    end
end
axis off


