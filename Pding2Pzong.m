function P9 = Pding2Pzong(p4, JMV)
for i = 1: length((JMV(:, 1)))   % 遍历所有速度单元
    %%%%%%% 构建单元压力
    pp = [p4(JMV(i, 1)), p4(JMV(i, 3)), p4(JMV(i, 9)), p4(JMV(i, 7))];
    %%%%%%% 构建单元压力

    %%%%%%% 结点1压力
    P9(JMV(i, 1)) = p4(JMV(i, 1));
    %%%%%%% 结点1压力
    %%%%%%% 结点2压力
    kesi = 0; ita = -1;
    Fyp = [ 1/4 * (1 - kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 + ita)
        1/4 * (1 - kesi) * (1 + ita)];
    P9(JMV(i, 2)) = pp * Fyp;
    %%%%%%% 结点2压力
    %%%%%%% 结点3压力
    P9(JMV(i, 3)) = p4(JMV(i, 3));
    %%%%%%% 结点3压力
    %%%%%%% 结点4压力
    kesi = -1; ita = 0;
    Fyp = [ 1/4 * (1 - kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 + ita)
        1/4 * (1 - kesi) * (1 + ita)];
    P9(JMV(i, 4)) = pp * Fyp;
    %%%%%%% 结点4压力
    %%%%%%% 结点5压力
    kesi = 0; ita = 0;
    Fyp = [ 1/4 * (1 - kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 + ita)
        1/4 * (1 - kesi) * (1 + ita)];
    P9(JMV(i, 5)) = pp * Fyp;
    %%%%%%% 结点5压力
    %%%%%%% 结点6压力
    kesi = 1; ita = 0;
    Fyp = [ 1/4 * (1 - kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 + ita)
        1/4 * (1 - kesi) * (1 + ita)];
    P9(JMV(i, 6)) = pp * Fyp;
    %%%%%%% 结点6压力
    %%%%%%% 结点7压力
    P9(JMV(i, 7)) = p4(JMV(i, 7));
    %%%%%%% 结点7压力
    %%%%%%% 结点8压力
    kesi = 0; ita = 1;
    Fyp = [ 1/4 * (1 - kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 - ita)
        1/4 * (1 + kesi) * (1 + ita)
        1/4 * (1 - kesi) * (1 + ita)];
    P9(JMV(i, 8)) = pp * Fyp;
    %%%%%%% 结点8压力
    %%%%%%% 结点9压力
    P9(JMV(i, 9)) = p4(JMV(i, 9));
    %%%%%%% 结点9压力
end