clc
clear
load result_of_n1

%%%%%% 出口速度提取
for i = 1: length(BP2)
    UB2(i, 1) = ux_k_1(BP2(i), 1);
end
%%%%%% 出口速度提取