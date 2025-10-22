%%=========================================================
%  随机信号分析：经典功率谱估计（直接法 vs 自相关法）
%=========================================================
clear; clc; close all;

%%-------------------- 参数设置 --------------------%%
N = 1024;                      % 信号长度
n = 0:N-1;                     % 时间索引
f1 = 0.1; f2 = 0.3;            % 两个信号频率
A1 = 10; A2 = 5;               % 信号幅度

%%-------------------- 信号生成 --------------------%%
x_clean = A1*sin(2*pi*f1*n + pi/3) + A2*sin(2*pi*f2*n + pi/4); % 原信号
noise = randn(1, N);           % 加性高斯白噪声
x_noisy = x_clean + noise;     % 含噪信号

%%-------------------- (1) 直接法功率谱估计 --------------------%%
Xk = fft(x_noisy, N);          % DFT 计算
freq_axis = (0:N-1)/N;         % 归一化频率轴（0~1 对应 0~2π）
P_direct = abs(Xk/N).^2;       % 功率谱密度（幅度平方归一化）

subplot(2,1,1);
plot(freq_axis, 10*log10(P_direct), 'b', 'LineWidth', 1.2);
title('直接法功率谱估计（DFT法）', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('归一化频率 (×π rad/sample)');
ylabel('功率谱密度 (dB)');
grid on;
xlim([0 1]);
%set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

%%-------------------- (2) 自相关函数法功率谱估计 --------------------%%
Rxx = xcorr(x_noisy, 'biased');    % 估计自相关函数
P_corr = abs(fft(Rxx, N));         % 自相关函数的 FFT 即功率谱估计
P_corr = P_corr / max(P_corr);     % 归一化

subplot(2,1,2);
plot(freq_axis, 10*log10(P_corr), 'r', 'LineWidth', 1.2);
title('自相关函数法功率谱估计', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('归一化频率 (×π rad/sample)');
ylabel('功率谱密度 (dB)');
grid on;
xlim([0 1]);
%set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);


figure;
plot(freq_axis, 10*log10(P_direct/max(P_direct)), 'b', 'LineWidth', 1.5); hold on;
plot(freq_axis, 10*log10(P_corr/max(P_corr)), 'r--', 'LineWidth', 1.5);
legend('直接法', '自相关法');
xlabel('归一化频率 (×π rad/sample)');
ylabel('功率谱密度 (dB)');
title('经典功率谱估计方法对比');
grid on;


%%-------------------- 效果对比说明 --------------------%%
% 两种方法结果比较说明：
% 1. 直接法计算简单，但存在频谱泄漏（因有限长度截断）。
% 2. 自相关法通过平滑功率谱减少噪声波动，但分辨率略低。
