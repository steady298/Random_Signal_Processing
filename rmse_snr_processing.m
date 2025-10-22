%%===============================================================
%  随机信号分析：经典功率谱估计下的 RMSE–SNR 关系（SNR控制版）
%  功能：研究频率估计精度随信噪比变化的规律
%  方法：通过指定 SNR(dB)，动态调整噪声功率
%  指标：频率估计的 RMSE（均方根误差）
%===============================================================

clear; clc; close all;

%%-------------------- 参数设置 --------------------%%
N = 1024;                     % 信号长度
f_true = [0.1, 0.3];          % 真实频率
A = [2, 5];                   % 幅度
SNR_dB_list = -5:1:70;        % SNR 取值范围 (dB)
num_trials = 300;             % 每个 SNR 下重复实验次数（蒙特卡洛次数）

% 初始化结果存储
rmse_f1 = zeros(1, length(SNR_dB_list));
rmse_f2 = zeros(1, length(SNR_dB_list));

%%-------------------- 原始信号生成 --------------------%%
n = 1:N;
x_clean = A(1)*sin(2*pi*f_true(1)*n + pi/3) + ...
          A(2)*sin(2*pi*f_true(2)*n + pi/4);   % 原始信号
R_x = xcorr(x_clean);
P_signal = sum(abs(fft(R_x, N)));              % 信号功率

%%-------------------- 主循环：不同 SNR --------------------%%
fprintf("开始仿真：经典功率谱估计 RMSE–SNR 分析\n");
for i = 1:length(SNR_dB_list)
    
    SNR_dB = SNR_dB_list(i);
    P_noise = P_signal / (10^(SNR_dB/10));  % 根据SNR计算噪声功率
    noise_sigma = sqrt(P_noise);            % 噪声标准差
    
    err1_sum = 0; err2_sum = 0;             % 误差累加
    
    for t = 1:num_trials
        % 生成高斯白噪声（满足当前 SNR）
        noise = noise_sigma * randn(1, N);
        y_noisy = x_clean + noise;
        
        % 功率谱估计
        R_y = xcorr(y_noisy, 'biased');
        P_y = abs(fft(R_y, N));
        
        % 仅取正频率部分
        spectrum_half = P_y(1:N/2);
        
        % 峰值检测
        [~, locs] = findpeaks(spectrum_half, 'SortStr', 'descend');
        if numel(locs) < 2
            continue; % 若噪声太强导致无法检测到两个峰，跳过
        end
        
        % 频率估计值
        f1_est = locs(2)/N;
        f2_est = locs(1)/N;
        
        % 累加误差平方
        err1_sum = err1_sum + (f1_est - f_true(1))^2;
        err2_sum = err2_sum + (f2_est - f_true(2))^2;
    end
    
    % 求平均均方根误差（RMSE）
    rmse_f1(i) = sqrt(err1_sum / num_trials);
    rmse_f2(i) = sqrt(err2_sum / num_trials);
    
    fprintf("SNR = %2d dB 完成\n", SNR_dB);
end

fprintf("全部仿真完成 ✅\n");

%%-------------------- 绘制结果 --------------------%%
figure('Name','RMSE–SNR关系图','Color','w');
plot(SNR_dB_list, rmse_f1, '-ob', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName','RMSE(f_1)');
hold on;
plot(SNR_dB_list, rmse_f2, '-sr', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName','RMSE(f_2)');
hold off;

% 图表美化
grid on;
title('经典功率谱估计下频率估计 RMSE 随 SNR 变化', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Root Mean Square Error (RMSE)', 'FontSize', 12);
legend('Location','northeast');

