%%=========================================================
%  经典功率谱估计（直接法 vs 自相关法）
%  ——不同 SNR 下检测频率 f1、f2（修正版）
%=========================================================
clear; clc; close all;

%%-------------------- 参数设置 --------------------%%
N = 1024;
n = 0:N-1;
f1 = 0.1; f2 = 0.3;
A1 = 10; A2 = 5;
numMonte = 100;                 % 每个 SNR 下的平均次数
SNR_dB = -15:1:10;               % SNR范围(dB)

%%-------------------- 结果存储 --------------------%%
f1_est_direct = zeros(size(SNR_dB));
f2_est_direct = zeros(size(SNR_dB));
f1_est_corr   = zeros(size(SNR_dB));
f2_est_corr   = zeros(size(SNR_dB));

%%-------------------- 主循环 --------------------%%
for si = 1:length(SNR_dB)
    snr = SNR_dB(si);
    est_f1_d = zeros(1,numMonte);
    est_f2_d = zeros(1,numMonte);
    est_f1_c = zeros(1,numMonte);
    est_f2_c = zeros(1,numMonte);

    for mc = 1:numMonte
        %%--- 信号生成 ---%%
        x_clean = A1*sin(2*pi*f1*n + pi/3) + A2*sin(2*pi*f2*n + pi/4);
        Px = mean(abs(x_clean).^2);
        Pn = Px / (10^(snr/10));
        noise = randn(1,N);
        noise = sqrt(Pn) * noise / std(noise);
        x_noisy = x_clean + noise;

        %%--- (1) 直接法 ---%%
        Xk = fft(x_noisy, N);
        P_direct = abs(Xk/N).^2;
        % 只保留前半段 (0~0.5)
        P_direct = P_direct(1:N/2);
        freq_axis = (0:N/2-1)/N;

        [pks, locs] = findpeaks(P_direct, freq_axis, ...
            'SortStr','descend', 'NPeaks', 2);
        est_freqs_d = sort(locs);
        if numel(est_freqs_d) < 2
            est_freqs_d = [NaN NaN];
        end
        est_f1_d(mc) = est_freqs_d(1);
        est_f2_d(mc) = est_freqs_d(2);

        %%--- (2) 自相关法 ---%%
        Rxx = xcorr(x_noisy, 'biased');
        P_corr = abs(fft(Rxx, N));
        P_corr = P_corr(1:N/2);
        [pks2, locs2] = findpeaks(P_corr, freq_axis, ...
            'SortStr','descend', 'NPeaks', 2);
        est_freqs_c = sort(locs2);
        if numel(est_freqs_c) < 2
            est_freqs_c = [NaN NaN];
        end
        est_f1_c(mc) = est_freqs_c(1);
        est_f2_c(mc) = est_freqs_c(2);
    end

    % 平均去噪
    f1_est_direct(si) = mean(est_f1_d,'omitnan');
    f2_est_direct(si) = mean(est_f2_d,'omitnan');
    f1_est_corr(si)   = mean(est_f1_c,'omitnan');
    f2_est_corr(si)   = mean(est_f2_c,'omitnan');
end

%%-------------------- 输出结果表 --------------------%%
T = table(SNR_dB', f1_est_direct', f2_est_direct', f1_est_corr', f2_est_corr', ...
    'VariableNames', {'SNR_dB','f1_direct','f2_direct','f1_corr','f2_corr'});
disp('不同SNR下估计的频率结果（平均值）：');
disp(T);

%%-------------------- 绘图 --------------------%%
figure;
subplot(2,1,1);
plot(SNR_dB, f1_est_direct, 'bo-', 'LineWidth',1.4); hold on;
plot(SNR_dB, f1_est_corr, 'rs--', 'LineWidth',1.4);
yline(f1, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)');
ylabel('估计频率 f_1');
title('f_1 估计结果 vs SNR');
legend('直接法','自相关法','真实值','Location','best');

subplot(2,1,2);
plot(SNR_dB, f2_est_direct, 'bo-', 'LineWidth',1.4); hold on;
plot(SNR_dB, f2_est_corr, 'rs--', 'LineWidth',1.4);
yline(f2, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)');
ylabel('估计频率 f_2');
title('f_2 估计结果 vs SNR');
legend('直接法','自相关法','真实值','Location','best');
