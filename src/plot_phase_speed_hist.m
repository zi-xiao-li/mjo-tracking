function [h, ax] = plot_phase_speed_hist(speed)
%PLOT_PHASE_SPEED_HIST  Plot histogram of MJO phase speeds.
%
%   [h, ax] = PLOT_PHASE_SPEED_HIST(speed)
%
%   INPUT
%     speed : 列向量或行向量，MJO 相速度 (m/s)
%
%   OUTPUT
%     h  : histogram 对象句柄
%     ax : 坐标轴句柄
%


figure;

h = histogram(speed);
h.FaceColor = [255 129 113]./255;  % 填充颜色
h.BinWidth  = 0.5;                 % 直方图宽度

ax = gca;
ax.LineWidth = 2;
ax.FontSize  = 16;
ax.FontName  = 'Arial';

% y 轴范围与刻度
yticks = 0:2:20;
ax.YTick      = yticks;
% 原脚本最后一个标签为空字符串
ax.YTickLabel = [arrayfun(@num2str, yticks(1:end-1), 'UniformOutput', false), {''}];
ax.YLim       = [0 19];

xlabel('MJO Phase Speed (m/s)', 'FontSize', 16, 'FontName', 'Arial');
ylabel('Number of Events',      'FontSize', 16, 'FontName', 'Arial');

title('a.', 'FontSize', 16);  % 如果你想添加更长的标题，可自行修改

% colorbar
c = colorbar;
c.LineWidth = 1;
c.FontName  = 'Arial';
c.FontSize  = 12;

end

