function plot_mjo_hovmoller(data_fast, data_slow, t_fast, t_slow, ...
    Lon_Label, I6, I7, significance_level)
%PLOT_MJO_HOVMOLLER  Plot fast/slow MJO composite Hovmöller diagrams.
%
%   plot_mjo_hovmoller(data_fast, data_slow, t_fast, t_slow, ...
%       Lon_Label, I6, I7, significance_level)
%
%   INPUT
%     data_fast  : [nTime x nLonSeg]  — Fast_MJO_Composite(I6:I7,:)' 的结果
%     data_slow  : [nTime x nLonSeg]  — Slow_MJO_Composite(I6:I7,:)' 的结果
%     t_fast     : [nTime x nLonSeg]  — t_fast(I6:I7,:)' 的结果
%     t_slow     : [nTime x nLonSeg]  — t_slow(I6:I7,:)' 的结果
%     Lon_Label  : cell
%     I6, I7     : 整数，经度索引 20°E–140°W
%     significance_level : t 阈值，例如 2.75 表示 |t| >= 2.75 画显著性标记
%
%   说明
%     - data_fast / data_slow 要已经完成：
%           data_fast = Fast_MJO_Composite(I6:I7, :)';
%           data_slow = Slow_MJO_Composite(I6:I7, :)';
%       即维度为 [time x lonSeg]，time 对应 [-30:+30]。
%     - t_fast / t_slow 同样需要转成 [time x lonSeg]。
%     - Lon_Label 与 I6:I7 一起用来生成 "20°E, 40°E, ..., 140°W" 这样的标签；
%

if nargin < 8 || isempty(significance_level)
    significance_level = 2.75;
end

% -------------------- 基本检查 --------------------
if isempty(data_fast) || isempty(data_slow)
    error('data_fast / data_slow 不能为空。请先计算合成场。');
end

[nTime, nLonSeg] = size(data_fast);
if size(data_slow,1) ~= nTime || size(data_slow,2) ~= nLonSeg
    error('data_fast 与 data_slow 尺寸不一致。');
end

if size(t_fast,1) ~= nTime || size(t_fast,2) ~= nLonSeg || ...
   size(t_slow,1) ~= nTime || size(t_slow,2) ~= nLonSeg
    error('t_fast / t_slow 尺寸必须与 data_fast / data_slow 相同。');
end

% -------------------- 配置参数 --------------------
% 两幅子图的位置
position = [ ...
    0.1401  0.1799  0.3625  0.7428; ...
    0.5703  0.1811  0.3547  0.7439];

titles = ["b. Fast (24 cases)                                           5.4 m/s", ...
          "c. Slow (25 cases)                                        2.7 m/s" ];

% 自定义 colormap
colortable = [ ...
      0   23   69; ...
      0   86  159; ...
      0  113  184; ...
     23  161  207; ...
    145  220  231; ...
    251  251  251; ...
    251  251  251; ...
    255  188  170; ...
    255  106   70; ...
    237   16   24; ...
    194    6   12; ...
    111    0    0];
colortable = colortable ./ 255;

% -------------------- 准备数据容器 --------------------
% data_save(:,:,1) 对应 fast, data_save(:,:,2) 对应 slow
data_save = zeros(nTime, nLonSeg, 2);
t_save    = zeros(nTime, nLonSeg, 2);

data_save(:,:,1) = data_fast;
data_save(:,:,2) = data_slow;
t_save(:,:,1)    = t_fast;
t_save(:,:,2)    = t_slow;

% -------------------- 生成经度标签索引 --------------------
% 先生成 X_Tick，再映射到 Lon_Label(I6:I7)
X_Tick = 5:12:81;   % 和你原脚本保持一致
X_Tick_Label_Pre = cell(1, I7-I6+1);
for i = I6:I7
    X_Tick_Label_Pre{i-I6+1} = Lon_Label{i};
end
X_Tick_Label = cell(1, numel(X_Tick));
for i = 1:numel(X_Tick)
    X_Tick_Label{i} = X_Tick_Label_Pre{X_Tick(i)};
end

% -------------------- 开始画图 --------------------
figure;
set(gcf, 'Color', 'w');
colormap(colortable);

for ii = 1:2

    subplot(1,2,ii);

    z  = data_save(:,:,ii);
    z1 = z;          % 负值部分
    z2 = z;          % 正值部分
    z1(z1>0)  = 0;
    z2(z2<=0) = 0;

    % ---------------- 等值面填色 ----------------
    contourf(z, 'LineStyle', 'none');
    hold on;

    % 负值（虚线）、正值（实线）等值线
    contour(z1, 'LineStyle', '--', 'LineWidth', 1.5, ...
            'LevelList', -30:5:-5, 'LineColor', 'k');
    hold on;
    contour(z2, 'LineStyle', '-',  'LineWidth', 1.5, ...
            'LevelList', 5:5:30,  'LineColor', 'k');
    hold on;

    % ---------------- 显著性打点 (可选) ----------------
    [XX, YY] = meshgrid(1:nLonSeg, 1:nTime);
    FLAG = find(abs(t_save(:,:,ii)) >= significance_level);
    if ~isempty(FLAG)
        Xdot = XX(FLAG);
        Ydot = YY(FLAG);
        % 如果需要打点，在这里解除注释：
        % scatter(Xdot, Ydot, 5, 'o', ...
        %     'MarkerFaceColor','k', 'MarkerEdgeColor','k');
    end

    % ---------------- 坐标轴与 colorbar 设置 ----------------
    c = colorbar;
    caxis([-30 30]);
    set(c, 'YTick', -30:5:30);
    c.LineWidth = 1;
    c.FontName  = 'Arial';
    c.FontSize  = 12;

    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize  = 16;
    ax.FontName  = 'Arial';
    ax.Position  = position(ii,:);

    title(titles(ii), 'FontSize', 16);

    % Y 轴：Days (-30…30)
    ylabel('Days');
    ax.YTick      = [1 11 21 31 41 51 61];
    ax.YTickLabel = {'-30','-20','-10','0','10','20','30'};

    % X 轴：经度（使用 Lon_Label）
    ax.XTick      = X_Tick;
    ax.XTickLabel = X_Tick_Label;

end

end
