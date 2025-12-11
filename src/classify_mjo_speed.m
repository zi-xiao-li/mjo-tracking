function [data_all, data_fast, data_slow, speed] = classify_mjo_speed( ...
    inFile, outFileAll, outFileFast, outFileSlow, ...
    speedSlow, speedFast, lonMin, lonMax)
%CLASSIFY_MJO_SPEED  Classify MJO events into slow/fast groups by phase speed.
%
%   [data_all, data_fast, data_slow, speed] = CLASSIFY_MJO_SPEED( ...
%       inFile, outFileAll, outFileFast, outFileSlow, ...
%       speedSlow, speedFast, lonMin, lonMax)
%
%   INPUT
%     inFile      : 输入的 Excel 文件
%     outFileAll  : 所有事件输出 Excel 路径
%     outFileFast : fast MJO 输出 Excel 路径
%     outFileSlow : slow MJO 输出 Excel 路径
%     speedSlow   : 慢速 MJO 阈值（例如 3.3 m/s，取 speed < speedSlow）
%     speedFast   : 快速 MJO 阈值（例如 4.5 m/s，取 speed > speedFast）
%     lonMin      : 起始经度筛选条件
%     lonMax      : 终止经度筛选条件
%
%   OUTPUT
%     data_all    : 经过经度筛选后的所有事件
%     data_fast   : 快速 MJO 子集（speed > speedFast）
%     data_slow   : 慢速 MJO 子集（speed < speedSlow）
%     speed       : 所有事件的相速度列向量（第 10 列）
%


%% 默认参数（如果未提供）
if nargin < 5 || isempty(speedSlow)
    speedSlow = 3.3;   % 原脚本默认
end
if nargin < 6 || isempty(speedFast)
    speedFast = 4.5;   % 原脚本默认
end
if nargin < 7 || isempty(lonMin)
    lonMin = 80;       % StartLon <= 80
end
if nargin < 8 || isempty(lonMax)
    lonMax = 120;      % EndLon   >= 120
end

%% 读取数据
fprintf('Reading MJO tracking data from: %s\n', inFile);
data_input = xlsread(inFile);   % 原始数据
data       = data_input;        % 保留一份原始副本

% 列意义
% 1-3: Start date [Y M D]
% 4-6: End   date [Y M D]
% 7-9: Day0  date [Y M D]
% 10 : Speed (m/s)
% 11 : StartLon (deg)
% 12 : EndLon   (deg)
% 13 : Bm

%% 经度过滤
idxLon = (data(:,11) <= lonMin) & (data(:,12) >= lonMax);
data   = data(idxLon,:);

%% 按速度分类
speed     = data(:,10);
data_slow = data(speed < speedSlow, :);   % 严格 <
data_fast = data(speed > speedFast, :);   % 严格 >
data_all  = data;

%% 写出三个 Excel
if ~isempty(outFileAll)
    fprintf('Writing ALL events to:   %s\n', outFileAll);
    xlswrite(outFileAll, data_all);
end

if ~isempty(outFileFast)
    fprintf('Writing FAST events to:  %s\n', outFileFast);
    xlswrite(outFileFast, data_fast);
end

if ~isempty(outFileSlow)
    fprintf('Writing SLOW events to:  %s\n', outFileSlow);
    xlswrite(outFileSlow, data_slow);
end

%% 打印统计信息
fprintf('---------------------------------------------\n');
fprintf('Average speed of SLOW MJO: %.2f m/s (N = %d)\n', ...
    mean(data_slow(:,10)), size(data_slow,1));
fprintf('Average speed of FAST MJO: %.2f m/s (N = %d)\n', ...
    mean(data_fast(:,10)), size(data_fast,1));
fprintf('Average speed of ALL  MJO: %.2f m/s (N = %d)\n', ...
    mean(data_all(:,10)),  size(data_all,1));
fprintf('---------------------------------------------\n');

end

