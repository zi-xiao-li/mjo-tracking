function data_all = track_mjo_all_years(Seg, K_real, ref_lon)
%TRACK_MJO_ALL_YEARS  Fit optimal propagation lines for all years.
%
%   data_all = TRACK_MJO_ALL_YEARS(Seg, K_real, ref_lon)
%
%   Seg     : struct from build_mjo_segments，包含字段
%             Seg.OLR   [nLonSeg x nDay x nYear]  — Hovmöller OLR
%             Seg.Time  [nDay x 3 x nYear]        — [year month day]
%             Seg.t0_year{iyr} [nEvent_yr x 3]    — 每年所有 t0 日期
%             Seg.Lon   [nLonSeg x 1]             — 分段经度（度）
%             Seg.std   [nLonSeg x 1]             — 每个经度的 std
%             Seg.mean  [nLonSeg x 1]             — 每个经度的 mean
%
%   K_real  : 试探相速 (m/s) 向量，例如 1:0.1:25
%   ref_lon : t0 所在的参考经度 (deg)，例如 90
%
%   OUTPUT:
%       data_all : [nEvent x 11] 矩阵
%          [StartY StartM StartD  EndY EndM EndD ...
%           Day0Y  Day0M  Day0D  Speed(m/s) StartLon EndLon  Bm]
%
%   注意：这里不再写 Excel，由调用者自行 writematrix/xlswrite。
%
%   Author: Xinyu Li
% ---------------------------------------------------------------------

%% 默认参数
if nargin < 2 || isempty(K_real)
    K_real = 1:0.1:25;
end
if nargin < 3 || isempty(ref_lon)
    ref_lon = 90;
end

%% 从 Seg 中取出数据
Segment_OLR  = Seg.OLR;      % lonSeg x nDay x nYear
Segment_Time = Seg.Time;     % nDay x 3 x nYear
Segment_t0   = Seg.t0_year;  % cell, 每年一个 cell
Segment_Lon  = Seg.Lon(:);   % lonSeg x 1
Segment_mean = Seg.mean(:);  % lonSeg x 1

nLonSeg = numel(Segment_Lon);
nYear   = size(Segment_OLR, 3);

% Segment_std 第二列保存“经度索引”
Segment_std          = zeros(nLonSeg, 2);
Segment_std(:,1)     = Seg.std(:);
Segment_std(:,2)     = (1:nLonSeg)';    % index_x

%% 试探斜率 K 的经验公式转换
% 经验公式：K = 111001 / (34560 * V + 1)
K = 111001 ./ (34560 .* K_real + 1);

data_print = [];

%% 对每一年做循环
for Z = 1:nYear

    %—————————— 初始化这一年的变量 ——————————%
    ex_time    = []; 
    ex_OLR     = []; 
    ex_t0      = []; 
    ex_lon     = []; 
    ex_y0      = []; 
    ex_x0      = [];
    EX_t0      = [];
    ex_ref_y0  = [];

    Propagate_Start     = [];
    Propagate_End       = [];
    Propagate_Day_0     = [];
    Propagate_Start_Lon = [];
    Propagate_End_Lon   = [];
    Propagate_Speed     = [];
    Propagate_Bm        = [];

    % 当年冬季对应的切片、时间、t0
    ex_OLR  = Segment_OLR(:,:,Z)';          % time x lon
    ex_time = Segment_Time(:,1:3,Z);       % [nDay x 3]
    ex_lon  = Segment_Lon;                 % lonSeg
    EX_t0   = Segment_t0{Z};               % 该年所有 t0（年月日）

    ex_x0 = ref_lon;                       % 参考经度（度）

    % 如果该年没有发生 MJO 事件
    if isempty(EX_t0)
        continue
    end

    %===============================================================
    %            对这一年中的每一个 MJO 事件进行 tracking
    %===============================================================
    for H = 1:size(EX_t0,1)

        ex_t0(1,:) = EX_t0(H,:);   %#ok<AGROW>

        % 找到 t0 在该年冬季 hov 图中的行号（时间索引）
        ex_ref_y0 = find( ex_time(:,1)==ex_t0(1,1) & ...
                          ex_time(:,2)==ex_t0(1,2) & ...
                          ex_time(:,3)==ex_t0(1,3) );

        % 周围 ±12 天的时间窗口（-0.5 --坐标系习惯）
        ex_y0 = (ex_ref_y0-12 : ex_ref_y0+12) - 0.5;

        %-----------------------------------------------------------
        % 建立 trial lines：从参考经度 ref_lon 出发，扫描所有 slope K
        %-----------------------------------------------------------
        [~, x0] = min(abs(ex_lon(:) - ex_x0));   % 参考经度对应的列索引
        x0      = x0 - 0.5;                      

        % 右侧交点 (x>=x0)
        coords_right = zeros(length(ex_y0), length(K), 2); 
        for ii = 1:length(ex_y0)
            y0 = ex_y0(ii);
            for jj = 1:length(K)
                k  = K(jj);
                y  = @(x) k*(x-x0) + y0;
                if y(81) <= 243      % 未超过右边界
                    coords_right(ii,jj,1) = 81;
                    coords_right(ii,jj,2) = y(81);
                else                 % 超过下边界
                    y2 = @(x) (1./k)*(x-y0) + x0;
                    coords_right(ii,jj,1) = y2(243);
                    coords_right(ii,jj,2) = 243;
                end
            end
        end

        % 左侧交点 (x<=x0)
        coords_left = zeros(length(ex_y0), length(K), 2); 
        for ii = 1:length(ex_y0)
            y0 = ex_y0(ii);
            for jj = 1:length(K)
                k  = K(jj);
                y  = @(x) k*(x-x0) + y0;
                if y(0) >= 0         % 未超过上边界
                    coords_left(ii,jj,1) = 0;
                    coords_left(ii,jj,2) = y(0);
                else                 % 超过上边界
                    y2 = @(x) (1./k)*(x-y0) + x0;
                    coords_left(ii,jj,1) = y2(0);
                    coords_left(ii,jj,2) = 0;
                end
            end
        end

        %-----------------------------------------------------------
        % 对每一条 trial line，求其经过的所有格点（segs_info）
        %-----------------------------------------------------------
        segs_info = cell(size(coords_left,1), size(coords_left,2));

        for ii = 1:size(coords_right,1)
            y0_local = ex_y0(ii);
            for jj = 1:size(coords_right,2)

                k   = K(jj);
                P1  = squeeze(coords_left(ii,jj,:))';
                P2  = squeeze(coords_right(ii,jj,:))';

                P1_x = P1(1);  P1_y = P1(2);
                P2_x = P2(1);  P2_y = P2(2);

                xmin = min(P1_x, P2_x);
                xmax = max(P1_x, P2_x);
                ymin = min(P1_y, P2_y);
                ymax = max(P1_y, P2_y);

                % 沿 x 方向取整
                yx = @(x) k*(x-x0) + y0_local;
                xdx = []; xdy = [];
                for i = ceil(xmin):floor(xmax)
                    xdx(end+1) = i;                  %#ok<AGROW>
                    yi        = yx(i);
                    if (yi < 0 && yi > -1)
                        yi = 0;
                    end
                    xdy(end+1) = yi;                 %#ok<AGROW>
                end

                % 沿 y 方向取整
                xy = @(x) (1./k)*(x-y0_local) + x0;
                ydx = []; ydy = [];
                for j = ceil(ymin):floor(ymax)
                    ydy(end+1) = j;                  %#ok<AGROW>
                    xi        = xy(j);
                    if (xi < 0 && xi > -1)
                        xi = 0;
                    end
                    ydx(end+1) = xi;                 %#ok<AGROW>
                end

                % 所有端点与交点
                SP = unique([P1_x,P2_x,xdx,ydx; P1_y,P2_y,xdy,ydy]', 'rows');

                segs = struct([]);
                for t = 1:size(SP,1)-1
                    segs(t).index_x = max(ceil(SP(t+1,1)), ceil(SP(t,1))); % 列（x）
                    segs(t).index_y = max(ceil(SP(t+1,2)), ceil(SP(t,2))); % 行（y）
                    segs(t).OLRA    = ex_OLR(segs(t).index_y, segs(t).index_x);
                end
                segs_info{ii,jj} = segs;
            end
        end

        %===========================================================
        % 1) 识别 OLRA < (mean - 1*STD) 的格点
        % 2) 合并经向跨度 <=10° 的碎片
        % 3) 找到最长的 segment，计算 A 和 L
        %===========================================================
        for ii = 1:size(segs_info,1)
            for jj = 1:size(segs_info,2)
                var1   = segs_info{ii,jj};
                Series = zeros(1, numel(var1));

                % 识别所有 < mean - 1SD 的格点
                for t = 1:length(var1)
                    xidx = var1(t).index_x;  % 经向格点索引
                    mean_x = Segment_mean(xidx);
                    std_x  = Segment_std(xidx,1);

                    if var1(t).OLRA <= mean_x - std_x
                        var1(t).track = 1;
                    else
                        var1(t).track = 0;
                    end
                    Series(t) = var1(t).track;
                end

                % 合并间隔不超过 10°（~4 格点）的 0 段
                temp    = diff(Series);
                aa_start= find(temp==-1) + 1; % 0 段开始
                aa_end  = find(temp==1);      % 0 段结束

                if Series(1)==0
                    aa_start = [1 aa_start];
                end
                if Series(end)==0
                    aa_end = [aa_end length(Series)];
                end

                for i = 1:length(aa_end)
                    if var1(aa_end(i)).index_x - var1(aa_start(i)).index_x <= 4
                        Series(aa_start(i):aa_end(i)) = 1;
                    end
                end

                for t = 1:length(var1)
                    var1(t).track = Series(t);
                end
                segs_info{ii,jj} = var1;

                % 找到最长的 1 序列（Longest Segment）
                temp    = diff(Series);
                aa_start= find(temp==1) + 1;  % 1 段开始
                aa_end  = find(temp==-1);     % 1 段结束

                if Series(1)==1
                    aa_start = [1 aa_start];
                end
                if Series(end)==1
                    aa_end = [aa_end length(Series)];
                end

                if isempty(aa_end)
                    segs_info{ii,jj}(1).StartLong = 0;
                    segs_info{ii,jj}(1).EndLong   = 0;
                    segs_info{ii,jj}(1).StartDate = 0;
                    segs_info{ii,jj}(1).EndDate   = 0;
                    segs_info{ii,jj}(1).A         = 0;
                    segs_info{ii,jj}(1).L         = 0;
                else
                    [~,b]    = max(aa_end-aa_start);
                    bb_start = aa_start(b);
                    bb_end   = aa_end(b);

                    cc = 0;
                    for t = bb_start:bb_end
                        cc = cc + segs_info{ii,jj}(t).OLRA;
                    end

                    segs_info{ii,jj}(1).StartLong = segs_info{ii,jj}(bb_start).index_x;
                    segs_info{ii,jj}(1).EndLong   = segs_info{ii,jj}(bb_end).index_x;
                    segs_info{ii,jj}(1).StartDate = segs_info{ii,jj}(bb_start).index_y;
                    segs_info{ii,jj}(1).EndDate   = segs_info{ii,jj}(bb_end).index_y;
                    segs_info{ii,jj}(1).A         = cc;
                    segs_info{ii,jj}(1).L         = 2.5 * ...
                        (segs_info{ii,jj}(bb_end).index_x - ...
                         segs_info{ii,jj}(bb_start).index_x);
                end
            end
        end

        %===========================================================
        % 计算 A(t,c)、L(t,c)、B(t,c) 并找出 B 最大的那条轨迹
        %===========================================================
        A = zeros(size(segs_info,1), size(segs_info,2));
        L = zeros(size(segs_info,1), size(segs_info,2));

        for ii = 1:size(segs_info,1)
            for jj = 1:size(segs_info,2)
                A(ii,jj) = abs(segs_info{ii,jj}(1).A);
                L(ii,jj) = abs(segs_info{ii,jj}(1).L);
            end
        end

        A_m = max(A(:));
        L_m = max(L(:));
        B   = A./A_m + L./L_m;

        B_m        = max(B(:));
        [pos_x,pos_y] = find(B == B_m);

        % 处理多解（取 pos_x / pos_y 最大值）
        if numel(pos_x) ~= 1
            pos_x = max(pos_x);
        end
        if numel(pos_y) ~= 1
            pos_y = max(pos_y);
        end

        Final_K    = K_real(pos_y);          % 对应的相速 (m/s)
        Final_info = segs_info{pos_x,pos_y}; % 对应线段详细信息

        % 起止日期 & 经度
        start_date = ex_time(Final_info(1).StartDate,:);
        end_date   = ex_time(Final_info(1).EndDate,:);
        day0_date  = ex_t0;

        % 经度
        start_lon  = Final_info(1).StartLong*2.5 + 20;
        end_lon    = Final_info(1).EndLong*2.5 + 20;

        % 也可选择更通用的：
        % start_lon = Segment_Lon(Final_info(1).StartLong);
        % end_lon   = Segment_Lon(Final_info(1).EndLong);

        fprintf('Start: %4d-%02d-%02d   End: %4d-%02d-%02d   Day0: %4d-%02d-%02d   V=%.2f m/s\n',...
            start_date(1),start_date(2),start_date(3), ...
            end_date(1),end_date(2),end_date(3), ...
            day0_date(1), day0_date(2), day0_date(3), ...
            Final_K);

        % 累积到这一年
        Propagate_Start     = [Propagate_Start;     start_date];
        Propagate_End       = [Propagate_End;       end_date];
        Propagate_Day_0     = [Propagate_Day_0;     day0_date];
        Propagate_Start_Lon = [Propagate_Start_Lon; start_lon];
        Propagate_End_Lon   = [Propagate_End_Lon;   end_lon];
        Propagate_Speed     = [Propagate_Speed;     Final_K];
        Propagate_Bm        = [Propagate_Bm;        B_m];

    end % H (event loop)

    % 把这一年的所有事件整理成矩阵并追加到 data_print
    if ~isempty(Propagate_Start)
        data_year = [Propagate_Start, ...
                     Propagate_End, ...
                     Propagate_Day_0, ...
                     Propagate_Speed, ...
                     Propagate_Start_Lon, ...
                     Propagate_End_Lon, ...
                     Propagate_Bm];
        data_print = [data_print; data_year];
    end

end % Z (year loop)

data_all = data_print;
end

