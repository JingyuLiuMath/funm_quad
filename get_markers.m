function markers = get_markers(m)
%GET_MARKERS 生成 m 种有区分度的标记样式
%   markers = get_markers(m) 返回一个 1×m 的 cell 数组，每个元素是一种标记
%
%   示例：
%       markers = get_markers(5);
%       for i = 1:5
%           plot(i, rand(), 'Marker', markers{i}, 'MarkerSize', 10);
%           hold on;
%       end
%
%   输入：
%       m - 正整数，需要的标记数量
%
%   输出：
%       markers - 1×m 的 cell 数组，每个元素是标记字符串

    if ~isscalar(m) || ~isnumeric(m) || m < 1 || mod(m, 1) ~= 0
        error('输入必须是正整数');
    end
    
    % 按区分度排序的标记列表
    % 优先选择形状差异大的标记
    marker_list = {
        'o'           % 圆圈
        's'           % 方块
        'd'           % 菱形
        '^'           % 上三角
        'v'           % 下三角
        '*'           % 星号
        'p'           % 五角星
        '+'           % 加号
        'x'           % 叉号
        '<'           % 左三角
        '>'           % 右三角
        'h'           % 六角星
        '.'           % 点
    };
    
    num_basic = length(marker_list);
    
    if m <= num_basic
        % 如果 m 在基本标记范围内，直接返回
        markers = marker_list(1:m);
    else
        % 需要更多标记时，循环使用但添加填充变体
        markers = cell(1, m);
        for i = 1:m
            idx = mod(i-1, num_basic) + 1;
            markers{i} = marker_list{idx};
        end
        
        % 注意：MATLAB 没有直接的 "filled" 字符串标记
        % 如果需要在同一个图中区分超过13种标记，建议结合颜色使用
        warning('请求了 %d 种标记，但 MATLAB 基本标记只有 %d 种。建议结合 get_colors 函数使用颜色+标记组合。', m, num_basic);
    end
end
