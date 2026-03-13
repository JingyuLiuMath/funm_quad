function colors = get_colors(m)
%GET_COLORS 生成 m 种有区分度的颜色
%   colors = get_colors(m) 返回一个 m×3 的 RGB 颜色矩阵
%   每行代表一种颜色，值域为 [0, 1]
%
%   示例：
%       colors = get_colors(5);
%       plot(rand(10, 5), 'LineWidth', 2);  % 使用 5 种颜色绘制线条
%       colormap(colors);
%
%   输入：
%       m - 正整数，需要的颜色数量
%
%   输出：
%       colors - m×3 的 RGB 颜色矩阵

    if ~isscalar(m) || ~isnumeric(m) || m < 1 || mod(m, 1) ~= 0
        error('输入必须是正整数');
    end
    
    % 优先使用 MATLAB 默认的线条颜色（前 7 种区分度很好）
    default_colors = [
        0    0.4470    0.7410   % 蓝色
        0.8500    0.3250    0.0980   % 橙色
        0.9290    0.6940    0.1250   % 黄色
        0.4940    0.1840    0.5560   % 紫色
        0.4660    0.6740    0.1880   % 绿色
        0.3010    0.7450    0.9330   % 浅蓝
        0.6350    0.0780    0.1840   % 红色
    ];
    
    if m <= size(default_colors, 1)
        % 如果 m 小于等于默认颜色数，直接使用
        colors = default_colors(1:m, :);
    else
        % 需要更多颜色时，使用 HSV 色彩空间均匀分布
        colors = hsv(m);
    end
end
