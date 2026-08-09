function colors = get_colors(m)
%GET_COLORS 生成 m 种有区分度的深色
%   colors = get_colors(m) 返回一个 m×3 的 RGB 颜色矩阵
%   每行代表一种颜色，值域为 [0, 1]，整体偏深以保证白色背景上的可读性
%
%   示例：
%       colors = get_colors(5);
%       plot(rand(10, 5), 'LineWidth', 2);
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

    % 深色默认调色板：降低亮度、提高饱和度
    default_colors = [
        0.0000    0.2000    0.5000   % 深蓝
        0.6000    0.1500    0.0000   % 深橙/砖红
        0.5000    0.3500    0.0000   % 深金/橄榄
        0.3500    0.0000    0.4500   % 深紫
        0.0000    0.4500    0.1500   % 深绿
        0.0000    0.3500    0.4000   % 深青
        0.5500    0.0000    0.1000   % 深红
        0.4000    0.2500    0.1000   % 深褐
        0.2500    0.2500    0.2500   % 深灰
    ];

    if m <= size(default_colors, 1)
        colors = default_colors(1:m, :);
    else
        % 扩展颜色：HSV 生成后压暗、加饱和
        h = linspace(0, 1, m+1); h(end) = [];   % 均匀色相
        s = 0.85 * ones(m, 1);                   % 高饱和
        v = 0.60 * ones(m, 1);                   % 低亮度（深色）
        colors = hsv2rgb([h(:), s, v]);
    end
end