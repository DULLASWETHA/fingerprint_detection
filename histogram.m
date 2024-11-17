clear all;
clc;

img = imread('clear silicone (1).bmp');

if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

img_gray = double(img_gray);
[rows, cols] = size(img_gray);
result_lbp = zeros(rows, cols);
result_sym_diff = zeros(rows, cols);
mask_size = 3;
half_mask = floor(mask_size / 2);

uniform_patterns = [0; 255];
for i = 1:8
    pattern = zeros(1, 8);
    pattern(1:i) = 1;
    uniform_patterns = [uniform_patterns; bi2de(pattern)];
    
    pattern = ones(1, 8);
    pattern(1:i) = 0;
    uniform_patterns = [uniform_patterns; bi2de(pattern)];
end

for i = 1:6
    for j = i+1:7
        pattern = zeros(1, 8);
        pattern(1:i) = 0;
        pattern(i+1:j) = 1;
        pattern(j+1:end) = 0;
        uniform_patterns = [uniform_patterns; bi2de(pattern)];
        
        pattern = ones(1, 8);
        pattern(1:i) = 1;
        pattern(i+1:j) = 0;
        pattern(j+1:end) = 1;
        uniform_patterns = [uniform_patterns; bi2de(pattern)];
    end
end

mask_positions = [128, 64, 32, 16, 8, 4, 2, 1];

for r = 1 + half_mask : rows - half_mask
    for c = 1 + half_mask : cols - half_mask
        neighborhood = img_gray(r-half_mask:r+half_mask, c-half_mask:c+half_mask);
        center_pixel = img_gray(r, c);
        
        if center_pixel == 0
            result_lbp(r, c) = 58;
            continue;
        end
        
        abs_difference = abs(neighborhood - center_pixel);
        normalized_difference = abs_difference / center_pixel;

        binary_values = normalized_difference >= 0.1;
        binary_code = [binary_values(1) binary_values(2) binary_values(3) ...
                       binary_values(6) binary_values(9) binary_values(8) ...
                       binary_values(7) binary_values(4)];
        
        decimal_value = binary_code * mask_positions';
        
        idx = find(uniform_patterns == decimal_value);
        if ~isempty(idx) && numel(idx) == 1
            result_lbp(r, c) = idx - 1;
        else
            result_lbp(r, c) = 58;
        end
    end
end

for r = 1 + half_mask : rows - half_mask
    for c = 1 + half_mask : cols - half_mask
        neighborhood = img_gray(r-half_mask:r+half_mask, c-half_mask:c+half_mask);
        
        top_left = neighborhood(1, 1);
        top_center = neighborhood(1, 2);
        top_right = neighborhood(1, 3);
        right_center = neighborhood(2, 3);
        bottom_right = neighborhood(3, 3);
        bottom_center = neighborhood(3, 2);
        bottom_left = neighborhood(3, 1);
        left_center = neighborhood(2, 1);
        
        diffs = [top_left - bottom_right, top_center - bottom_center, ...
                 top_right - bottom_left, right_center - left_center];
        binary_values = diffs >= 0;
        binary_code = binary_values;
        
        decimal_value = binary_code * [8; 4; 2; 1];
        result_sym_diff(r, c) = decimal_value;
    end
end
pairs = [];

for r = 1:rows
    for c = 1:cols
        if result_lbp(r, c) ~= 58
            pairs = [pairs; result_lbp(r, c), result_sym_diff(r, c)];
        end
    end
end

x_values = pairs(:, 1);
y_values = pairs(:, 2);

edges_x = 0:58;
edges_y = 0:max(y_values);

[~, x_bins] = histc(x_values, edges_x);
[~, y_bins] = histc(y_values, edges_y);

x_bins(x_bins == 0) = 1;
y_bins(y_bins == 0) = 1;

hist2d = accumarray([x_bins, y_bins], 1, [numel(edges_x), numel(edges_y)]);

figure;
imagesc(edges_x, edges_y, hist2d');
colorbar;
xlabel('Uniform Pattern Index');
ylabel('Center-Symmetric Value');
title('2D Histogram of Uniform Pattern Index vs. Center-Symmetric Value');
axis xy;