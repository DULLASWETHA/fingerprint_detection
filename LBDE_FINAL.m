clear all;
clc;
img = imread('clear silicone (1).bmp');
scaleFactor = 0.5;
img = imresize(img, scaleFactor);
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

img_gray = double(img_gray);
[rows, cols] = size(img_gray);
result = zeros(rows, cols);
mask_positions = [128, 64, 32, 16, 8, 4, 2, 1];
mask_size = 3;
half_mask = floor(mask_size / 2);
for r = 1 + half_mask : rows - half_mask
    for c = 1 + half_mask : cols - half_mask
        neighborhood = img_gray(r-half_mask:r+half_mask, c-half_mask:c+half_mask);
        center_pixel = img_gray(r, c);
        if center_pixel == 0
            result(r, c) = NaN;
            continue;
        end
        abs_difference = abs(neighborhood - center_pixel);
        normalized_difference = abs_difference / center_pixel;
        binary_values = normalized_difference >= 0.1;
        binary_code = [binary_values(1) binary_values(2) binary_values(3) binary_values(6) binary_values(9) binary_values(8) binary_values(7) binary_values(4)];
        decimal_value = binary_code * mask_positions';
        result(r, c) = decimal_value;
    end
end
disp('LBDE FIGURE');
disp(result);
figure;
imshow(result, [], 'DisplayRange', []);
title('LBDE FIGURE');  
numBins = 256;
histogramValues = zeros(1, numBins);
for i = 1:size(result, 1)
    for j = 1:size(result, 2)
        value = result(i, j);
        if value >= 1 && value <= numBins
            histogramValues(value + 1) = histogramValues(value + 1) + 1;
        end
    end
end
figure;
bar(0:numBins-1, histogramValues);
title('Histogram of Result Values');
xlabel('Value');
ylabel('Frequency');