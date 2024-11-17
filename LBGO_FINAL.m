clear all;
clc;
img = imread('unsigned_live (1).bmp');
%scaleFactor = 0.5; % Scale down to 50%
%img = imresize(img, scaleFactor);
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end
img_gray = double(img_gray);
[rows, cols] = size(img_gray);
result = zeros(rows, cols);
mask_size = 3;
half_mask = floor(mask_size / 2);
for r = 1 + half_mask : rows - half_mask
    for c = 1 + half_mask : cols - half_mask
        neighborhood = img_gray(r-half_mask:r+half_mask, c-half_mask:c+half_mask);
        center_pixel = img_gray(r, c);
        top_left = neighborhood(1, 1);
        top_center = neighborhood(1, 2);
        top_right = neighborhood(1, 3);
        right_center = neighborhood(2, 3);
        bottom_right = neighborhood(3, 3);
        bottom_center = neighborhood(3, 2);
        bottom_left = neighborhood(3, 1);
        left_center = neighborhood(2, 1);
        diffs = [top_left - bottom_right, top_center - bottom_center, top_right - bottom_left, right_center - left_center];
        binary_values = diffs >= 0;
        binary_code = binary_values;
        decimal_value = binary_code * [8; 4; 2; 1];
        result(r, c) = decimal_value;
    end
end
disp('Decimal values based on center-symmetric pixel differences:');
disp(result);
figure;
imshow(result, [], 'DisplayRange', []);
title('Decimal Values from Center-Symmetric Pixel Differences');

% Assuming 'result' is your 59x16 matrix

% Define the number of bins
numBins = 20; % For values ranging from 0 to 255
histogramValues = zeros(1, numBins);

% Count the frequency of each value in the result matrix
for i = 1:size(result, 1)
    for j = 1:size(result, 2)
        value = result(i, j);
        if value >= 1 && value <= numBins % Ensure value is within expected range
            histogramValues(value + 1) = histogramValues(value + 1) + 1; % +1 for MATLAB indexing
        end
    end
end

% Create the histogram plot
figure;
bar(0:numBins-1, histogramValues);
title('Histogram of Result Values');
xlabel('Value');
ylabel('Frequency');