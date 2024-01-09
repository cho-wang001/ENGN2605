clc;
clear;
close all;

%Run the 'vl_setup' file before running it
image1 = imread('Images\problem_3_and_4\set1\img1.png');
image2 = imread('Images\problem_3_and_4\set1\img2.png');

image3 = imread('Images\problem_3_and_4\set2\img1.png');
image4 = imread('Images\problem_3_and_4\set2\img2.png');

image5 = imread('Images\problem_3_and_4\set3\img1.png');
image6 = imread('Images\problem_3_and_4\set3\img2.png');

image7 = imread('Images\problem_3_and_4\set4\img1.png');
image8 = imread('Images\problem_3_and_4\set4\img2.png');

if size(image1, 3) == 3
    image1 = rgb2gray(image1);
end
if size(image2, 3) == 3
    image2 = rgb2gray(image2);
end
if size(image3, 3) == 3
    image3 = im2gray(image3);
end
if size(image4, 3) == 3
    image4 = im2gray(image4);
end
if size(image5, 3) == 3
    image5 = im2gray(image5);
end
if size(image6, 3) == 3
    image6 = im2gray(image6);
end
if size(image7, 3) == 3
    image7 = im2gray(image7);
end
if size(image8, 3) == 3
    image8 = im2gray(image8);
end

PeakThresh = 0.01;
PeakThresh1 = 10;

EdgeThresh = 10;
EdgeThresh1 = 100;

[f1, d1] = vl_sift(single(image1), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);
[f2, d2] = vl_sift(single(image2), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);

[f3, d3] = vl_sift(single(image3), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);
[f4, d4] = vl_sift(single(image4), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);

[f5, d5] = vl_sift(single(image5), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);
[f6, d6] = vl_sift(single(image6), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);

[f7, d7] = vl_sift(single(image7), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);
[f8, d8] = vl_sift(single(image8), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh);

%Play with the input parameters, e.g., PeakThresh and edgethresh, and plot SIFT keypoints (features) on four pairs
%of images in Figure 3.
[f9, d9] = vl_sift(single(image1), 'PeakThresh', PeakThresh, 'EdgeThresh', EdgeThresh1);
[f10, d10] = vl_sift(single(image1), 'PeakThresh', PeakThresh1, 'EdgeThresh', EdgeThresh);



x1 = f1(1,:);
y1 = f1(2,:);
x2 = f2(1,:);
y2 = f2(2,:);
x3 = f3(1,:);
y3 = f3(2,:);
x4 = f4(1,:);
y4 = f4(2,:);
x5 = f5(1,:);
y5 = f5(2,:);
x6 = f6(1,:);
y6 = f6(2,:);
x7 = f7(1,:);
y7 = f7(2,:);
x8 = f8(1,:);
y8 = f8(2,:);
x9 = f9(1,:);
y9 = f9(2,:);
x10 = f10(1,:);
y10 = f10(2,:);


figure();
subplot(2, 2, 1);
imshow(image1);
title('Image 1');

subplot(2, 2, 2);
imshow(image2);
title('Image 2');

subplot(2, 2, 3);
imshow(image1);
hold on;
scatter(x1, y1, 3, 'og');

subplot(2, 2, 4);
imshow(image2);
hold on;
scatter(x2, y2, 3, 'sb');


figure();
subplot(2, 2, 1);
imshow(image3);
title('Image 3');

subplot(2, 2, 2);
imshow(image4);
title('Image 4');

subplot(2, 2, 3);
imshow(image3);
hold on;
scatter(x3, y3, 3, 'og');

subplot(2, 2, 4);
imshow(image4);
hold on;
scatter(x4, y4, 3, 'sb');


figure();
subplot(2, 2, 1);
imshow(image5);
title('Image 5');

subplot(2, 2, 2);
imshow(image6);
title('Image 6');

subplot(2, 2, 3);
imshow(image5);
hold on;
scatter(x5, y5, 3, 'og');

subplot(2, 2, 4);
imshow(image6);
hold on;
scatter(x6, y6, 3, 'sb');


figure();
subplot(2, 2, 1);
imshow(image7);
title('Image 7');

subplot(2, 2, 2);
imshow(image8);
title('Image 8');

subplot(2, 2, 3);
imshow(image7);
hold on;
scatter(x7, y7, 3, 'og');

subplot(2, 2, 4);
imshow(image8);
hold on;
scatter(x8, y8, 3, 'sb');

figure();
subplot(1, 3, 1);
imshow(image1);
title('Image 1');

subplot(1, 3, 2);
imshow(image1);
hold on;
scatter(x9, y9, 3, 'or');

subplot(1, 3, 3);
imshow(image1);
hold on;
scatter(x10, y10, 3, 'or');
