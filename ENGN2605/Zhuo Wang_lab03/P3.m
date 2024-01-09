close all;
clc;
clear;

I1 = imread("bds1.jpg");
I2 = imread("bds2.jpg");
I3 = imread("bds3.jpg");
I4 = imread("bds4.jpg");
I5 = imread("bds5.jpg");
I6 = imread("bds6.jpg");

I1 = double(rgb2gray(I1));
I2 = double(rgb2gray(I2));
I3 = double(rgb2gray(I3));
I4 = double(rgb2gray(I4));
I5 = double(rgb2gray(I5));
I6 = double(rgb2gray(I6));

ei1 = edge(I1, "canny");
ei2 = edge(I2, "canny");
ei3 = edge(I3, "canny");
ei4 = edge(I4, "canny");
ei5 = edge(I5, "canny");
ei6 = edge(I6, "canny");

figure();
imshow(ei1);
title("Canny Edge Detector bd1");

figure();
imshow(ei2);
title("Canny Edge Detector bd2");

figure();
imshow(ei3);
title("Canny Edge Detector bd3");

figure();
imshow(ei4);
title("Canny Edge Detector bd4");

figure();
imshow(ei5);
title("Canny Edge Detector bd5");

figure();
imshow(ei6);
title("Canny Edge Detector bd6");

low_threshold = 0.0005;
high_threshold = 0.08;

ei1dt = edge(I1, 'canny', [low_threshold, high_threshold]);

figure();
subplot(1,2,1);
imshow(ei1);
title("Canny Edge Detector bd1");
subplot(1,2,2);
imshow(ei1dt);
title("Canny Edge Detector bd1 with double thresholds");

gaussian_sigma = 2.5;

ei2gs = edge(I2, 'canny', [], gaussian_sigma);

figure();
subplot(1,2,1);
imshow(ei2);
title("Canny Edge Detector bd1");
subplot(1,2,2);
imshow(ei2gs);
title("Canny Edge Detector bd1 with different gaussian sigma");

