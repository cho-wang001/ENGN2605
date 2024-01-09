clc;
clear;

L = double(imread('Lenna.png'));
T = double(imread('traffic.jpg'));

ak5 = (1/25).*ones(5, 5);

gk35 = gauss_kernel(3, 0.5);
gk55 = gauss_kernel(5, 0.5);
gk75 = gauss_kernel(7, 0.5);
gk510 = gauss_kernel(5, 1.0);
gk515 = gauss_kernel(5, 1.5);


fs35 = fspecial('gaussian', 3, 0.5);
fs55 = fspecial('gaussian', 5, 0.5);
fs75 = fspecial('gaussian', 7, 0.5);

result = isequal(gk35, fs35);

if result
    disp('The kernel is equal. (Window Size of 3, σ = 0.5)');
else
    disp('The kernel is not equal. (Window Size of 3, σ = 0.5)');
end

result = isequal(gk55, fs55);

if result
    disp('The kernel is equal. (Window Size of 5, σ = 0.5)');
else
    disp('The kernel is not equal. (Window Size of 5, σ = 0.5)');
end

result = isequal(gk75, fs75);

if result
    disp('The kernel is equal. (Window Size of 7, σ = 0.5)');
else
    disp('The kernel is not equal. (Window Size of 7, σ = 0.5)');
end

imlgk35 = imfilter(L, gk35);
imlgk55 = imfilter(L, gk55);
imlgk75 = imfilter(L, gk75);
imlgk510 = imfilter(L, gk510);
imlgk515 = imfilter(L, gk515);
imlak5 = imfilter(L, ak5);


imtgk35 = imfilter(T, gk35);
imtgk55 = imfilter(T, gk55);
imtgk75 = imfilter(T, gk75);
imtgk510 = imfilter(T, gk510);
imtgk515 = imfilter(T, gk515);
imtak5 = imfilter(T, ak5);

figure(1);
subplot(1, 2, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 2, 2);
imshow(imlgk35/256);
title('A very sharp image after gaussian filter. (Window Size of 3, σ = 0.5)');

figure(2);
subplot(1, 2, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 2, 2);
imshow(imlgk55/256);
title('A very sharp image after gaussian filter. (Window Size of 5, σ = 0.5)');

figure(3);
subplot(1, 2, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 2, 2);
imshow(imlgk75/256);
title('A very sharp image after gaussian filter. (Window Size of 7, σ = 0.5)');

figure(4);
subplot(1, 2, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 2, 2);
imshow(imtgk35/256);
title('An image degraded by noise after gaussian filter. (Window Size of 3, σ = 0.5)');

figure(5);
subplot(1, 2, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 2, 2);
imshow(imtgk55/256);
title('An image degraded by noise after gaussian filter. (Window Size of 5, σ = 0.5)');

figure(6);
subplot(1, 2, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 2, 2);
imshow(imtgk75/256);
title('An image degraded by noise after gaussian filter. (Window Size of 7, σ = 0.5)');

figure(7);
subplot(1, 3, 1);
imshow(imlgk35/256);
title('A very sharp image after gaussian filter. (Window Size of 3, σ = 0.5)');
subplot(1, 3, 2);
imshow(imlgk55/256);
title('A very sharp image after gaussian filter. (Window Size of 5, σ = 0.5)');
subplot(1, 3, 3);
imshow(imlgk75/256);
title('A very sharp image after gaussian filter. (Window Size of 7, σ = 0.5)');

figure(8);
subplot(1, 3, 1);
imshow(imtgk35/256);
title('An image degraded by noise after gaussian filter. (Window Size of 3, σ = 0.5)');
subplot(1, 3, 2);
imshow(imtgk55/256);
title('An image degraded by noise after gaussian filter. (Window Size of 5, σ = 0.5)');
subplot(1, 3, 3);
imshow(imtgk75/256);
title('An image degraded by noise after gaussian filter. (Window Size of 7, σ = 0.5)');

figure(9);
subplot(1, 3, 1);
imshow(imlgk55/256);
title('A very sharp image after gaussian filter. (Window Size of 5, σ = 0.5)');
subplot(1, 3, 2);
imshow(imlgk510/256);
title('A very sharp image after gaussian filter. (Window Size of 5, σ = 1.0)');
subplot(1, 3, 3);
imshow(imlgk515/256);
title('A very sharp image after gaussian filter. (Window Size of 5, σ = 1.5)');

figure(10);
subplot(1, 3, 1);
imshow(imtgk55/256);
title('An image degraded by noise after gaussian filter. (Window Size of 5, σ = 0.5)');
subplot(1, 3, 2);
imshow(imtgk510/256);
title('An image degraded by noise after gaussian filter. (Window Size of 5, σ = 1.0)');
subplot(1, 3, 3);
imshow(imtgk515/256);
title('An image degraded by noise after gaussian filter. (Window Size of 5, σ = 1.5)');

figure(11);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(imlgk515/256);
title('A very sharp image after gaussian filter. (Window Size of 5, σ = 1.5)');
subplot(1, 3, 3);
imshow(imlak5/256);
title('A very sharp image after 5 × 5 box Filter filter. (Window Size of 5)');

figure(12);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(imtgk515/256);
title('An image degraded by noise after gaussian filter. (Window Size of 5, σ = 1.5)');
subplot(1, 3, 3);
imshow(imtak5/256);
title('An image degraded by noise after 5 × 5 box Filter filter. (Window Size of 5)');

function [filter]=gauss_kernel(filter_size,sigma)

    filter = zeros(filter_size);
    
    [x, y] = meshgrid(-(filter_size-1)/2:(filter_size-1)/2, -(filter_size-1)/2:(filter_size-1)/2);

    for i = 1:filter_size    
        for j = 1:filter_size    
            filter(i,j) = exp(-((x(i,j))^2 + (y(i,j))^2)/(2*sigma^2));
        end
    end

    filter = filter / sum(filter(:));

end


  





