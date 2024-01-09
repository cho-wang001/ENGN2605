close all;
clc;
clear;

I1 = imread("bds1.jpg");
I2 = imread("bds2.jpg");
I3 = imread("bds3.jpg");
I4 = imread("bds4.jpg");
I5 = imread("bds5.jpg");
I6 = imread("bds6.jpg");



i1 = intensity_edge_detector(I1, 3);
i2 = intensity_edge_detector(I2, 4);
i3 = intensity_edge_detector(I3, 3);
i4 = intensity_edge_detector(I4, 8);
i5 = intensity_edge_detector(I5, 18.5);
i6 = intensity_edge_detector(I6, 3);

figure()
subplot(1,2,1);
imshow(I1);
subplot(1,2,2);
imshow(i1);
title('edge img (binary edge by thresholding)');

figure()
subplot(1,2,1);
imshow(I2);
subplot(1,2,2);
imshow(i2);
title('edge img (binary edge by thresholding)');

figure()
subplot(1,2,1);
imshow(I3);
subplot(1,2,2);
imshow(i3);
title('edge img (binary edge by thresholding)');

figure()
subplot(1,2,1);
imshow(I4);
subplot(1,2,2);
imshow(i4);
title('edge img (binary edge by thresholding)');

figure()
subplot(1,2,1);
imshow(I5);
subplot(1,2,2);
imshow(i5);
title('edge img (binary edge by thresholding)');

figure()
subplot(1,2,1);
imshow(I6);
subplot(1,2,2);
imshow(i6);
title('edge img (binary edge by thresholding)');



function edge_img = intensity_edge_detector(image,threshold)
    if size(image, 3) == 3
        image = double(rgb2gray(image));
    else
        image = double(image);
    end

    fs31 = fspecial('gaussian', 3, 1.0);

    image = imfilter(image, fs31, 'same');

    [dy,dx] = gradient(image);

    M = sqrt(dx.^2+dy.^2);

    figure();
    imshow(M);
    title('Magnitude Map M');
    
    figure();
    quiver(dy, dx);
    
    gdm = atan2(dy, dx);
    directions = [0, pi/4, pi/2, 3*pi/4, pi, -3*pi/4, -pi/2, -pi/4];

    [rows, cols] = size(M);

    for i = 2:rows-1
        for j = 2:cols-1

            [~, index] = min(abs(gdm(i, j) - directions));
        
            if index == 1 || index == 5
                q = M(i, j+1);
                q_bar = M(i, j-1);
            elseif index == 2 || index == 6
                q = M(i-1, j+1);
                q_bar = M(i+1, j-1);
            elseif index == 3 || index == 7
                q = M(i-1, j);
                q_bar = M(i+1, j);
            else
                q = M(i-1, j-1);
                q_bar = M(i+1, j+1);
            end
        
            if M(i, j) < max(q, q_bar)
                M(i, j) = 0;
            end
        end
    end
    
    figure();
    imshow(M);

    edge_img = M > threshold;

end

