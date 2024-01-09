clc;
clear;

L = double(imread('Lenna.png'));
T = double(imread('traffic.jpg'));
A = double(imread('ImageA.jpg'));
B = double(imread('ImageB.jpg'));

f1 = (1/9).*ones(3,3);
f2 = (1/16).*[1, 2, 1; 2, 4, 2; 1, 2, 1];
f3 = (1/3).*[0, 0, 0; 1, 1, 1; 0, 0, 0];
f4 = (1/3).*[0, 1, 0; 0, 1, 0; 0, 1, 0];
f5 = (1/5).*[0, 1, 0; 1, 1, 1; 0, 1, 0];
f6 = (1/25).*ones(5,5);
f7 = (1/25).*[0, 1, 1, 1, 0; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 0, 1, 1, 1, 0];
f8 = [0, 1, 0; 1, -4, 1; 0, 1, 0];
f9 = [1, 1, 1; 1, -8, 1; 1, 1, 1];
f10 = [1, 0, -1; 1, 0, -1; 1, 0, -1];
f11 = [1, 1, 1; 0, 0, 0; -1, -1, -1];
f12 = [1, 0, -1; 2, 0, -2; 1, 0, -1];
f13 = [1, 2, 1; 0, 0, 0; -1, -2, -1];

of1l1 = imfilter(L, f1);
of1cl1 = my_conv(L, f1);
of1l = int8(imfilter(L, f1));
of1cl = int8(my_conv(L, f1));

of2l1 = imfilter(L, f2);
of2cl1 = my_conv(L, f2);
of2l = int8(imfilter(L, f2));
of2cl = int8(my_conv(L, f2));

of3l1 = imfilter(L, f3);
of3cl1 = my_conv(L, f3);
of3l = int8(imfilter(L, f3));
of3cl = int8(my_conv(L, f3));

of4l1 = imfilter(L, f4);
of4cl1 = my_conv(L, f4);
of4l = int8(imfilter(L, f4));
of4cl = int8(my_conv(L, f4));

of5l1 = imfilter(L, f5);
of5cl1 = my_conv(L, f5);
of5l = int8(imfilter(L, f5));
of5cl = int8(my_conv(L, f5));

of6l1 = imfilter(L, f6);
of6cl1 = my_conv(L, f6);
of6l = int8(imfilter(L, f6));
of6cl = int8(my_conv(L, f6));

of7l1 = imfilter(L, f7);
of7cl1 = my_conv(L, f7);
of7l = int8(imfilter(L, f7));
of7cl = int8(my_conv(L, f7));

result = isequal(of1l, of1cl);

if result
    disp('The results is equal. (Box Filter, a very sharp image)');
else
    disp('The results is not equal. (Box Filter, a very sharp image)');
end

result = isequal(of2l, of2cl);

if result
    disp('The results is equal. (Weighted Box Filter, a very sharp image)');
else
    disp('The results is not equal. (Weighted Box Filter, a very sharp image)');
end

result = isequal(of3l, of3cl);

if result
    disp('The results is equal. (Horizontal Filter, a very sharp image)');
else
    disp('The results is not equal. (Horizontal Filter, a very sharp image)');
end

result = isequal(of4l, of4cl);

if result
    disp('The results is equal. (Vertical Filter, a very sharp image)');
else
    disp('The results is not equal. (Vertical Filter, a very sharp image)');
end

result = isequal(of5l, of5cl);

if result
    disp('The results is equal. (Circular-Shape Filter, a very sharp image)');
else
    disp('The results is not equal. (Circular-Shape Filter, a very sharp image)');
end

result = isequal(of6l, of6cl);

if result
    disp('The results is equal. (5 × 5 Box Filter, a very sharp image)');
else
    disp('The results is not equal. (5 × 5 Box Filter, a very sharp image)');
end

result = isequal(of7l, of7cl);

if result
    disp('The results is equal. (5 × 5 Circular Shape Filter, a very sharp image)');
else
    disp('The results is not equal. (5 × 5 Circular Shape Filter, a very sharp image)');
end

of1t1 = imfilter(T, f1);
of1ct1 = my_conv(T, f1);
of1t = int8(imfilter(T, f1));
of1ct = int8(my_conv(T, f1));

of2t1 = imfilter(T, f2);
of2ct1 = my_conv(T, f2);
of2t = int8(imfilter(T, f2));
of2ct = int8(my_conv(T, f2));

of3t1 = imfilter(T, f3);
of3ct1 = my_conv(T, f3);
of3t = int8(imfilter(T, f3));
of3ct = int8(my_conv(T, f3));

of4t1 = imfilter(T, f4);
of4ct1 = my_conv(T, f4);
of4t = int8(imfilter(T, f4));
of4ct = int8(my_conv(T, f4));

of5t1 = imfilter(T, f5);
of5ct1 = my_conv(T, f5);
of5t = int8(imfilter(T, f5));
of5ct = int8(my_conv(T, f5));

of6t1 = imfilter(T, f6);
of6ct1 = my_conv(T, f6);
of6t = int8(imfilter(T, f6));
of6ct = int8(my_conv(T, f6));

of7t1 = imfilter(T, f7);
of7ct1 = my_conv(T, f7);
of7t = int8(imfilter(T, f7));
of7ct = int8(my_conv(T, f7));

result = isequal(of1t, of1ct);

if result
    disp('The results is equal. (Box Filter, an image degraded by noise)');
else
    disp('The results is not equal. (Box Filter, an image degraded by noise)');
end

result = isequal(of2t, of2ct);

if result
    disp('The results is equal. (Weighted Box Filter, an image degraded by noise)');
else
    disp('The results is not equal. (Weighted Box Filter, an image degraded by noise)');
end

result = isequal(of3t, of3ct);

if result
    disp('The results is equal. (Horizontal Filter, an image degraded by noise)');
else
    disp('The results is not equal. (Horizontal Filter, an image degraded by noise)');
end

result = isequal(of4t, of4ct);

if result
    disp('The results is equal. (Vertical Filter, an image degraded by noise)');
else
    disp('The results is not equal. (Vertical Filter, an image degraded by noise)');
end

result = isequal(of5t, of5ct);

if result
    disp('The results is equal. (Circular-Shape Filter, an image degraded by noise)');
else
    disp('The results is not equal. (Circular-Shape Filter, an image degraded by noise)');
end

result = isequal(of6t, of6ct);

if result
    disp('The results is equal. (5 × 5 Box Filter, an image degraded by noise)');
else
    disp('The results is not equal. (5 × 5 Box Filter, an image degraded by noise)');
end

result = isequal(of7t, of7ct);

if result
    disp('The results is equal. (5 × 5 Circular Shape Filter, an image degraded by noise)');
else
    disp('The results is not equal. (5 × 5 Circular Shape Filter, an image degraded by noise)');
end

figure(1);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of1cl1)/256);
title('A very sharp image after box filter by my conv');
subplot(1, 3, 3);
imshow(double(of1l1)/256);
title('A very sharp image after box filter by imfilter');

figure(2);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of2cl1)/256);
title('A very sharp image after weighted box filter by my conv');
subplot(1, 3, 3);
imshow(double(of2l1)/256);
title('A very sharp image after weighted box filter by imfilter');

figure(3);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of3cl1)/256);
title('A very sharp image after horizontal filter by my conv');
subplot(1, 3, 3);
imshow(double(of3l1)/256);
title('A very sharp image after horizontal filter by imfilter');

figure(4);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of4cl1)/256);
title('A very sharp image after vertical filter by my conv');
subplot(1, 3, 3);
imshow(double(of4l1)/256);
title('A very sharp image after vertical filter by imfilter');

figure(5);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of5cl1)/256);
title('A very sharp image after circular-shape filter by my conv');
subplot(1, 3, 3);
imshow(double(of5l1)/256);
title('A very sharp image after circular-shape filter by imfilter');

figure(6);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of6cl1)/256);
title('A very sharp image after 5 × 5 box filter by my conv');
subplot(1, 3, 3);
imshow(double(of6l1)/256);
title('A very sharp image after 5 × 5 box filter by imfilter');

figure(7);
subplot(1, 3, 1);
imshow(L/256);
title('A very sharp image');
subplot(1, 3, 2);
imshow(double(of7cl1)/256);
title('A very sharp image after 5 × 5 circular shape filter by my conv');
subplot(1, 3, 3);
imshow(double(of7l1)/256);
title('A very sharp image after 5 × 5 circular shape filter by imfilter');

figure(8);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of1ct1)/256);
title('An image degraded by noise after box filter by my conv');
subplot(1, 3, 3);
imshow(double(of1t1)/256);
title('An image degraded by noise after box filter by imfilter');

figure(9);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of2ct1)/256);
title('An image degraded by noise after weighted box filter by my conv');
subplot(1, 3, 3);
imshow(double(of2t1)/256);
title('An image degraded by noise after weighted box filter by imfilter');

figure(10);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of3ct1)/256);
title('An image degraded by noise after horizontal filter by my conv');
subplot(1, 3, 3);
imshow(double(of3t1)/256);
title('An image degraded by noise after horizontal filter by imfilter');

figure(11);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of4ct1)/256);
title('An image degraded by noise after vertical filter by my conv');
subplot(1, 3, 3);
imshow(double(of4t1)/256);
title('An image degraded by noise after vertical filter by imfilter');

figure(12);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of5ct1)/256);
title('An image degraded by noise after circular-shape filter by my conv');
subplot(1, 3, 3);
imshow(double(of5t1)/256);
title('An image degraded by noise after circular-shape filter by imfilter');

figure(13);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of6ct1)/256);
title('An image degraded by noise after 5 × 5 box filter by my conv');
subplot(1, 3, 3);
imshow(double(of6t1)/256);
title('An image degraded by noise after 5 × 5 box filter by imfilter');

figure(14);
subplot(1, 3, 1);
imshow(T/256);
title('An image degraded by noise');
subplot(1, 3, 2);
imshow(double(of7ct1)/256);
title('An image degraded by noise after 5 × 5 circular shape filter by my conv');
subplot(1, 3, 3);
imshow(double(of7t1)/256);
title('An image degraded by noise after 5 × 5 circular shape filter by imfilter');

of8a1 = imfilter(A, f8);
of8ca1 = my_conv(A, f8);
of8a = int8(imfilter(A, f8));
of8ca = int8(my_conv(A, f8));

of9a1 = imfilter(A, f9);
of9ca1 = my_conv(A, f9);
of9a = int8(imfilter(A, f9));
of9ca = int8(my_conv(A, f9));

result = isequal(of8a, of8ca);

if result
    disp('The results is equal. (4 × 4 Laplacian Filter, images to test on for enhancement and edge filters (A))');
else
    disp('The results is not equal. (4 × 4 Laplacian Filter, images to test on for enhancement and edge filters (A))');
end

result = isequal(of9a, of9ca);

if result
    disp('The results is equal. (8 × 8 Laplacian Filter, images to test on for enhancement and edge filters (A))');
else
    disp('The results is not equal. (8 × 8 Laplacian Filter, images to test on for enhancement and edge filters (A))');
end

of8b1 = imfilter(B, f8);
of8cb1 = my_conv(B, f8);
of8b = int8(imfilter(B, f8));
of8cb = int8(my_conv(B, f8));

of9b1 = imfilter(B, f9);
of9cb1 = my_conv(B, f9);
of9b = int8(imfilter(B, f9));
of9cb = int8(my_conv(B, f9));

result = isequal(of8b, of8cb);

if result
    disp('The results is equal. (4 × 4 Laplacian Filter, images to test on for enhancement and edge filters (B))');
else
    disp('The results is not equal. (4 × 4 Laplacian Filter, images to test on for enhancement and edge filters (B))');
end

result = isequal(of9b, of9cb);

if result
    disp('The results is equal. (8 × 8 Laplacian Filter, images to test on for enhancement and edge filters (B))');
else
    disp('The results is not equal. (8 × 8 Laplacian Filter, images to test on for enhancement and edge filters (B))');
end

figure(15);
subplot(1, 3, 1);
imshow(A/256);
title('Images to test on for enhancement and edge filters (A)');
subplot(1, 3, 2);
imshow(double(of8ca1)/256);
title('Images to test on for enhancement and edge filters (A) after 4 × 4 Laplacian Filter by my conv');
subplot(1, 3, 3);
imshow(double(of8a1)/256);
title('Images to test on for enhancement and edge filters (A) after 4 × 4 Laplacian Filter by imfilter');

figure(16);
subplot(1, 3, 1);
imshow(A/256);
title('Images to test on for enhancement and edge filters (A)');
subplot(1, 3, 2);
imshow(double(of9ca1)/256);
title('Images to test on for enhancement and edge filters (A) after 8 × 8 Laplacian Filter by my conv');
subplot(1, 3, 3);
imshow(double(of9a1)/256);
title('Images to test on for enhancement and edge filters (A) after 8 × 8 Laplacian Filter by imfilter');

figure(17);
subplot(1, 3, 1);
imshow(B/256);
title('Images to test on for enhancement and edge filters (B)');
subplot(1, 3, 2);
imshow(double(of8cb1)/256);
title('Images to test on for enhancement and edge filters (B) after 4 × 4 Laplacian Filter by my conv');
subplot(1, 3, 3);
imshow(double(of8b1)/256);
title('Images to test on for enhancement and edge filters (B) after 4 × 4 Laplacian Filter by imfilter');

figure(18);
subplot(1, 3, 1);
imshow(B/256);
title('Images to test on for enhancement and edge filters (B)');
subplot(1, 3, 2);
imshow(double(of9cb1)/256);
title('Images to test on for enhancement and edge filters (B) after 8 × 8 Laplacian Filter by my conv');
subplot(1, 3, 3);
imshow(double(of9b1)/256);
title('Images to test on for enhancement and edge filters (B) after 8 × 8 Laplacian Filter by imfilter');

of10a1 = imfilter(A, f10);
of10ca1 = my_conv(A, f10);
of10a = int8(imfilter(A, f10));
of10ca = int8(my_conv(A, f10));

of11a1 = imfilter(A, f11);
of11ca1 = my_conv(A, f11);
of11a = int8(imfilter(A, f11));
of11ca = int8(my_conv(A, f11));

of12a1 = imfilter(A, f12);
of12ca1 = my_conv(A, f12);
of12a = int8(imfilter(A, f12));
of12ca = int8(my_conv(A, f12));

of13a1 = imfilter(A, f13);
of13ca1 = my_conv(A, f13);
of13a = int8(imfilter(A, f13));
of13ca = int8(my_conv(A, f13));

result = isequal(of10a, of10ca);

if result
    disp('The results is equal. (Prewitt X Filter, images to test on for enhancement and edge filters (A))');
else
    disp('The results is not equal. (Prewitt X Filter, images to test on for enhancement and edge filters (A))');
end

result = isequal(of11a, of11ca);

if result
    disp('The results is equal. (Prewitt Y Filter, images to test on for enhancement and edge filters (A))');
else
    disp('The results is not equal. (Prewitt Y Filter, images to test on for enhancement and edge filters (A))');
end

result = isequal(of12a, of12ca);

if result
    disp('The results is equal. (Sobel X Filter, images to test on for enhancement and edge filters (A))');
else
    disp('The results is not equal. (Sobel X Filter, images to test on for enhancement and edge filters (A))');
end

result = isequal(of13a, of13ca);

if result
    disp('The results is equal. (Sobel Y Filter, images to test on for enhancement and edge filters (A))');
else
    disp('The results is not equal. (Sobel Y Filter, images to test on for enhancement and edge filters (A))');
end

figure(19);
subplot(1, 3, 1);
imshow(A/256);
title('Images to test on for enhancement and edge filters (A)');
subplot(1, 3, 2);
imshow(double(of10ca1)/256);
title('Images to test on for enhancement and edge filters (A) after Prewitt X Filter by my conv');
subplot(1, 3, 3);
imshow(double(of10a1)/256);
title('Images to test on for enhancement and edge filters (A) after Prewitt X Filter by imfilter');

figure(20);
subplot(1, 3, 1);
imshow(A/256);
title('Images to test on for enhancement and edge filters (A)');
subplot(1, 3, 2);
imshow(double(of11ca1)/256);
title('Images to test on for enhancement and edge filters (A) after Prewitt Y Filter by my conv');
subplot(1, 3, 3);
imshow(double(of11a1)/256);
title('Images to test on for enhancement and edge filters (A) after Prewitt Y Filter by imfilter');

figure(21);
subplot(1, 3, 1);
imshow(A/256);
title('Images to test on for enhancement and edge filters (A)');
subplot(1, 3, 2);
imshow(double(of12ca1)/256);
title('Images to test on for enhancement and edge filters (A) after Sobel X Filter by my conv');
subplot(1, 3, 3);
imshow(double(of12a1)/256);
title('Images to test on for enhancement and edge filters (A) after Sobel X Filter by imfilter');

figure(22);
subplot(1, 3, 1);
imshow(A/256);
title('Images to test on for enhancement and edge filters (A)');
subplot(1, 3, 2);
imshow(double(of13ca1)/256);
title('Images to test on for enhancement and edge filters (A) after Sobel Y Filter by my conv');
subplot(1, 3, 3);
imshow(double(of13a1)/256);
title('Images to test on for enhancement and edge filters (A) after Sobel Y Filter by imfilter');

of10b1 = imfilter(B, f10);
of10cb1 = my_conv(B, f10);
of10b = int8(imfilter(B, f10));
of10cb = int8(my_conv(B, f10));

of11b1 = imfilter(B, f11);
of11cb1 = my_conv(B, f11);
of11b = int8(imfilter(B, f11));
of11cb = int8(my_conv(B, f11));

of12b1 = imfilter(B, f12);
of12cb1 = my_conv(B, f12);
of12b = int8(imfilter(B, f12));
of12cb = int8(my_conv(B, f12));

of13b1 = imfilter(B, f13);
of13cb1 = my_conv(B, f13);
of13b = int8(imfilter(B, f13));
of13cb = int8(my_conv(B, f13));

result = isequal(of10b, of10cb);

if result
    disp('The results is equal. (Prewitt X Filter, images to test on for enhancement and edge filters (B))');
else
    disp('The results is not equal. (Prewitt X Filter, images to test on for enhancement and edge filters (B))');
end

result = isequal(of11b, of11cb);

if result
    disp('The results is equal. (Prewitt Y Filter, images to test on for enhancement and edge filters (B))');
else
    disp('The results is not equal. (Prewitt Y Filter, images to test on for enhancement and edge filters (B))');
end

result = isequal(of12b, of12cb);

if result
    disp('The results is equal. (Sobel X Filter, images to test on for enhancement and edge filters (B))');
else
    disp('The results is not equal. (Sobel X Filter, images to test on for enhancement and edge filters (B))');
end

result = isequal(of13b, of13cb);

if result
    disp('The results is equal. (Sobel Y Filter, images to test on for enhancement and edge filters (B))');
else
    disp('The results is not equal. (Sobel Y Filter, images to test on for enhancement and edge filters (B))');
end

figure(23);
subplot(1, 3, 1);
imshow(B/256);
title('Images to test on for enhancement and edge filters (B)');
subplot(1, 3, 2);
imshow(double(of10cb1)/256);
title('Images to test on for enhancement and edge filters (B) after Prewitt X Filter by my conv');
subplot(1, 3, 3);
imshow(double(of10b1)/256);
title('Images to test on for enhancement and edge filters (B) after Prewitt X Filter by imfilter');

figure(24);
subplot(1, 3, 1);
imshow(B/256);
title('Images to test on for enhancement and edge filters (B)');
subplot(1, 3, 2);
imshow(double(of11cb1)/256);
title('Images to test on for enhancement and edge filters (B) after Prewitt Y Filter by my conv');
subplot(1, 3, 3);
imshow(double(of11b1)/256);
title('Images to test on for enhancement and edge filters (B) after Prewitt Y Filter by imfilter');

figure(25);
subplot(1, 3, 1);
imshow(B/256);
title('Images to test on for enhancement and edge filters (B)');
subplot(1, 3, 2);
imshow(double(of12cb1)/256);
title('Images to test on for enhancement and edge filters (B) after Sobel X Filter by my conv');
subplot(1, 3, 3);
imshow(double(of12b1)/256);
title('Images to test on for enhancement and edge filters (B) after Sobel X Filter by imfilter');

figure(26);
subplot(1, 3, 1);
imshow(B/256);
title('Images to test on for enhancement and edge filters (B)');
subplot(1, 3, 2);
imshow(double(of13cb1)/256);
title('Images to test on for enhancement and edge filters (B) after Sobel Y Filter by my conv');
subplot(1, 3, 3);
imshow(double(of13b1)/256);
title('Images to test on for enhancement and edge filters (B) after Sobel Y Filter by imfilter');

function [output] = my_conv(image, filter)
    
    [rows, cols] = size(image);
    
    [frows, fcols] = size(filter);
    
    filterCenterX = floor(fcols / 2) + 1;
    filterCenterY = floor(frows / 2) + 1;

    output = zeros(rows, cols);

    for x = 1:cols
        for y = 1:rows

            foutput = 0;
           
            for filterX = 1:fcols
                for filterY = 1:frows
                 
                    inputX = x + filterX - filterCenterX;
                    inputY = y + filterY - filterCenterY;

                    if inputX >= 1 && inputX <= cols && inputY >= 1 && inputY <= rows
                        foutput = foutput + image(inputY, inputX) * filter(filterY, filterX);
                    end
                end
            end

            output(y, x) = foutput;
        end
    end
end