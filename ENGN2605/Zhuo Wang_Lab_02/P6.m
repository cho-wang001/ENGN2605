clc;
clear;

L = double(imread('Lenna.png'));

f1 = (1/25).*ones(5,5);
f2 = (1/21).*[0, 1, 1, 1, 0; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 1, 1, 1, 1, 1; 0, 1, 1, 1, 0];
fs75 = fspecial('gaussian', 7, 0.5);
fs175 = fspecial('gaussian', 17, 0.5);
fs3520 = fspecial('gaussian', 35, 2.0);

tic
lcsf1 = conv_separable(L, f1);
toc
disp('5 × 5 Box Filter by Separable Convolution');

tic
lmcf1 = my_conv(L, f1);
toc
disp('5 × 5 Box Filter by Traditional Convolution');

tic
lcsf2 = conv_separable(L, f2);
toc
disp('5 × 5 circular shape filter cannot by Separable Convolution');

tic
lmcf2 = my_conv(L, f2);
toc
disp('5 × 5 circular shape filter by Traditional Convolution');

tic
lcsfs75 = conv_separable(L, fs75);
toc
disp('7 × 7 Gaussian with sigma of 0.5 by Separable Convolution');

tic
lmcfs75 = my_conv(L, fs75);
toc
disp('7 × 7 Gaussian with sigma of 0.5 by Traditional Convolution');

tic
lcsfs175 = conv_separable(L, fs175);
toc
disp('17 × 17 Gaussian with sigma of 0.5 by Separable Convolution');

tic
lmcfs175 = my_conv(L, fs175);
toc
disp('17 × 17 Gaussian with sigma of 0.5 by Traditional Convolution');

tic
lcsfs3520 = conv_separable(L, fs3520);
toc
disp('35 × 35 Gaussian with sigma of 2.0 by Separable Convolution');

tic
lmcfs3520 = my_conv(L, fs3520);
toc
disp('35 × 35 Gaussian with sigma of 2.0 by Traditional Convolution');



figure(1);
imshow(lcsf1/256);
title('5 × 5 Box Filter by Separable Convolution');

figure(2);
imshow(lmcf1/256);
title('5 × 5 Box Filter by Traditional Convolution');

figure(3);
imshow(lcsf2/256);
title('5 × 5 circular shape filter by Separable Convolution');

figure(4);
imshow(lmcf2/256);
title('5 × 5 circular shape filter by Traditional Convolution');

figure(5);
imshow(lcsfs75/256);
title('7 × 7 Gaussian with sigma of 0.5 by Separable Convolution');

figure(6);
imshow(lmcfs75/256);
title('7 × 7 Gaussian with sigma of 0.5 by Traditional Convolution');

figure(7);
imshow(lcsfs175/256);
title('17 × 17 Gaussian with sigma of 0.5 by Separable Convolution');

figure(8);
imshow(lmcfs175/256);
title('17 × 17 Gaussian with sigma of 0.5 by Traditional Convolution');

figure(9);
imshow(lcsfs3520/256);
title('35 × 35 Gaussian with sigma of 2.0 by Separable Convolution');

figure(10);
imshow(lmcfs3520/256);
title('35 × 35 Gaussian with sigma of 2.0 by Traditional Convolution');




function [output] = conv_separable(image, filter)
    [U, S, V] = svd(filter);
    k = 1;
    row_vector = -(U(:, k) * sqrt(S(k, k)));
    column_vector = -(sqrt(S(k, k)) * V(:, k)');

    yz = round(row_vector*column_vector, 3);
    test = round(filter, 3);

    r = isequal(yz, test);
    if r == 0
        disp('Filter inseparable. Return 0.');
        output = 0;
    else
        [image_height, image_width] = size(image);
        col_vector_length = length(column_vector);
        row_vector_length = length(row_vector);

        filterCenter = floor(col_vector_length / 2) + 1;

        pimage = zeros(image_height + filterCenter, image_width + filterCenter);

        pimage(1:image_height, 1:image_width) = image;

        buffer = zeros(image_height+filterCenter, image_width+filterCenter);
        
        for col = 1:image_width+filterCenter
            for row = 1:image_height+filterCenter
    
                col_result = 0;
                for m = 1:col_vector_length
                    y = row - m + 1;
                    if y >= 1
                        col_result = col_result + column_vector(m) * pimage(y, col);
                    end
                end
                buffer(row, col) = col_result;
            end
        end
        
        output = zeros(image_height+filterCenter, image_width+filterCenter);
    
        for row = 1:image_height
            for col = 1:image_width
    
                row_result = 0;
                for n = 1:row_vector_length
                    x = col - n + 1;
                    if x >= 1
                        row_result = row_result + row_vector(n) * buffer(row, x);
                    end
                end
                output(row, col) = row_result;
            end
        end
    end
end

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



