clc;
clear;

C = double(imread('salt_and_pepper_coins.png'));
T = double(imread('traffic.jpg'));


ak5 = (1/25)./ones(5, 5);
gk515 = gauss_kernel(5, 1.5);

imcak5 = imfilter(C, ak5);
imcgk515 = imfilter(C, gk515);

figure(1);
subplot(1, 3, 1);
imshow(C/256);
title('salt and pepper coins');
subplot(1, 3, 2);
imshow(imcak5/256);
title('salt and pepper coins after 5 × 5 box Filter filter. (Window Size of 5)');
subplot(1, 3, 3);
imshow(imcgk515/256);
title('salt and pepper coins after gaussian filter. (Window Size of 5, σ = 1.5)');

imtmed3 = my_medfilter(T, 3);
imtmax3 = my_maxfilter(T, 3);
imtmin3 = my_minfilter(T, 3);

figure(2);
subplot(2, 2, 1);
imshow(T/256);
title('The noisy traffic image');
subplot(2, 2, 2);
imshow(imtmed3/256);
title('The noisy traffic image after median filter. (Window Size of 3)');
subplot(2, 2, 3);
imshow(imtmax3/256);
title('The noisy traffic image after max filter. (Window Size of 3)');
subplot(2, 2, 4);
imshow(imtmin3/256);
title('The noisy traffic image after min filter. (Window Size of 3)');


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

function [output] = my_minfilter(image, filter_size)  
      
    [rows, cols] = size(image);
    frows = filter_size;
    fcols = filter_size;
      
    filterCenterX = floor(fcols / 2) + 1;  
    filterCenterY = floor(frows / 2) + 1;  
  
    output = zeros(rows, cols);  
  
    for x = 1:cols  
        for y = 1:rows  
  
            foutput = [];  
             
            for filterX = 1:fcols  
                for filterY = 1:frows  
                   
                    inputX = x + filterX - filterCenterX;  
                    inputY = y + filterY - filterCenterY;  
  
                    if inputX >= 1 && inputX <= cols && inputY >= 1 && inputY <= rows  
                        foutput = [foutput; image(inputY, inputX)];  
                    end  
                end  
            end  
     
            output(y, x) = min(foutput);  
        end  
    end  
end

function [output] = my_maxfilter(image, filter_size)  
      
    [rows, cols] = size(image);
    frows = filter_size;
    fcols = filter_size;
      
    filterCenterX = floor(fcols / 2) + 1;  
    filterCenterY = floor(frows / 2) + 1;  
  
    output = zeros(rows, cols);  
  
    for x = 1:cols  
        for y = 1:rows  
  
            foutput = [];  
             
            for filterX = 1:fcols  
                for filterY = 1:frows  
                   
                    inputX = x + filterX - filterCenterX;  
                    inputY = y + filterY - filterCenterY;  
  
                    if inputX >= 1 && inputX <= cols && inputY >= 1 && inputY <= rows  
                        foutput = [foutput; image(inputY, inputX)];  
                    end  
                end  
            end  
     
            output(y, x) = max(foutput);  
        end  
    end  
end

function [output] = my_medfilter(image, filter_size)  
      
    [rows, cols] = size(image);
    frows = filter_size;
    fcols = filter_size;
      
    filterCenterX = floor(fcols / 2) + 1;  
    filterCenterY = floor(frows / 2) + 1;  
  
    output = zeros(rows, cols);  
  
    for x = 1:cols  
        for y = 1:rows  
  
            foutput = [];  
             
            for filterX = 1:fcols  
                for filterY = 1:frows  
                   
                    inputX = x + filterX - filterCenterX;  
                    inputY = y + filterY - filterCenterY;  
  
                    if inputX >= 1 && inputX <= cols && inputY >= 1 && inputY <= rows  
                        foutput = [foutput; image(inputY, inputX)];  
                    end  
                end  
            end  
   
            foutput = sort(foutput);  
            output(y, x) = foutput(floor(length(foutput)/2) + 1);  
        end  
    end  
end
