clc;
clear;

C = double(imread('salt_and_pepper_coins.png'));

cmin3 = my_minfilter(C, 3);
cmin5 = my_minfilter(C, 5);
cmin7 = my_minfilter(C, 7);

cmax3 = my_maxfilter(C, 3);
cmax5 = my_maxfilter(C, 5);
cmax7 = my_maxfilter(C, 7);

cmed3 = my_medfilter(C, 3);
cmed5 = my_medfilter(C, 5);
cmed7 = my_medfilter(C, 7);

figure(1);
subplot(2, 2, 1);
imshow(C/256);
title('salt and pepper coins');
subplot(2, 2, 2);
imshow(cmin3/256);
title('salt and pepper coins after min filter. (Window Size of 3)');
subplot(2, 2, 3);
imshow(cmin5/256);
title('salt and pepper coins after min filter. (Window Size of 5)');
subplot(2, 2, 4);
imshow(cmin7/256);
title('salt and pepper coins after min filter. (Window Size of 7)');

figure(2);
subplot(2, 2, 1);
imshow(C/256);
title('salt and pepper coins');
subplot(2, 2, 2);
imshow(cmax3/256);
title('salt and pepper coins after max filter. (Window Size of 3)');
subplot(2, 2, 3);
imshow(cmax5/256);
title('salt and pepper coins after max filter. (Window Size of 5)');
subplot(2, 2, 4);
imshow(cmax7/256);
title('salt and pepper coins after max filter. (Window Size of 7)');

figure(3);
subplot(2, 2, 1);
imshow(C/256);
title('salt and pepper coins');
subplot(2, 2, 2);
imshow(cmed3/256);
title('salt and pepper coins after median filter. (Window Size of 3)');
subplot(2, 2, 3);
imshow(cmed5/256);
title('salt and pepper coins after median filter. (Window Size of 5)');
subplot(2, 2, 4);
imshow(cmed7/256);
title('salt and pepper coins after median filter. (Window Size of 7)');

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