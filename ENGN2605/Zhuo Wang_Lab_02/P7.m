clc;
clear;


load psf_eye_chart.mat;
load psf_circuit_board.mat;

E = double(imread('eyechart_degraded.jpg'));
C = double(imread('Circuit_Board_degraded.png'));

ep = lucy_richardson(E, psf_eye_chart);
ep1 = deconvlucy(E, psf_eye_chart);

cp = lucy_richardson(C, psf_circuit_board);
cp1 = deconvlucy(C, psf_circuit_board);

epm3 = my_medfilter(E, 3);
cpm3 = my_medfilter(C, 3);

epm7 = my_medfilter(E, 7);
cpm7 = my_medfilter(C, 7);

rep = ep;

f1 = [0, 1, 0; 1, -4, 1; 0, 1, 0];
f2 = [1, 1, 1; 1, -8, 1; 1, 1, 1];
f3 = (1/9).*ones(3,3);

epf1 = imfilter(E, f1);
epf2 = imfilter(E, f2);

cpf1 = imfilter(C, f1);
cpf2 = imfilter(C, f2);

rep1 = my_maxfilter(rep, 5);


figure(1)
subplot(1,2,1);
imshow(ep/256);
title('Lucy-Richardson implementation');
subplot(1,2,2);
imshow(ep1/256);
title('deconvlucy implementation');

figure(2);
subplot(1,2,1);
imshow(cp/256);
title('Lucy-Richardson implementation');
subplot(1,2,2);
imshow(cp1/256);
title('deconvlucy implementation');

figure(3);
subplot(2,2,1);
imshow(E/256);
title('Degraded image');
subplot(2,2,2);
imshow(ep/256);
title('Lucy-Richardson implementation');
subplot(2,2,3);
imshow(epm3/256);
title('Degraded image after median filter. (Window Size of 3)');
subplot(2,2,4);
imshow(epm7/256);
title('Degraded image after median filter. (Window Size of 7)');

figure(4);
subplot(2,2,1);
imshow(C/256);
title('Degraded image');
subplot(2,2,2);
imshow(cp/256);
title('Lucy-Richardson implementation');
subplot(2,2,3);
imshow(cpm3/256);
title('Degraded image after median filter. (Window Size of 3)');
subplot(2,2,4);
imshow(cpm7/256);
title('Degraded image after median filter. (Window Size of 7)');

figure(5);
subplot(2,2,1);
imshow(E/256);
title('Degraded image');
subplot(2,2,2);
imshow(ep/256);
title('Lucy-Richardson implementation');
subplot(2,2,3);
imshow(epf1);
title('Degraded image after  4 × 4 Laplacian');
subplot(2,2,4);
imshow(epf2);
title('Degraded image after 8 × 8 Laplacian');

figure(6);
subplot(2,2,1);
imshow(C/256);
title('Degraded image');
subplot(2,2,2);
imshow(cp/256);
title('Lucy-Richardson implementation');
subplot(2,2,3);
imshow(cpf1/256);
title('Degraded image after  4 × 4 Laplacian');
subplot(2,2,4);
imshow(cpf2/256);
title('Degraded image after 8 × 8 Laplacian');

figure(7);
imshow(rep1/256);
title('Result');

function [restored_image] = lucy_richardson(degraded_image, psf)

    estimate = degraded_image;
    
    for iteration = 1:100
        reblurred_image = conv2(estimate, psf, 'same');
        error_image = degraded_image ./ reblurred_image;
        correction_image = filter2(psf,error_image);
        estimate = estimate .* correction_image;
    end

    restored_image = estimate;
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




