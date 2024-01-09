clc;
clear;

ei = magic(10);
f = ones(3,3);

o = imfilter(ei,f);
o1 = my_conv(ei,f);

o
o1

result = isequal(o, o1);

if result
    disp('The results is equal.');
else
    disp('The results is not equal.');
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

