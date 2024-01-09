p5 = imread('white_blood_cells.png');

p51 = rgb2gray(p5);
p52 = rgb2gray(p5);

for i=1:480
    for j=1:640
        if (0 <= p51(i,j)) && (p51(i,j) <= 15)
            p51(i,j) = 8;
        elseif (16 <= p51(i,j)) && (p51(i,j) <= 31)
            p51(i,j) = 24;
        elseif (32 <= p51(i,j)) && (p51(i,j) <= 47)
            p51(i,j) = 40;
        elseif (48 <= p51(i,j)) && (p51(i,j) <= 63)
            p51(i,j) = 56;     
        elseif (64 <= p51(i,j)) && (p51(i,j) <= 79)
            p51(i,j) = 72;
        elseif (80 <= p51(i,j)) && (p51(i,j) <= 95)
            p51(i,j) = 88;
        elseif (96 <= p51(i,j)) && (p51(i,j) <= 111)
            p51(i,j) = 104;      
        elseif (112 <= p51(i,j)) && (p51(i,j) <= 127)
            p51(i,j) = 120;
        elseif (128 <= p51(i,j)) && (p51(i,j) <= 143)
            p51(i,j) = 136;
        elseif (144 <= p51(i,j)) && (p51(i,j) <= 159)
            p51(i,j) = 152;     
        elseif (160 <= p51(i,j)) && (p51(i,j) <= 175)
            p51(i,j) = 168;
        elseif (176 <= p51(i,j)) && (p51(i,j) <= 191)
            p51(i,j) = 184;
        elseif (192 <= p51(i,j)) && (p51(i,j) <= 207)
            p51(i,j) = 200;
        elseif (208 <= p51(i,j)) && (p51(i,j) <= 223)
            p51(i,j) = 216;     
        elseif (224 <= p51(i,j)) && (p51(i,j) <= 239)
            p51(i,j) = 232;
        elseif (240 <= p51(i,j)) && (p51(i,j) <= 255)
            p51(i,j) = 248;           
        end
    end
end

figure(1)
subplot(1,3,1)
imshow(p5)
title('Original image')
subplot(1,3,2)
imshow(p52)
title('Gray original image')
subplot(1,3,3)
imshow(p51)
title('Reduce to 16 levels image')