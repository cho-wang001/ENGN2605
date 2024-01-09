p6 = imread('white_blood_cells.png');

p6 = rgb2gray(p6);

p61 = ones(69,92);
o = 1;
k = 1;


for i=1:7:480
    for j=1:7:640
        p61(o,k) = p6(i,j);
        k = k+1;
    end
    o = o+1;
    k = 1;
end

p61 = (p61./256);

figure(1)
subplot(1,2,1)
imshow(p6)
title('Original image')
subplot(1,2,2)
imshow(p61)
title('Down Sampling image')
