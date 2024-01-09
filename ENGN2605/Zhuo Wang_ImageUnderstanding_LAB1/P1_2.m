p2 = imread('snowman.jpg');
p21 = imread('snowman.jpg');

T0 = 100;
T1 = 170;

for i=1:200
    for j=1:255
        if p2(i,j) >= T0 && p2(i,j) < T1
            p2(i,j) = 255;
        else
            p2(i,j) = 0;
        end
    end
end


figure(1)
subplot(1,2,1)
imshow(p21)
title('Original image')
subplot(1,2,2)
imshow(p2)
title('Threshold segmentation image')