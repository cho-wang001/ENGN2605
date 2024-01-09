p1 = imread('white_blood_cells.png');
grayp1 = rgb2gray(p1);
gp1 = rgb2gray(p1);

T0 = 241;

for i=1:480
    for j=1:640
        if grayp1(i,j) >= T0
            grayp1(i,j) = 255;
        else
            grayp1(i,j) = 0;
        end
    end
end

figure(1)
subplot(1,3,1)
imshow(p1)
title('Original image')
subplot(1,3,2)
imshow(gp1)
title('Gray original image')
subplot(1,3,3)
imshow(grayp1)
title('Threshold segmentation image')