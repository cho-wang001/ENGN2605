p3 = imread('ImageA.jpg');
p31 = imread('ImageA.jpg');

p32 = p3 + 25;

p33 = p3 + 50;

figure(1)
subplot(1,3,1)
imshow(p31)
title('Original image')
subplot(1,3,2)
imshow(p32)
title('Increase bright image I0=25')
subplot(1,3,3)
imshow(p33)
title('Increase bright image I0=50')