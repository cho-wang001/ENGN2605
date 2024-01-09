clc
clear

p71 = double(rgb2gray(imread("1.jpg")));
p72 = double(rgb2gray(imread("2.jpg")));
p73 = double(rgb2gray(imread("3.jpg")));
p74 = double(rgb2gray(imread("4.jpg")));
p75 = double(rgb2gray(imread("5.jpg")));
p76 = double(rgb2gray(imread("6.jpg")));
p77 = double(rgb2gray(imread("7.jpg")));
p78 = double(rgb2gray(imread("8.jpg")));
p79 = double(rgb2gray(imread("9.jpg")));
p710 = double(rgb2gray(imread("10.jpg")));
p711 = double(rgb2gray(imread("11.jpg")));
p712 = double(rgb2gray(imread("12.jpg")));
p713 = double(rgb2gray(imread("13.jpg")));
p714 = double(rgb2gray(imread("14.jpg")));
p715 = double(rgb2gray(imread("15.jpg")));
p716 = double(rgb2gray(imread("16.jpg")));
p717 = double(rgb2gray(imread("17.jpg")));
p718 = double(rgb2gray(imread("18.jpg")));
p719 = double(rgb2gray(imread("19.jpg")));
p720 = double(rgb2gray(imread("20.jpg")));

imageList = {p71, p72, p73, p74, p75, p76, p77, p78, p79, p710, p711, p712, p713, p714, p715, p716, p717, p718, p719, p720};

averageImage = mean(cat(4, imageList{:}), 4);
figure(1)
imshow(averageImage/256);
title('averageImage')

imageStack = cat(4, imageList{:});
stdImage = std(imageStack, 0, 4);
figure(2)
imshow(stdImage, []);
title('stdImage')

alist = [];

for i=1:20
    ac = imageList{i} - averageImage;
    a = sum(ac(:));
    alist = [alist; a];
end

absalist = abs(alist);

[maxAbsValue, maxAbsIndex] = max(absalist(:));

disp(['The absolute value of the maximum difference between the mean images：', num2str(maxAbsValue)]);
disp(['Position in the list (which image)：', num2str(maxAbsIndex)]);

hlist = [];

for i=1:20
    h = imageList{i}(666,666);
    hlist = [hlist; h];
end

figure(3);
hist(hlist);
