clc;
clear;
close all;


%Run the 'vl_setup' file before running it
I1 = imread("Images\problem_3_and_4\set1\img1.png");
I2 = imread("Images\problem_3_and_4\set1\img2.png");

peakthres = 10;
edgethres = 20;

I1 = single(im2gray(I1));
I2 = single(im2gray(I2));

[f1,d1] = vl_sift(I1,'PeakThresh',peakthres,'edgethresh',edgethres);
[f2,d2] = vl_sift(I2,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [f1(1,:);f1(2,:);ones(size(f1(1,:)))];

h = load("Images\problem_3_and_4\set1\H_1to2.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I2);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr = [];
[~,num] = size(r_f);
[~,num1] = size(f2);

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f2(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 1
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr = [cr,f2(1:2,min_ind)];
    end
end
[~,num2] = size(cr);
rp = num/num2;
disp(rp);

figure();
subplot(1, 3, 1);
imshow(I1/256);
hold on
scatter(f1(1,:),f1(2,:),'r+');

subplot(1, 3, 2);
imshow(I2/256);
hold on
scatter(f2(1,:),f2(2,:),'bo');

subplot(1, 3, 3);
imshow(I2/256);
hold on
scatter(cr(1,:), cr(2,:),'y*');





I3 = imread("Images\problem_3_and_4\set2\img1.png");
I4 = imread("Images\problem_3_and_4\set2\img2.png");

peakthres = 10;
edgethres = 20;

I3 = single(im2gray(I3));
I4 = single(im2gray(I4));

[f3,d3] = vl_sift(I3,'PeakThresh',peakthres,'edgethresh',edgethres);
[f4,d4] = vl_sift(I4,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [f3(1,:);f3(2,:);ones(size(f3(1,:)))];

h = load("Images\problem_3_and_4\set2\H_1to2.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I4);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr = [];
[~,num] = size(r_f);
[~,num1] = size(f4);

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f4(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 1
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr = [cr,f4(1:2,min_ind)];
    end
end
[~,num2] = size(cr);
rp = num/num2;
disp(rp);

figure();
subplot(1, 3, 1);
imshow(I3/256);
hold on
scatter(f3(1,:),f3(2,:),'r+');

subplot(1, 3, 2);
imshow(I4/256);
hold on
scatter(f4(1,:),f4(2,:),'bo');

subplot(1, 3, 3);
imshow(I4/256);
hold on
scatter(cr(1,:), cr(2,:),'y*');





I5 = imread("Images\problem_3_and_4\set3\img1.png");
I6 = imread("Images\problem_3_and_4\set3\img2.png");

peakthres = 10;
edgethres = 20;

I5 = single(im2gray(I5));
I6 = single(im2gray(I6));

[f5,d5] = vl_sift(I5,'PeakThresh',peakthres,'edgethresh',edgethres);
[f6,d6] = vl_sift(I6,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [f5(1,:);f5(2,:);ones(size(f5(1,:)))];

h = load("Images\problem_3_and_4\set3\H_1to2.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I6);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr = [];
[~,num] = size(r_f);
[~,num1] = size(f6);

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f6(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 2
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr = [cr,f6(1:2,min_ind)];
    end
end
[~,num2] = size(cr);
rp = num/num2;
disp(rp);

figure();
subplot(1, 3, 1);
imshow(I5/256);
hold on
scatter(f5(1,:),f5(2,:),'r+');

subplot(1, 3, 2);
imshow(I6/256);
hold on
scatter(f6(1,:),f6(2,:),'bo');

subplot(1, 3, 3);
imshow(I6/256);
hold on
scatter(cr(1,:), cr(2,:),'y*');




I7 = imread("Images\problem_3_and_4\set4\img1.png");
I8 = imread("Images\problem_3_and_4\set4\img2.png");

peakthres = 10;
edgethres = 20;

I7 = single(im2gray(I7));
I8 = single(im2gray(I8));

[f7,d7] = vl_sift(I7,'PeakThresh',peakthres,'edgethresh',edgethres);
[f8,d8] = vl_sift(I8,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [f7(1,:);f7(2,:);ones(size(f7(1,:)))];

h = load("Images\problem_3_and_4\set4\H_1to2.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I8);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr = [];
[~,num] = size(r_f);
[~,num1] = size(f8);

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f8(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 1
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr = [cr,f8(1:2,min_ind)];
    end
end
[~,num2] = size(cr);
rp = num/num2;
disp(rp);

figure();
subplot(1, 3, 1);
imshow(I7/256);
hold on
scatter(f7(1,:),f7(2,:),'r+');

subplot(1, 3, 2);
imshow(I8/256);
hold on
scatter(f8(1,:),f8(2,:),'bo');

subplot(1, 3, 3);
imshow(I8/256);
hold on
scatter(cr(1,:), cr(2,:),'y*');



%Use my corner detector in problem 1 and repeat steps 1 to 4.

cd = corner_detector(I1);
cd1 = corner_detector(I2);

[row, col] = find(cd == 1);
f9 = [col, row];
f9 = f9';


[row, col] = find(cd1 == 1);
f10 = [col, row];
f10 = f10';

pc = [f9(1,:);f9(2,:);ones(size(f9(1,:)))];

h = load("Images\problem_3_and_4\set1\H_1to2.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I2);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr = [];
[~,num] = size(r_f);
[~,num1] = size(f10);

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f10(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 2
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr = [cr,f10(1:2,min_ind)];
    end
end
[~,num2] = size(cr);
rp = num/num2;
disp(rp);

figure();
subplot(1, 3, 1);
imshow(I1/256);
hold on
scatter(f9(1,:),f9(2,:),'r+');

subplot(1, 3, 2);
imshow(I2/256);
hold on
scatter(f10(1,:),f10(2,:),'bo');

subplot(1, 3, 3);
imshow(I2/256);
hold on
scatter(cr(1,:), cr(2,:),'y*');


function corners = corner_detector(image)  

    if size(image, 3) == 3
        image = im2double(rgb2gray(image));
    else
        image = im2double(image);
    end

    sigma1 = 0.7;
    % For the window size of the Gaussians, use ⌈6σ⌉. Make it an odd number by adding 1 to it if necessary.
    ws = int8(6.*sigma1)+1;
    ws = double(ws);
    [x, y] = meshgrid(-(ws-1)/2:(ws-1)/2, -(ws-1)/2:(ws-1)/2);

    Gx = exp(-(x.^2 + y.^2) / (2*sigma1^2)) .* (-x / (2*pi*sigma1^4));
    fx = conv2(image, Gx, 'same');

    Gy = exp(-(x.^2 + y.^2) / (2*sigma1^2)) .* (-y / (2*pi*sigma1^4));
    fy = conv2(image, Gy, 'same');
    
    fxp2 = fx.^2;

    fyp2 = fy.^2;
    
    fxfy = fx.*fx;

    sigma2 = 2.0;

    sfxp2 = imgaussfilt(fxp2, sigma2);
    sfyp2 = imgaussfilt(fyp2, sigma2);
    sfxfy = imgaussfilt(fxfy, sigma2);
    
    M = [sfxp2, sfxfy; sfxfy, sfyp2];
    
    a = 0.04;

    R = (sfxp2 .* sfyp2 - sfxfy.^2) - 0.04 * (sfxp2 + sfyp2).^2;
    
    [rows, cols] = size(R);

    corner_markers = zeros(rows, cols);

    neighborhood_offsets = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 0; 0, 1; 1, -1; 1, 0; 1, 1];
    
    [maxr, ~] = max(R(:));

    R0 = maxr.*0.01;

    for x = 2:rows-1
        for y = 2:cols-1
            current_R = R(x, y);
        
            is_local_max = true;
        
            for i = 1:9

                neighbor_x = x + neighborhood_offsets(i, 1);
                neighbor_y = y + neighborhood_offsets(i, 2);
            
                neighbor_R = R(neighbor_x, neighbor_y);
            
                if neighbor_R > current_R
                    is_local_max = false;
                    break;
                end
            end
            
            if is_local_max && (current_R > R0)
                corner_markers(x, y) = 1;
            end

        end
    end

    corners = corner_markers;

end




