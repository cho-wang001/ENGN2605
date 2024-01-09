clc;
clear;
close all;

%Run the 'vl_setup' file before running it

I1 = imread("Images\problem_5\1.png");
I2 = imread("Images\problem_5\2.png");
I3 = imread("Images\problem_5\3.png");
I4 = imread("Images\problem_5\4.png");
I5 = imread("Images\problem_5\5.png");

i1 = imread("Images\problem_5\1.png");


peakthres = 10;
edgethres = 20;

I1 = single(im2gray(I1));
I2 = single(im2gray(I2));
I3 = single(im2gray(I3));
I4 = single(im2gray(I4));
I5 = single(im2gray(I5));

[f1,d1] = vl_sift(I1,'PeakThresh',peakthres,'edgethresh',edgethres);
[f2,d2] = vl_sift(I2,'PeakThresh',peakthres,'edgethresh',edgethres);


pc = [f1(1,:);f1(2,:);ones(size(f1(1,:)))];

h = load("Images\problem_5\H_1to2.txt");

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

cr12 = [];
[~,num] = size(r_f);
[~,num1] = size(f2);
index = [];
j = 1;

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f2(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 20
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr12 = [cr12,f2(1:2,min_ind)];
        index = [index, j];
        j = j+1;
    end
end



[f3,d3] = vl_sift(I3,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [cr12(1,:);cr12(2,:);ones(size(cr12(1,:)))];

h = load("Images\problem_5\H_2to3.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I3);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr23 = [];
[~,num] = size(r_f);
[~,num1] = size(f3);
index2 = [];
j = 1;

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f3(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 20
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr23 = [cr23,f3(1:2,min_ind)];
        index2 = [index2, j];
        j = j+1;
    end
end



[f4,d4] = vl_sift(I4,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [cr23(1,:);cr23(2,:);ones(size(cr23(1,:)))];

h = load("Images\problem_5\H_3to4.txt");

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

cr34 = [];
[~,num] = size(r_f);
[~,num1] = size(f4);
index3 = [];
j = 1;

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f4(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 20
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr34 = [cr34,f4(1:2,min_ind)];
        index3 = [index3, j];
        j = j + 1;
    end
end


[f5,d5] = vl_sift(I5,'PeakThresh',peakthres,'edgethresh',edgethres);

pc = [cr34(1,:);cr34(2,:);ones(size(cr34(1,:)))];

h = load("Images\problem_5\H_4to5.txt");

c_pc = h*pc;

w = c_pc(3, :);

r = c_pc(1:2,:)./w;

[~,num] = size(r);
[height,width] = size(I5);
r_f = [];

for i = 1:num
    if r(1,i) > 0 && r(1,i) <= width && r(2,i) > 0 && r(2,i) <= height
        r_f(:, end+1) = r(:, i);
    end
end

cr45 = [];
[~,num] = size(r_f);
[~,num1] = size(f5);
index4 = [];
j = 1;

for i = 1:num
    p1 = r_f(:,i);
    min_ind = 0;
    min_dist = Inf;
    for m = 1:num1
        p2 = f5(1:2,m);
        dist = norm(p1-p2);
        if dist < min_dist && dist <= 20
            min_dist = dist;
            min_ind = m;
        end
    end
    if min_ind ~= 0
        cr45 = [cr45,f5(1:2,min_ind)];
        index4 = [index4, j];
        j = j + 1;
    end
end


common_ind1 = ismember(index, index2);
common_ind2 = ismember(index, index3);
common_ind3 = ismember(index, index4);

common_ind = common_ind1 & common_ind2 & common_ind3;
  
figure();
imshow(i1);
hold on  

scatter(cr12(1, :), cr12(2, :), 'filled', 'r');
scatter(cr23(1, :), cr23(2, :), 'filled', 'b');
scatter(cr34(1, :), cr34(2, :), 'filled', 'y');
scatter(cr45(1, :), cr45(2, :), 'filled', 'g');

for i = 1:size(common_ind, 2)  
    if common_ind(i)

        p1 = [cr12(1, i), cr12(2, i)];  
        p2 = [cr23(1, i), cr23(2, i)]; 
        p3 = [cr34(1, i), cr34(2, i)];
        p4 = [cr45(1, i), cr45(2, i)]; 

        plot([p1(1), p2(1)], [p1(2), p2(2)], 'g-');
        plot([p2(1), p3(1)], [p2(2), p3(2)], 'g-');
        plot([p3(1), p4(1)], [p3(2), p4(2)], 'g-');

    end  
end  
hold off

% figure();
% scatter3(cr12(1,:), cr12(2,:), 10); hold on;
% scatter3(cr23(1,:), cr23(2,:), 20);
% scatter3(cr34(1,:), cr34(2,:), 30);
% scatter3(cr45(1,:), cr45(2,:), 40);
% 
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% 
% 
% view(3);

