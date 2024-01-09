clc;
clear;
close all;
format long;

%Run the 'vl_setup' file before running it
I1 = rgb2gray(imread('Q2_data\img_pair_1\1.png'));
I2 = rgb2gray(imread('Q2_data\img_pair_1\2.png'));
I3 = rgb2gray(imread('Q2_data\img_pair_2\1.png'));
I4 = rgb2gray(imread('Q2_data\img_pair_2\2.png'));

k1 = load("Q2_data\img_pair_1\Calibration_Matrix.mat").K;
k2 = load("Q2_data\img_pair_2\Calibration_Matrix.mat").K;

e12 = load("Q2_data\img_pair_1\E_1to2.mat").E;
e34 = load("Q2_data\img_pair_2\E_1to2.mat").E;

i1 = single(I1);
i2 = single(I2);
i3 = single(I3);
i4 = single(I4);

peak_thresh = 1;
edgethresh = 5;

[f1, d1] = vl_sift(i1,'PeakThresh',peak_thresh,'EdgeThresh',edgethresh);
[f2, d2] = vl_sift(i2,'PeakThresh',peak_thresh,'EdgeThresh',edgethresh);
[f3, d3] = vl_sift(i3,'PeakThresh',peak_thresh,'EdgeThresh',edgethresh);
[f4, d4] = vl_sift(i4,'PeakThresh',peak_thresh,'EdgeThresh',edgethresh);

[i1, i2, s1, s2] = find_matches(d1, d2,"sift", 0.75);
mm = [i1, s1];
smm = sortrows(mm, 3);
ii1 = smm(:, 1:2);
ss1 = smm(:, 3);
mm = smm(1:20, :);
i1 = mm(:, 1:2);
s1 = mm(:, 3);
[c1, c2] = coordinates_find_sift(i1, s1, f1, f2);
[cc1, cc2] = coordinates_find_sift(ii1, ss1, f1, f2);

c1 = [c1, ones(20,1)];
c2 = [c2, ones(20,1)];

a = [];
b = [];
c = [];

invk1 = inv(k1);

figure();
I12 = [I1, I2];
imshow(I12);
hold on;

s = size(I12, 2)./2;

for i = 1:20
    c13d = invk1*[c1(i, 1);c1(i, 2);1];

    A = c13d(1, 1)*e12(1,1) + c13d(2, 1)*e12(1,2) + e12(1,3);
    B = c13d(1, 1)*e12(2,1) + c13d(2, 1)*e12(2,2) + e12(2,3);
    C = c13d(1, 1)*e12(3,1) + c13d(2, 1)*e12(3,2) + e12(3,3);
    
    aa = A;
    bb = B;
    cc = C;

    A = aa * invk1(1,1);
    B = bb * invk1(2,2);
    C = aa * invk1(1,3) + bb * invk1(2,3) + cc;

    distance = abs(A * c2(i, 1) + B * c2(i, 2) + C) / sqrt(A^2 + B^2);
  
    a = [a; A];
    b = [b; B];
    c = [c; C];

    x = 1:s;
    y = (-A * x - C) / B;
    plot((x + s), y, 'g', 'LineWidth', 1);
    scatter((c2(i, 1)+s),c2(i, 2), 'b', 'LineWidth', 1);
    scatter(c1(i, 1),c1(i, 2), 'b', 'LineWidth', 1);
    
end

o = 0;

for i = 1:size(cc1, 1)
    cc13d = invk1*[cc1(i, 1);cc1(i, 2);1];

    A = cc13d(1, 1)*e12(1,1) + cc13d(2, 1)*e12(1,2) + e12(1,3);
    B = cc13d(1, 1)*e12(2,1) + cc13d(2, 1)*e12(2,2) + e12(2,3);
    C = cc13d(1, 1)*e12(3,1) + cc13d(2, 1)*e12(3,2) + e12(3,3);
    
    aa = A;
    bb = B;
    cc = C;

    A = aa * invk1(1,1);
    B = bb * invk1(2,2);
    C = aa * invk1(1,3) + bb * invk1(2,3) + cc;

    distance = abs(A * cc2(i, 1) + B * cc2(i, 2) + C) / sqrt(A^2 + B^2);
  
    if distance >= 2
        o = o+1;
    end
    
end

r = o/(size(cc1, 1));
disp('Outlier Ratio (Pair_1):')
disp(r);





[i3, i4, s3, s4] = find_matches(d3, d4,"sift", 0.8);
mm = [i3, s3];
smm = sortrows(mm, 3);
ii1 = smm(:, 1:2);
ss1 = smm(:, 3);
mm = smm(1:20, :);
i1 = mm(:, 1:2);
s1 = mm(:, 3);
[c1, c2] = coordinates_find_sift(i1, s1, f3, f4);
[cc1, cc2] = coordinates_find_sift(ii1, ss1, f3, f4);

c1 = [c1, ones(20,1)];
c2 = [c2, ones(20,1)];

a = [];
b = [];
c = [];

invk1 = inv(k2);

figure();
I34 = [I3, I4];
imshow(I34);
hold on;

s = size(I34, 2)./2;

for i = 1:20
    c13d = invk1*[c1(i, 1);c1(i, 2);1];

    A = c13d(1, 1)*e34(1,1) + c13d(2, 1)*e34(1,2) + e34(1,3);
    B = c13d(1, 1)*e34(2,1) + c13d(2, 1)*e34(2,2) + e34(2,3);
    C = c13d(1, 1)*e34(3,1) + c13d(2, 1)*e34(3,2) + e34(3,3);
    
    aa = A;
    bb = B;
    cc = C;

    A = aa * invk1(1,1);
    B = bb * invk1(2,2);
    C = aa * invk1(1,3) + bb * invk1(2,3) + cc;

    distance = abs(A * c2(i, 1) + B * c2(i, 2) + C) / sqrt(A^2 + B^2);
  
    a = [a; A];
    b = [b; B];
    c = [c; C];

    x = 1:s;
    y = (-A * x - C) / B;
    plot((x + s), y, 'g', 'LineWidth', 1);
    scatter((c2(i, 1)+s),c2(i, 2), 'b', 'LineWidth', 1);
    scatter(c1(i, 1),c1(i, 2), 'b', 'LineWidth', 1);
    
end

o = 0;

for i = 1:size(cc1, 1)
    cc13d = invk1*[cc1(i, 1);cc1(i, 2);1];

    A = cc13d(1, 1)*e34(1,1) + cc13d(2, 1)*e34(1,2) + e34(1,3);
    B = cc13d(1, 1)*e34(2,1) + cc13d(2, 1)*e34(2,2) + e34(2,3);
    C = cc13d(1, 1)*e34(3,1) + cc13d(2, 1)*e34(3,2) + e34(3,3);
    
    aa = A;
    bb = B;
    cc = C;

    A = aa * invk1(1,1);
    B = bb * invk1(2,2);
    C = aa * invk1(1,3) + bb * invk1(2,3) + cc;

    distance = abs(A * cc2(i, 1) + B * cc2(i, 2) + C) / sqrt(A^2 + B^2);
  
    if distance >= 2
        o = o+1;
    end
    
end

r = o/(size(cc1, 1));
disp('Outlier Ratio (Pair_2):')
disp(r);



function [coord1,coord2] = coordinates_find_sift(best_match,simi_match,input_f1,input_f2)
    in_best_match = [best_match, simi_match];
    sort_best_match = sortrows(in_best_match,3);
    length_best_match = size(sort_best_match,1);
    coord1=zeros(length_best_match,2);
    coord2=zeros(length_best_match,2);
    
    max_sort_best_match1 = max(sort_best_match(:,1));
    max_sort_best_match2 = max(sort_best_match(:,2));

    if max_sort_best_match1 >= max_sort_best_match2
        large_match = sort_best_match(:,1);
        sec_match = sort_best_match(:,2);
    else
        large_match = sort_best_match(:,2);
        sec_match = sort_best_match(:,1);
    end
    
    num_corner1 = size(input_f1,2);
    num_corner2 = size(input_f2,2);

    if num_corner1 > num_corner2
        ind_f1_temp = large_match;
        ind_f2_temp = sec_match;
    else
        ind_f1_temp = sec_match;
        ind_f2_temp = large_match;
    end
    
    for i_match=1:length_best_match
        ind_f1=ind_f1_temp(i_match);
        ind_f2=ind_f2_temp(i_match);
        f1 = input_f1(1:2,ind_f1); 
        f2 = input_f2(1:2,ind_f2);
        coord1(i_match,:)=[f1(1), f1(2)];
        coord2(i_match,:)=[f2(1), f2(2)];
    end

end

function [indx_1st_best, indx_2nd_best, simi_1st_best, simi_2nd_best] = find_matches(D1, D2, similarity_type, alpha)

    if similarity_type == "SSD"
        
        num_D1 = size(D1,3);

        indx_1st_best = zeros(num_D1,2);
        indx_2nd_best = zeros(num_D1,2);
        simi_1st_best = zeros(num_D1,1);
        simi_2nd_best = zeros(num_D1,1);
        
        for i_D1=1:num_D1
            corr_D1 = D1(:,:,i_D1);
            diff = D2-corr_D1;
            diff = power(diff,2);
            diff = sum(diff,[1 2]);
            [first_min, first_min_ind] = min(diff,[],3);
            second_min= min(diff(diff>first_min));
            second_min_ind = find(diff==second_min,1);
            indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
            indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
            simi_1st_best(i_D1)=first_min;
            simi_2nd_best(i_D1)=second_min;

            if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                indx_1st_best(i_D1,:)=0;
                indx_2nd_best(i_D1,:)=0;
                simi_1st_best(i_D1)=0;
                simi_2nd_best(i_D1)=0;
            end
        end

        indx_1st_best(~any(simi_1st_best,2),:)=[];
        indx_2nd_best(~any(simi_2nd_best,2),:)=[];
        simi_1st_best(~any(simi_1st_best,2),:)=[];
        simi_2nd_best(~any(simi_2nd_best,2),:)=[];

        

    elseif similarity_type == "NCC"
        num_D1 = size(D1,3);

        indx_1st_best = zeros(num_D1,2);
        indx_2nd_best = zeros(num_D1,2);
        simi_1st_best = zeros(num_D1,1);
        simi_2nd_best = zeros(num_D1,1);

        for i_D1=1:num_D1
            corr_D1_temp = D1(:,:,i_D1);
            mean_corr_D1 = mean(corr_D1_temp,[1 2]);
            abs_cent = norm(corr_D1_temp);
            corr_D1 = corr_D1_temp - mean_corr_D1;
            corr_D1 = corr_D1./abs_cent;
            corr_D2_temp = D2;
            mean_corr_D2 = mean(corr_D2_temp,[1 2]);
            corr_D2 = corr_D2_temp - mean_corr_D2;
            abs_cent = pagenorm(corr_D2_temp);
            corr_D2 = corr_D2./abs_cent;
            diff = corr_D2 - corr_D1;
            diff = power(diff,2);
            diff = sum(diff,[1 2]);
            [first_min, first_min_ind] = min(diff,[],3);
            second_min = min(diff(diff>first_min));
            second_min_ind = find(diff==second_min,1);
            indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
            indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
            simi_1st_best(i_D1)=first_min;
            simi_2nd_best(i_D1)=second_min;

            if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                indx_1st_best(i_D1,:)=0;
                indx_2nd_best(i_D1,:)=0;
                simi_1st_best(i_D1)=0;
                simi_2nd_best(i_D1)=0;
            end
        end

        indx_1st_best(~any(simi_1st_best,2),:)=[];
        indx_2nd_best(~any(simi_2nd_best,2),:)=[];
        simi_1st_best(~any(simi_1st_best,2),:)=[];
        simi_2nd_best(~any(simi_2nd_best,2),:)=[];


    elseif similarity_type == "Chi−Square"
        num_D1 = size(D1,1);

        indx_1st_best = zeros(num_D1,2);
        indx_2nd_best = zeros(num_D1,2);
        simi_1st_best = zeros(num_D1,1);
        simi_2nd_best = zeros(num_D1,1);

        for i_D1=1:num_D1
            corr_D1 = D1(i_D1,:);
            diff = D2-corr_D1;
            diff = power(diff,2);
            summ = D2+corr_D1;
            out = diff./summ;
            out(find(isnan(out)))=0;
            out = sum(out, 2);
            out = 0.5*out;
            [first_min, first_min_ind] = min(out,[],1);
            [second_min, second_min_ind]= min(out(out>first_min));
            indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
            indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
            simi_1st_best(i_D1)=double(first_min);
            simi_2nd_best(i_D1)=double(second_min);

            if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                indx_1st_best(i_D1,:)=0;
                indx_2nd_best(i_D1,:)=0;
                simi_1st_best(i_D1)=0;
                simi_2nd_best(i_D1)=0;
            end
        end

        indx_1st_best(~any(simi_1st_best,2),:)=[];
        indx_2nd_best(~any(simi_2nd_best,2),:)=[];
        simi_1st_best(~any(simi_1st_best,2),:)=[];
        simi_2nd_best(~any(simi_2nd_best,2),:)=[];

    elseif similarity_type == "sift"
        num_D1 = size(D1,2);

        indx_1st_best = zeros(num_D1,2);
        indx_2nd_best = zeros(num_D1,2);
        simi_1st_best = zeros(num_D1,1);
        simi_2nd_best = zeros(num_D1,1);

        for i_D1=1:num_D1
            corr_D1 = D1(:,i_D1);
            diff = D2-corr_D1;
            diff = power(diff,2);
            diff = sum(diff,1);
            diff = sqrt(diff);

            [first_min, first_min_ind] = min(diff,[],2);
            second_min= min(diff(diff>first_min));
            second_min_ind = find(diff==second_min,1);

            indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
            indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
            simi_1st_best(i_D1)=first_min;
            simi_2nd_best(i_D1)=second_min;

            if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                indx_1st_best(i_D1,:)=0;
                indx_2nd_best(i_D1,:)=0;
                simi_1st_best(i_D1)=0;
                simi_2nd_best(i_D1)=0;
            end
        end
    
        indx_1st_best(~any(simi_1st_best,2),:)=[];
        indx_2nd_best(~any(simi_2nd_best,2),:)=[];
        simi_1st_best(~any(simi_1st_best,2),:)=[];
        simi_2nd_best(~any(simi_2nd_best,2),:)=[];

    else
        error("Unknown input. Please type SSD, NCC, Chi-Square or SIFT.");
    end

end