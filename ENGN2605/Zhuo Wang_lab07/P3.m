clc;
clear;
close all;
format long;


%Run the 'vl_setup' file before running it
I1 = rgb2gray(imread('Q3_data\1.jpg'));
I2 = rgb2gray(imread('Q3_data\2.jpg'));

rg = load("Q3_data\GT_pose.mat").GT_pose.R;
tg = load("Q3_data\GT_pose.mat").GT_pose.T;

k = load("Q3_data\Calibration_Matrix.mat").K;

i1 = single(I1);
i2 = single(I2);

peak_thresh = 2.5;
edgethresh = 5;

[f1, d1] = vl_sift(i1,'PeakThresh',peak_thresh,'EdgeThresh',edgethresh);
[f2, d2] = vl_sift(i2,'PeakThresh',peak_thresh,'EdgeThresh',edgethresh);

[i1, i2, s1, s2] = find_matches(d1, d2,"sift", 0.65);

[c1, c2] = coordinates_find_sift(i1, s1, f1, f2);

ik = inv(k);

EE = Ransac4Essential(c1, c2, ik);

[U, S, V] = svd(EE);
W = [0, -1, 0; 1, 0, 0; 0, 0, 1];

R1 = U * W * V';
T1 = U(:, 3);
R2 = U * W' * V';
T2 = -U(:, 3);

if det(R1) < 0
    R1 = -R1;
    T1 = -T1;
end

if det(R2) < 0
    R2 = -R2;
    T2 = -T2;
end

rotation_error1 = real(acos(0.5 * (trace(rg * R1' - eye(3)))));
translation_error1 = abs(dot(tg, T1) - 1);
disp("rotation_error1:");
disp(rotation_error1);
disp("translation_error1:");
disp(translation_error1);

rotation_error2 = real(acos(0.5 * (trace(rg * R2' - eye(3)))));
translation_error2 = abs(dot(tg, T2) - 1);
disp("rotation_error2:");
disp(rotation_error2);
disp("translation_error2:");
disp(translation_error2);



function finalE = Ransac4Essential(gamma1, gamma2, K)
    I = 1000;

    E = [];

    ml = 0;

    for i = 1:I
        c1m = [];
        c2m = [];
        nm = size(gamma1, 1);
        sc = randsample(nm, 5,'true');

        for i=sc
            x_coord1 = gamma1(i,1);
            y_coord1 = gamma1(i,2);
            x_coord2 = gamma2(i,1);
            y_coord2 = gamma2(i,2);
        end

        cc1 = [x_coord1, y_coord1];
        cc2 = [x_coord2, y_coord2];
        
        for i = 1:size(cc1, 1)
            g1 = K*[cc1(i, 1);cc1(i, 2);1];
            g1 = g1';
            c1m = [c1m; g1];
        
            g2 = K*[cc2(i, 1);cc2(i, 2);1];
            g2 = g2';
            c2m = [c2m; g2];
        end
        
        mi = cat(3, c1m, c2m);
        
        es = fivePointAlgorithmSelf(mi);
        
        for j = 1:size(es, 3)
            e = es{:,:,j};
            
            l = 0;

            for n = 1:size(gamma1, 1)
                c13d = K*[gamma1(n, 1); gamma1(n, 2); 1];

                A = c13d(1, 1)*e(1,1) + c13d(2, 1)*e(1,2) + e(1,3);
                B = c13d(1, 1)*e(2,1) + c13d(2, 1)*e(2,2) + e(2,3);
                C = c13d(1, 1)*e(3,1) + c13d(2, 1)*e(3,2) + e(3,3);
    
                aa = A;
                bb = B;
                cc = C;

                A = aa * K(1,1);
                B = bb * K(2,2);
                C = aa * K(1,3) + bb * K(2,3) + cc;

                distance = abs(A * gamma2(n, 1) + B * gamma2(n, 2) + C) / sqrt(A^2 + B^2);

                if distance<2
                    l = l+1;
                end

            end

            if l >= ml
                ml = l;
                E = e;
            end

        end
        
    end
    

    disp("max_inlier:")
    disp(ml)
    finalE = E;

end

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