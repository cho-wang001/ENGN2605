clc;
clear;
close all;
format long;


image_dir1 = "images\illumination_change\";
image_dir2 = "images\planar_rotation_change\";
image_dir3 = "images\scale_rotation_change\";
image_dir4 = "images\view_change\";

I1 = rgb2gray(imread("images\illumination_change\img1.png"));
I2 = rgb2gray(imread("images\illumination_change\img2.png"));
I3 = rgb2gray(imread("images\planar_rotation_change\img1.png"));
I4 = rgb2gray(imread("images\planar_rotation_change\img2.png"));
I5 = im2gray(imread("images\scale_rotation_change\img1.png"));
I6 = im2gray(imread("images\scale_rotation_change\img2.png"));
I7 = rgb2gray(imread("images\view_change\img1.png"));
I8 = rgb2gray(imread("images\view_change\img2.png"));

j1 = corner_detector(I1);
j2 = corner_detector(I2);
j3 = corner_detector(I3);
j4 = corner_detector(I4);
j5 = corner_detector(I5);
j6 = corner_detector(I6);
j7 = corner_detector(I7);
j8 = corner_detector(I8);


R1 = region_descriptors(I1, j1, 21, 'pixels');
R2 = region_descriptors(I2, j2, 21, 'pixels');

[i1, i2, s1, s2] = find_matches_P2(R1, R2,"SSD", 0.7);

[c1, c2] = coordinates_find(i1, s1, j1, j2);

visualize_matches(I1, I2, c1, c2, 20);

a1 = accuracy_compute(image_dir1, c1, c2);
disp("SSD accuracy rate is:");
disp(a1);


R3 = region_descriptors(I3, j3, 21, 'pixels');
R4 = region_descriptors(I4, j4, 21, 'pixels');

[i1, i2, s1, s2] = find_matches_P2(R3, R4,"NCC", 0.7);

[c3, c4] = coordinates_find(i1, s1, j3, j4);

visualize_matches(I3, I4, c3, c4, 20);

a2 = accuracy_compute(image_dir2, c3, c4);
disp("NCC accuracy rate is:");
disp(a2);


R5 = region_descriptors(I5, j5, 21, 'pixels');
R6 = region_descriptors(I6, j6, 21, 'pixels');

[i1, i2, s1, s2] = find_matches_P2(R5, R6,"NCC", 0.9);

[c5, c6] = coordinates_find(i1, s1, j5, j6);

visualize_matches(I5, I6, c5, c6, 20);

a3 = accuracy_compute(image_dir3, c5, c6);
disp("NCC accuracy rate is:");
disp(a3);


R7 = region_descriptors(I7, j7, 21, 'histogram');
R8 = region_descriptors(I8, j8, 21, 'histogram');

[i1, i2, s1, s2] = find_matches_P2(R7, R8, "Chi-Square", 0.75);

[c7, c8] = coordinates_find(i1, s1, j7, j8);

visualize_matches(I7, I8, c7, c8, 20);

a4 = accuracy_compute(image_dir4, c7, c8);
disp("Chi−Square accuracy rate is:");
disp(a4);



function accuracy_rate = accuracy_compute(image_dir,best_match1, best_match2)
    mat_dir = strcat(image_dir, "H_1to2.txt");
    hom_mat = readmatrix(mat_dir);
    length_corners = size(best_match1,1);
    success_time = 0;
    for i_len_corners=1:length_corners
        x_curr = best_match1(i_len_corners,1);
        y_curr = best_match1(i_len_corners,2);
        trans_vec = [x_curr; y_curr;1];
        out = hom_mat*trans_vec;
        x_proj = out(1);
        y_proj = out(2);
        scale = out(3);
        x_proj = x_proj./scale;
        y_proj = y_proj./scale;
        x_corr = best_match2(i_len_corners,1);
        y_corr = best_match2(i_len_corners,2);
        var_x = power((x_corr-x_proj),2);
        var_y = power((y_corr-y_proj),2);
        magnitude = sqrt(var_x+var_y);
        if magnitude <= 2
            success_time = success_time + 1;
        end
    end
    
    accuracy_rate = success_time/length_corners;

end


function [coord1,coord2] = coordinates_find(best_match,simi_match,input_corners1,input_corners2)
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

    num_corner1 = size(input_corners1,1);
    num_corner2 = size(input_corners2,1);
    if num_corner1 > num_corner2
        ind_corn1_temp = large_match;
        ind_corn2_temp = sec_match;
    else
        ind_corn1_temp = sec_match;
        ind_corn2_temp = large_match;
    end


    for i_match=1:length_best_match
        ind_corn1=ind_corn1_temp(i_match);
        ind_corn2=ind_corn2_temp(i_match);
        corners1 = input_corners1(ind_corn1,1:2); 
        corners2 = input_corners2(ind_corn2,1:2);
        coord1(i_match,:)=[corners1(2), corners1(1)];
        coord2(i_match,:)=[corners2(2), corners2(1)];
    end

end


function [indx_1st_best, indx_2nd_best, simi_1st_best, simi_2nd_best] = find_matches_P2(D1, D2, similarity_type, alpha)
    
    if similarity_type == "Chi-Square"
        num_D1 = size(D1,1);
        num_D2 = size(D2,1);

        if num_D1 <= num_D2
            D1_temp = D1;
            D2_temp = D2;
            num_D1 = num_D1;
        else
            D1_temp = D2;
            D2_temp = D1;
            num_D1 = num_D2;
        end
    else
        num_D1 = size(D1,3);
        num_D2 = size(D2,3);

        if num_D1 <= num_D2
            D1_temp = D1;
            D2_temp = D2;
            num_D1 = num_D1;
        else
            D1_temp = D2;
            D2_temp = D1;
            num_D1 = num_D2;
        end
    end    
    
    indx_1st_best = zeros(num_D1,2);
    indx_2nd_best = zeros(num_D1,2);
    simi_1st_best = zeros(num_D1,1);
    simi_2nd_best = zeros(num_D1,1);
    
    for i_D1=1:num_D1
        if similarity_type == "SSD"
            corr_D1 = D1_temp(:,:,i_D1);
            diff = D2_temp-corr_D1;
            diff = power(diff,2);
            diff = sum(diff,[1 2]);
            [first_min, first_min_ind] = min(diff,[],3);
            second_min= min(diff(diff>first_min));
            second_min_ind = find(diff==second_min,1);
            corr_D2 = D2_temp(:,:,first_min_ind);
            diff = D1_temp-corr_D2;
            diff = power(diff,2);
            diff = sum(diff,[1 2]);
            [first_min_back, first_min_ind_back] = min(diff,[],3);
            back_D1 = D1_temp(:,:,first_min_ind_back);

            if first_min_ind_back == i_D1
                if ~ismember(first_min_ind,indx_1st_best)
                    indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
                    indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
                    simi_1st_best(i_D1)=first_min;
                    simi_2nd_best(i_D1)=second_min;
                else
                    [row_pre_min,col_pre_min] = find(indx_1st_best(:,2)==first_min_ind);
                    pre_min_score = simi_1st_best(row_pre_min);
                    if first_min > pre_min_score
                       indx_1st_best(i_D1,:)=0;
                       indx_2nd_best(i_D1,:)=0;
                       simi_1st_best(i_D1)=0;
                       simi_2nd_best(i_D1)=0;
                    else
                       indx_1st_best(row_pre_min,:)=0;
                       indx_2nd_best(row_pre_min,:)=0;
                       simi_1st_best(row_pre_min)=0;
                       simi_2nd_best(row_pre_min)=0;
                       indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
                       indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
                       simi_1st_best(i_D1)=first_min;
                       simi_2nd_best(i_D1)=second_min;
                    end
                end
                if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                    indx_1st_best(i_D1,:)=0;
                    indx_2nd_best(i_D1,:)=0;
                    simi_1st_best(i_D1)=0;
                    simi_2nd_best(i_D1)=0;
                end
            end
    
        elseif similarity_type == "NCC"
            corr_D1_temp = D1_temp(:,:,i_D1);
            mean_corr_D1 = mean(corr_D1_temp,[1 2]);
            abs_cent = norm(corr_D1_temp);
            corr_D1 = corr_D1_temp - mean_corr_D1;
            corr_D1 = corr_D1./abs_cent;
            corr_D2_temp = D2_temp;
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
            corr_D2_bac = D2_temp(:,:,first_min_ind);
            mean_corr_D2_bac = mean(corr_D2_bac, [1 2]);
            abs_cent = norm(corr_D2_bac);
            corr_D2_bac = corr_D2_bac - mean_corr_D2_bac;
            corr_D2_bac = corr_D2_bac./abs_cent;
            corr_D1_bac = D1_temp;
            mean_corr_D1_bac = mean(corr_D1_bac, [1 2]);
            abs_cent = pagenorm(corr_D1_bac);
            corr_D1_bac = corr_D1_bac - mean_corr_D1_bac;
            corr_D1_bac = corr_D1_bac./abs_cent;
            diff = corr_D1_bac - corr_D2_bac;
            diff = power(diff,2);
            diff = sum(diff, [1 2]);
            [first_min_bac, first_min_ind_bac] = min(diff,[],3);
            back_D1 = D1_temp(:,:,first_min_ind_bac);
            if back_D1 == corr_D1_temp
                if ~ismember(first_min_ind,indx_1st_best)
                    indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
                    indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
                    simi_1st_best(i_D1)=first_min;
                    simi_2nd_best(i_D1)=second_min;
                else
                    [row_pre_min,col_pre_min] = find(indx_1st_best(:,2)==first_min_ind);
                    pre_min_score = simi_1st_best(row_pre_min);

                    if first_min >= pre_min_score
                       indx_1st_best(i_D1,:)=0;
                       indx_2nd_best(i_D1,:)=0;
                       simi_1st_best(i_D1)=0;
                       simi_2nd_best(i_D1)=0;
                    else
                       indx_1st_best(row_pre_min,:)=0;
                       indx_2nd_best(row_pre_min,:)=0;
                       simi_1st_best(row_pre_min)=0;
                       simi_2nd_best(row_pre_min)=0;
                       indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
                       indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
                       simi_1st_best(i_D1)=first_min;
                       simi_2nd_best(i_D1)=second_min;
                    end
                end

                if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                    indx_1st_best(i_D1,:)=0;
                    indx_2nd_best(i_D1,:)=0;
                    simi_1st_best(i_D1)=0;
                    simi_2nd_best(i_D1)=0;
                end

            end
    
        elseif similarity_type == "Chi-Square"
            corr_D1 = D1_temp(i_D1,:);
            diff = D2_temp-corr_D1;
            diff = power(diff,2);
            summ = D2_temp+corr_D1;
            out = diff./summ;
            out(find(isnan(out)))=0;
            out = sum(out, 2);
            out = 0.5*out;
            [first_min, first_min_ind] = min(out,[],1);
            [second_min, second_min_ind]= min(out(out>first_min));
            D2_bac = D2_temp(first_min_ind,:);
            diff = D1_temp-D2_bac;
            diff = power(diff,2);
            summ = D1_temp+D2_bac;
            out = diff./summ;
            out(find(isnan(out)))=0;
            out = sum(out,2);
            out = 0.5*out;
            [first_min_bac, first_min_ind_bac] = min(out,[],1);
            D1_bac = D1_temp(first_min_ind_bac,:);
            if D1_bac == corr_D1
               if ~ismember(first_min_ind,indx_1st_best)
                   indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
                   indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
                   simi_1st_best(i_D1)=double(first_min);
                   simi_2nd_best(i_D1)=double(second_min);
               else
                    [row_pre_min,col_pre_min] = find(indx_1st_best(:,2)==first_min_ind);
                    pre_min_score = simi_1st_best(row_pre_min);

                    if first_min >= pre_min_score
                       indx_1st_best(i_D1,:)=0;
                       indx_2nd_best(i_D1,:)=0;
                       simi_1st_best(i_D1)=0;
                       simi_2nd_best(i_D1)=0;
                    else
                       indx_1st_best(row_pre_min,:)=0;
                       indx_2nd_best(row_pre_min,:)=0;
                       simi_1st_best(row_pre_min)=0;
                       simi_2nd_best(row_pre_min)=0;
                       indx_1st_best(i_D1,:)=[i_D1, first_min_ind];
                       indx_2nd_best(i_D1,:)=[i_D1, second_min_ind];
                       simi_1st_best(i_D1)=first_min;
                       simi_2nd_best(i_D1)=second_min;
                    end   
               end

               if simi_1st_best(i_D1) >= alpha*simi_2nd_best(i_D1)
                  indx_1st_best(i_D1,:)=0;
                  indx_2nd_best(i_D1,:)=0;
                  simi_1st_best(i_D1)=0;
                  simi_2nd_best(i_D1)=0;
              end
        end
    
        else
            error("Unknown input. Please type SSD, NCC, or Chi-Square.");
        end

    end
    
    indx_1st_best(~any(simi_1st_best,2),:)=[];
    indx_2nd_best(~any(simi_2nd_best,2),:)=[];
    simi_1st_best(~any(simi_1st_best,2),:)=[];
    simi_2nd_best(~any(simi_2nd_best,2),:)=[];

end


function D = region_descriptors(image,corners, window_sz, descriptor_type)
    img_length = size(image,1);
    img_width = size(image,2);

    length_corners = size(corners,1);

    if descriptor_type == "pixels" 
        D = zeros(window_sz, window_sz, length_corners);
    elseif descriptor_type  == "histogram"
        num_bins = 32;
        D = zeros(length_corners,num_bins);
    else
        error("Unknown analysis method found. Please type in: pixels, or histogram.")
    end

    corners_x = corners(:,2);
    corners_y = corners(:,1);
    for i_corner=1:length_corners
        x_coord = corners_x(i_corner);
        y_coord = corners_y(i_corner);

        offset = floor(window_sz/2);
        x_range = x_coord-offset:x_coord+offset;
        y_range = y_coord-offset:y_coord+offset;

        if x_coord-offset < 1
           x_range = 1:x_coord+offset;
        elseif x_coord+offset > img_width
           x_range = x_coord-offset:img_width;
        end

        if y_coord-offset < 1
           y_range = 1:y_coord+offset;
        elseif y_coord+offset > img_length
           y_range = y_coord-offset:img_length;
        end

        if descriptor_type == "pixels"
           out = image(y_range,x_range);
           out = double(out);
           len_out = size(out,1);
           wid_out = size(out,2); 
           for i_len=1:len_out
               for i_wid=1:wid_out
                   D(i_len,i_wid,i_corner) = out(i_len,i_wid);
               end
           end
      
        elseif descriptor_type == "histogram"
            out = double(image(y_range,x_range));
            out = out-mean(out,[1 2]);
            out = double(out)./norm(double(image(y_range,x_range)));
            length_out = size(out,1);
            width_out = size(out,2);
            out = reshape(out,[1,length_out*width_out]);
            out = hist(double(out),num_bins);
            D(i_corner,:) = out;

        else
            error("Unknown descriptor type has been input. Please try pixels or histogram");
        end
    end
end


function corners = corner_detector(image)
    if size(image,3)==3
        gray_image = rgb2gray(image);
    else
        gray_image = image;
    end
    gray_image = double(gray_image);
    
    sigma=0.7;
    wind_size = ceil(6*sigma);
    gauss_filt = fspecial('gaussian',[wind_size wind_size],sigma);
    [Gx, Gy] = gradient(gauss_filt);
    
    f_x = conv2(gray_image,Gx,'same');
    f_y = conv2(gray_image,Gy,'same');
    
    f_x_2 = f_x.*f_x;
    f_y_2 = f_y.*f_y;
    f_x_y = f_x.*f_y;
    
    sigma2 = 2.0;
    gauss_filt2 = fspecial('gaussian',[wind_size wind_size],sigma2);
    s_x = conv2(f_x_2,gauss_filt2,'same');
    s_y = conv2(f_y_2,gauss_filt2,'same');
    s_xy = conv2(f_x_y,gauss_filt2,'same');
    
    length_s = size(s_x,1);
    width_s = size(s_x,2);
    
    R = zeros(length_s,width_s);
    
    for i_len=1:length_s
        for i_wid=1:width_s
           s_x_point = s_x(i_len,i_wid);
           s_xy_point = s_xy(i_len,i_wid);
           s_y_point = s_y(i_len,i_wid);
           M = [s_x_point s_xy_point; s_xy_point s_y_point];
           alpha = 0.04;
           R_val = det(M)-alpha*power(trace(M),2);
           R(i_len,i_wid) = R_val;
        end
    end
    

    non_max_M = zeros(length_s, width_s);

    for i_len=1:length_s
        for i_wid=1:width_s
           len_rang = i_len-1:i_len+1;
           wid_rang = i_wid-1:i_wid+1;
           if i_len-1 < 1 
               len_rang = 1:i_len+1;
           end
           if i_len+1 > length_s
               len_rang = i_len-1:length_s;
           end
           if i_wid-1 < 1
               wid_rang = 1:i_wid+1;
           end
           if i_wid+1 > width_s
               wid_rang = i_wid-1:width_s;
           end
           region_of_interest = R(len_rang,wid_rang);
           
           max_regofint = max(region_of_interest,[],'all');
           
           if R(i_len,i_wid) < max_regofint
               non_max_M(i_len,i_wid)=0;
           else
               non_max_M(i_len,i_wid) = R(i_len,i_wid);
           end
        end
    end
    

    R0 = 0.01*max(non_max_M,[],'all');
    non_max_M(non_max_M<R0) = 0;
    
    
    [corner_y, corner_x, strength] = find(non_max_M);
    
    corners = [corner_y, corner_x, strength];

end