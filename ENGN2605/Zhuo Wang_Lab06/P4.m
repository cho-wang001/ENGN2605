clc;
clear;
close all;
format long;

%Run the 'vl_setup' file before running it

I1 = double(im2gray(imread("images\problem_1\Golden_Gate\goldengate-02.png")));
I2 = double(im2gray(imread("images\problem_1\Golden_Gate\goldengate-03.png")));
I3 = double(rgb2gray(imread("images\problem_1\LEMS_lab\lems-01.png")));
I4 = double(rgb2gray(imread("images\problem_1\LEMS_lab\lems-02.png")));

im1 = imread("images\problem_1\Golden_Gate\goldengate-02.png");
im2 = imread("images\problem_1\Golden_Gate\goldengate-03.png");
im3 = imread("images\problem_1\LEMS_lab\lems-01.png");
im4 = imread("images\problem_1\LEMS_lab\lems-02.png");


peakthresh = 0.1;
edgethresh = 5;

[f1, d1] = vl_sift(single(I1),'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f2, d2] = vl_sift(single(I2),'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f3, d3] = vl_sift(single(I3),'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f4, d4] = vl_sift(single(I4),'PeakThresh',peakthresh,'edgethresh',edgethresh);


d1 = double(d1);
d2 = double(d2);
d3 = double(d3);
d4 = double(d4);

cd1 = corner_detector(I1);
cd2 = corner_detector(I2);
cd3 = corner_detector(I3);
cd4 = corner_detector(I4);

r1 = region_descriptors(I1, cd1, 21, 'pixels');
r2 = region_descriptors(I2, cd2, 21, 'pixels');
r3 = region_descriptors(I3, cd3, 21, 'pixels');
r4 = region_descriptors(I4, cd4, 21, 'pixels');



[i1, i2, s1, s2] = find_matches(r1, r2,"SSD", 0.7);
[c1, c2] = coordinates_find(i1, s1, cd1, cd2);
disp('SSD GATE ↓');
h1 = Ransac4Homography(c1, c2);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h1, im1, im2);
blendType = 'weightBlend';
[stitchedImage] = blendImgs(warpedImage1, im1, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage));



[i1, i2, s1, s2] = find_matches(d1, d2,"sift", 0.65);
[c1, c2] = coordinates_find_sift(i1, s1, f1, f2);
disp('SIFT GATE ↓');
h2 = Ransac4Homography(c1, c2);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h2, im1, im2);
blendType = 'weightBlend';
[stitchedImage] = blendImgs(warpedImage1, im1, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage));




[i1, i2, s1, s2] = find_matches(r3, r4,"NCC", 0.7);
[c3, c4] = coordinates_find(i1, s1, cd3, cd4);
disp('NCC LAB ↓');
h3 = Ransac4Homography(c3, c4);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h3, im3, im4);
blendType = 'weightBlend';
[stitchedImage] = blendImgs(warpedImage1, im3, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage));



[i1, i2, s1, s2] = find_matches(d3, d4,"sift", 0.65);
[c3, c4] = coordinates_find_sift(i1, s1, f3, f4);
disp('SIFT LAB ↓');
h4 = Ransac4Homography(c3, c4);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h4, im3, im4);
blendType = 'weightBlend';
[stitchedImage] = blendImgs(warpedImage1, im3, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage));




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

function finalH = Ransac4Homography(matches1, matches2)
    I = 1000;

    H = [];

    ml = 0;

    for i = 1:I
        nm = size(matches1,1);
        ns = 4;
        sc = randsample(nm,ns,'true');

        for i=sc
            x_coord1 = matches1(i,1);
            y_coord1 = matches1(i,2);
            x_coord2 = matches2(i,1);
            y_coord2 = matches2(i,2);
        end
        
        a = [];

        nc = size(x_coord1, 1);

        for i_ind=1:nc
               
            var1 = (-1)*x_coord1(i_ind);
            var2 = (-1)*y_coord1(i_ind);
            var3 = (-1)*1;
            var4 = 0;
            var5 = 0;
            var6 = 0;
            var7 = x_coord1(i_ind)*x_coord2(i_ind);
            var8 = x_coord2(i_ind)*y_coord1(i_ind);
            var9 = x_coord2(i_ind);
            a_vec1 = [var1 var2 var3 var4 var5 var6 var7 var8 var9];
            a = cat(1,a,a_vec1);
            
            var1 = 0;
            var2 = 0;
            var3 = 0;
            var4 = (-1)*x_coord1(i_ind);
            var5 = (-1)*y_coord1(i_ind);
            var6 = (-1)*1;
            var7 = x_coord1(i_ind)*y_coord2(i_ind);
            var8 = y_coord1(i_ind)*y_coord2(i_ind);
            var9 = y_coord2(i_ind);
            a_vec2 = [var1 var2 var3 var4 var5 var6 var7 var8 var9];
            a = cat(1,a,a_vec2);
                
        end
        
        [u, s, v] = svd(a);
        h = v(:, end);
        h = reshape(h,[3 3]);

        h = h';

        im = size(matches1, 1);
        om = zeros(im, 1);

        for j=1:im
            hom_vec = [matches1(j,1); matches1(j,2); 1];
            out_vec = h*hom_vec;
            om(j,1) = out_vec(1)./out_vec(3);
            om(j,2) = out_vec(2)./out_vec(3);
        end

        l = 0;

       for n=1:im
           out = om(n,:)-matches2(n,:);
           out = power(out(1),2)+power(out(2),2);
           out = sqrt(out);
           
           if out <= 1
              l = l + 1;
           end
       end
   
        
       if l >= ml
          ml = l;
          H = h;
       end

    end
    

    disp("max_inlier:")
    disp(ml)
    finalH = H;

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
        error("Unknown input. Please type SSD, NCC, or Chi-Square.");
    end

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