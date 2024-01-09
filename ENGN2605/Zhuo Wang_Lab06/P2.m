clc;
clear;
close all;
format long;

%Run the 'vl_setup' file before running it

I1 = single(rgb2gray(imread("images\problem_2\half_dome\halfdome-05.png")));
I2 = single(rgb2gray(imread("images\problem_2\half_dome\halfdome-06.png")));
I3 = single(rgb2gray(imread("images\problem_2\half_dome\halfdome-07.png")));
I4 = single(rgb2gray(imread("images\problem_2\hotel\hotel-00.png")));
I5 = single(rgb2gray(imread("images\problem_2\hotel\hotel-01.png")));
I6 = single(rgb2gray(imread("images\problem_2\hotel\hotel-02.png")));

im1 = imread("images\problem_2\half_dome\halfdome-05.png");
im2 = imread("images\problem_2\half_dome\halfdome-06.png");
im3 = imread("images\problem_2\half_dome\halfdome-07.png");
im4 = imread("images\problem_2\hotel\hotel-00.png");
im5 = imread("images\problem_2\hotel\hotel-01.png");
im6 = imread("images\problem_2\hotel\hotel-02.png");

peakthresh = 0.1;
edgethresh = 5;

[f1, d1] = vl_sift(I1,'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f2, d2] = vl_sift(I2,'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f3, d3] = vl_sift(I3,'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f4, d4] = vl_sift(I4,'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f5, d5] = vl_sift(I5,'PeakThresh',peakthresh,'edgethresh',edgethresh);
[f6, d6] = vl_sift(I6,'PeakThresh',peakthresh,'edgethresh',edgethresh);

d1 = double(d1);
d2 = double(d2);
d3 = double(d3);
d4 = double(d4);
d5 = double(d5);
d6 = double(d6);

[i1, i2, s1, s2] = find_matches(d1, d2,"sift", 0.65);
[c1, c2] = coordinates_find_sift(i1, s1, f1, f2);
h1 = Ransac4Homography(c1, c2);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h1, im1, im2);
blendType = 'weightBlend';
[stitchedImage] = blendImgs(warpedImage1, im1, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage));

stitchedImages = single(stitchedImage(:,:,1));

[f7, d7] = vl_sift(stitchedImages,'PeakThresh',peakthresh,'edgethresh',edgethresh);

d7 = double(d7);

[i1, i2, s1, s2] = find_matches(d7, d3,"sift", 0.65);
[c3, c4] = coordinates_find_sift(i1, s1, f7, f3);
h2 = Ransac4Homography(c3, c4);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h2, stitchedImage, im3);
blendType = 'weightBlend';
[stitchedImage_1] = blendImgs(warpedImage1, stitchedImage, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage_1));


[i1, i2, s1, s2] = find_matches(d4, d5,"sift", 0.65);
[c5, c6] = coordinates_find_sift(i1, s1, f4, f5);
h3 = Ransac4Homography(c5, c6);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h3, im4, im5);
blendType = 'weightBlend';
[stitchedImage] = blendImgs(warpedImage1, im4, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage));

stitchedImages = single(stitchedImage(:,:,1));

[f8, d8] = vl_sift(stitchedImages,'PeakThresh',peakthresh,'edgethresh',edgethresh);

d8 = double(d8);

[i1, i2, s1, s2] = find_matches(d8, d6,"sift", 0.65);
[c7, c8] = coordinates_find_sift(i1, s1, f8, f6);
h4 = Ransac4Homography(c7, c8);

[warpedImage1, leftTopUnwarpX1, leftTopUnwarpY1, warpImgWeight1] = getNewImg(h4, stitchedImage, im6);
blendType = 'weightBlend';
[stitchedImage_2] = blendImgs(warpedImage1, stitchedImage, leftTopUnwarpX1, leftTopUnwarpY1, blendType, warpImgWeight1);
figure();
imshow(uint8(stitchedImage_2));


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