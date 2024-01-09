close all;
clc;
clear;


I1 = imread("bds1.jpg");

i1316 = hist_edge_detector(I1, 3, 6, 16);
i1516 = hist_edge_detector(I1, 5, 6, 16);
i1332 = hist_edge_detector(I1, 3, 6, 32);




function [edge_map,orient_map]=hist_edge_detector(img,rad,num_orient,num_bins)
    if size(img, 3) == 3
        img = im2double(rgb2gray(img));
    else
        img = im2double(img);
    end

    [rows, cols] = size(img);
    
    orient_map = zeros(rows, cols);
    edge_map = zeros(rows, cols);

    [xm, ym] = meshgrid(1:cols, 1:rows);
    
    
    for x = rad+1:rows-rad-1
        for y = rad+1:cols-rad-1
            
            on = zeros(1, num_orient);
            hist_left = zeros(1, num_bins);
            hist_right = zeros(1, num_bins);

            onn = 1;

            distance_to_center = sqrt((x - ym).^2 + (y - xm).^2);

            linear_indices = find(distance_to_center <= rad);

            [mr, mc] = ind2sub(size(distance_to_center), linear_indices);

            for angle = linspace(0, pi, num_orient)
                hist_left(:) = 0;
                hist_right(:) = 0;
                for i = 1:size(mr)
                    a = atan2(mc(i) - y, mr(i) - x);
                    if a >= angle || a <= (angle - pi)
                        left_pixel = img(mr(i), mc(i));
                        bin_left = floor((left_pixel) * (num_bins-1)) + 1;
                        hist_left(bin_left) = hist_left(bin_left) + 1;
                    else
                        right_pixel = img(mr(i), mc(i));
                        bin_right = floor((right_pixel) * (num_bins-1)) + 1;
                        hist_right(bin_right) = hist_right(bin_right) + 1;
                    end
                end
                
                min_val_l = min(hist_left);
                max_val_l = max(hist_left);
                normalized_hist_l = (hist_left - min_val_l) / (max_val_l - min_val_l);

                min_val_r = min(hist_right);
                max_val_r = max(hist_right);
                normalized_hist_r = (hist_right - min_val_r) / (max_val_r - min_val_r);

                chi_squared_distance = 0.5.*sum((normalized_hist_l - normalized_hist_r).^2 ./ (normalized_hist_l + normalized_hist_r + eps));

                on(onn) = chi_squared_distance;

                onn = onn +1;
            end
            
            [~, index] = max(on);

            array = linspace(0, pi, num_orient);
            value = array(index);
            value = value/pi;
            orient_map(x,y)=value;
        end

        directions = [0, 1/4, 1/2, 3/4, 1];

        edge_map=orient_map;

        [rows, cols] = size(edge_map);

        for i = 2:rows-1
            for j = 2:cols-1

                [~, index] = min(abs(edge_map(i, j) - directions));
        
                if index == 1 || index == 5
                    q = edge_map(i, j+1);
                    q_bar = edge_map(i, j-1);
                elseif index == 2
                    q = edge_map(i-1, j+1);
                    q_bar = edge_map(i+1, j-1);
                elseif index == 3
                    q = edge_map(i-1, j);
                    q_bar = edge_map(i+1, j);
                else
                    q = edge_map(i-1, j-1);
                    q_bar = edge_map(i+1, j+1);
                end
        
                if edge_map(i, j) < max(q, q_bar)
                    edge_map(i, j) = 0;
                end
            end
        end

    end

    figure();
    imshow(edge_map);


end
