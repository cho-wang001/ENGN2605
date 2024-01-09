close all;
clc;
clear;

I1 = imread('bds1.jpg');
I1e = text_edge_detector(I1,5);

function [edge_img]=text_edge_detector(image,threshold)
    
    if size(image, 3) == 3
        image = double(rgb2gray(image));
    else
        image = double(image);
    end
    
    range = -24:1:24;

    [X, Y] = meshgrid(range, range);
    
    k_values = [1, 2, 3];
    scales_sigma1 = 2.^(1./k_values);
    scales_sigma2 = 3.*(2.^(1./k_values));
    
    filter_bank = cell(1, 46);

    orientations = linspace(0, -pi, 6);
    
    
    filter_bank_index_18_o = 1;
  
    for i = 1:length(k_values)
        for j = 1:6
            sigma1 = scales_sigma1(i);
            sigma2 = scales_sigma2(i);
            
            f18o = @(x, y) ((-1/(sigma1^3*sqrt(2*pi)))*y.*exp(-y.^2/(2*(sigma1^2)))) .*  ((1/(sigma2*sqrt(2*pi))).*exp(-x.^2/(2*(sigma2^2))));

            rf18o = f18o(X,Y);
            rf18ob = imrotate(rf18o, rad2deg(orientations(j)),'crop');

            filter_bank{filter_bank_index_18_o} = rf18ob;

            filter_bank_index_18_o = filter_bank_index_18_o + 1;
        end
    end

    filter_bank_index_18_e = 1;

    for i = 1:length(k_values)
        for j = 1:6
            sigma1 = scales_sigma1(i);
            sigma2 = scales_sigma2(i);
            
            f18e = @(x, y) ((1/(sigma1^5*sqrt(2*pi)))*(y.^2).*exp(-y.^2/(2*(sigma1^2)))) .*  ((1/(sigma2*sqrt(2*pi))).*exp(-x.^2/(2*(sigma2^2))));

            rf18e = f18e(X,Y);
            rf18eb = imrotate(rf18e, rad2deg(orientations(j)),'crop');

            filter_bank{18 + filter_bank_index_18_e} = rf18eb;

            filter_bank_index_18_e = filter_bank_index_18_e + 1;
        end
    end
    
    lk_values = [1, 2, 3, 4];
    lscales_sigma1 = 2.^(1./lk_values);
    lscales_sigma2 = 3.*(2.^(1./lk_values));

    filter_bank_index_L_1 = 1;
    filter_bank_index_L_2 = 2;

    for i = 1:length(lk_values)
        lsigma1 = lscales_sigma1(i);
        l = fspecial('log', 24, lsigma1);
        filter_bank{36 + filter_bank_index_L_1} = l;
        filter_bank_index_L_1 = filter_bank_index_L_1 + 2;


        lsigma2 = lscales_sigma2(i);
        l2 = fspecial('log', 24, lsigma2);
        filter_bank{36 + filter_bank_index_L_2} = l2;
        filter_bank_index_L_2 = filter_bank_index_L_2 + 2;
    end
    
    gk_values = [1, 2, 3, 4];
    gscales_sigma1 = 2.^(1./gk_values);

    filter_bank_index_g = 1;
    
    for i = 1:length(gk_values)
        gsigma = gscales_sigma1(i);
        g = fspecial('gaussian', 24, gsigma);
        filter_bank{44 + filter_bank_index_g} = g;
        filter_bank_index_g = filter_bank_index_g + 1;
    end

    num_filters = numel(filter_bank);

    num_rows = ceil(sqrt(num_filters));
    num_cols = ceil(num_filters / num_rows);

    figure();
    
    for i = 1:num_filters
        subplot(num_rows, num_cols, i);
        filter = filter_bank{i};
        imagesc(range, range, filter);
        axis off; 
        title(['Filter ' num2str(i)]);
       
    end
    
    if size(image, 3) == 3
        image = im2double(rgb2gray(image));
    else
        image = im2double(image);
    end
    
    [height, width] = size(image);

    num_filters = numel(filter_bank);
    filter_responses = zeros(height, width, num_filters);

    for i = 1:num_filters
        filter = filter_bank{i};
        filtered_response = imfilter(image, filter, 'conv', 'symmetric');
        filter_responses(:, :, i) = filtered_response;
    end
    
    edge_img = filtered_response;

end
