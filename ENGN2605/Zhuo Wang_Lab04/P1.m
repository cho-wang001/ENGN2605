close all;
clc;
clear;

P1I1 = imread("Images\problem_1\board.png");
P1I2 = imread("Images\problem_1\coliseum.jpg");
P1I3 = imread("Images\problem_1\gantrycrane.png");

fP1I1 = corner_detector(P1I1);
fP1I2 = corner_detector(P1I2);
fP1I3 = corner_detector(P1I3);

[row, col] = find(fP1I1 == 1);

figure();
imshow(rgb2gray(P1I1));
hold on;
plot(col, row,'o','MarkerSize', 3.5, 'MarkerFaceColor','auto');

[row, col] = find(fP1I2 == 1);

figure();
imshow(rgb2gray(P1I2));
hold on;
plot(col, row,'o','MarkerSize', 3.5, 'MarkerFaceColor','auto');

[row, col] = find(fP1I3 == 1);

figure();
imshow(rgb2gray(P1I3));
hold on;
plot(col, row,'o','MarkerSize', 3.5, 'MarkerFaceColor','auto');

I = checkerboard(50, 10, 10);
corners = corner_detector(I);

[row, col] = find(corners == 1);

figure;
imshow(I);
hold on;
plot(col, row, 'o', 'MarkerSize', 5, 'MarkerFaceColor','auto');

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

