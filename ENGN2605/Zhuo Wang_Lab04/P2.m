close all;
clc;
clear;

imPath = 'Images\problem_2\car'; imExt = 'jpg';

if isdir(imPath) == 0
    error('User error : The image directory does not exist');
end

filearray = dir([imPath filesep '*.' imExt]); 
NumImages = size(filearray,1);
if NumImages < 0
    error('No image in the directory');
end

disp('Loading image files from the video sequence, please be patient...');

imgname = [imPath filesep filearray(1).name];
I = imread(imgname);
if size(I, 3)==3
   I = rgb2gray(I);
end
VIDEO_WIDTH = size(I,2);
VIDEO_HEIGHT = size(I,1);

ImSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
for i=1:NumImages
    imgname = [imPath filesep filearray(i).name];
    img = imread(imgname);
    img = rgb2gray(img);
    ImSeq(:,:,i) = img;
end
disp(' ... OK!');

scc = strong_corner_detector(ImSeq(:,:,1));

[row, col] = find(scc == 1);

prow = row;
pcol = col;

for n = 1:NumImages-1
    Frame = ImSeq(:,:,n);
    [dx, dy] = gradient(Frame);

    figure(1);
    imshow(Frame/256);
    hold on;
    plot(pcol, prow,'ro','MarkerSize', 3.5, 'MarkerFaceColor','auto');

    for i = 1:length(row)
        px = prow(i);
        py = pcol(i);
    
        m = 5;
        [Y, X] = meshgrid(px - m:px + m, py - m:py + m);
    
        
        extrapval = 0;

        W0 = interp2(Frame, X, Y, 'linear', extrapval);
    
        Ix = interp2(dx, X, Y, 'linear', extrapval);
        Iy = interp2(dy, X, Y, 'linear', extrapval);
    
        M = [sum(Ix(:).^2), sum(Ix(:).*Iy(:)); sum(Ix(:).*Iy(:)), sum(Iy(:).^2)];

        v0 = [0; 0];
        
        K = 15;

        NextFrame = ImSeq(:,:,(n+1));

        for k=1:K
            vk = v0;

            X = X + vk(1,1);
            Y = Y + vk(2,1);

            wki = interp2(NextFrame, X, Y, 'linear', extrapval);

            oI = W0 - wki;

            b_bar = [sum(oI(:).*Ix(:)); sum(oI(:).*Iy(:))];

            nk = (inv(M))*b_bar;

            vk = vk + nk;
            
            if norm(nk) < 0.01
                break;
            end

        end
        
        px = px + vk(1,1);
        py = py + vk(2,1);
       
        prow(i,1) = px;
        pcol(i,1) = py;

    end
    
    
end

function corners = strong_corner_detector(image)  

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
        
    [rows, cols] = size(R);

    sc = zeros(rows,cols);

    [row, col] = find(corner_markers == 1);

    sc(row, col) = R(row, col);

    fsc = sc(:);

    rfsc = sort(fsc, 'descend');

    tw = rfsc(20);

    for i=1:rows
        for j=1:cols
            if sc(i, j)>=tw
                sc(i, j)=1;
            else
                sc(i,j)=0;
            end
        end
    end
    
    corners = sc;

end

