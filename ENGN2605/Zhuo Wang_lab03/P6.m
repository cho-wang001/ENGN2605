close all;
clc;
clear;

I1 = imread("bds2.jpg");

i1 = intensity_edge_detector(I1,4,0.1);

figure();
imshow(i1);


function edge_img = intensity_edge_detector(image, Th, Tl)
    if size(image, 3) == 3
        image = double(rgb2gray(image));
    else
        image = double(image);
    end

    fs31 = fspecial('gaussian', 3, 1.0);

    image = imfilter(image, fs31, 'same');

    [dy,dx] = gradient(image);

    M = sqrt(dx.^2+dy.^2);

    figure();
    imshow(M);
    title('Magnitude Map M');
    
    figure();
    quiver(dy, dx);
    
    gdm = atan2(dy, dx);
    directions = [0, pi/4, pi/2, 3*pi/4, pi, -3*pi/4, -pi/2, -pi/4];

    [rows, cols] = size(M);

    for i = 2:rows-1
        for j = 2:cols-1

            [~, index] = min(abs(gdm(i, j) - directions));
        
            if index == 1 || index == 5
                q = M(i, j+1);
                q_bar = M(i, j-1);
            elseif index == 2 || index == 6
                q = M(i-1, j+1);
                q_bar = M(i+1, j-1);
            elseif index == 3 || index == 7
                q = M(i-1, j);
                q_bar = M(i+1, j);
            else
                q = M(i-1, j-1);
                q_bar = M(i+1, j+1);
            end
        
            if M(i, j) < max(q, q_bar)
                M(i, j) = 0;
            end
        end
    end
    
    figure();
    imshow(M);


    Ih = (M >= Th);
    figure();
    imshow(Ih);
    Il = (M >= Tl);
    figure();
    imshow(Il);

    final_edge_map = zeros(size(Ih));

    [rows, cols] = size(Ih);


    for i = 1:rows
        for j = 1:cols
            if Ih(i, j) == 1
                edge_direction = atan2(dy(i,j), dx(i,j));
            
                [~, index] = min(abs(edge_direction - directions));
            
                switch index
                    case 1
                        x_step = -1;
                        y_step = 0;
                    case 2
                        x_step = -1;
                        y_step = -1;
                    case 3
                        x_step = 0;
                        y_step = -1;
                    case 4
                        x_step = 1;
                        y_step = -1;
                    case 5
                        x_step = 1;
                        y_step = 0;
                    case 6
                        x_step = 1;
                        y_step = 1;
                    case 7
                        x_step = 0;
                        y_step = 1;
                    case 8
                        x_step = -1;
                        y_step = 1;
                end
            
                x = i;
                y = j;
            
                while Ih(x, y) == 1 
                    final_edge_map(x, y) = 1;
                    
                    x = x + x_step;
                    y = y + y_step;
                    
                    if x < 1 || x > rows || y < 1 || y > cols
                        break;
                    end

                    final_edge_map(x, y) = Il(x,y);

                end
            end
        end
    end

    edge_img = final_edge_map;

end
