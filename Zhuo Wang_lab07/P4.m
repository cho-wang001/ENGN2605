clc;
clear;
close all;
format long;

d1 = load("Q4_data\Two_View_Correspondences.mat").img1_points;
d2 = load("Q4_data\Two_View_Correspondences.mat").img2_points;
k4 = load("Q4_data\Calibration_Matrix.mat").K;

N = 45;

r1 = [];
r2 = [];
t1 = [];
t2 = [];

for i = 1:45
    sl = size(d1, 1);
    sj = randsample(sl, 5,'true');
    
    for f = sj
        x_coord1 = d1(f,1);
        y_coord1 = d1(f,2);
        x_coord2 = d2(f,1);
        y_coord2 = d2(f,2);
    end
    
    cc1 = [x_coord1, y_coord1];
    cc2 = [x_coord2, y_coord2];
    
    c1m = [];
    c2m = [];
    
    for i = 1:size(cc1, 1)
        g1 = k4*[cc1(i, 1);cc1(i, 2);1]; 
        g1 = g1';
        c1m = [c1m; g1];
        
        g2 = k4*[cc2(i, 1);cc2(i, 2);1];
        g2 = g2';
        c2m = [c2m; g2];
    end
    
    mi = cat(3, c1m, c2m);
    
    es = fivePointAlgorithmSelf(mi);
    
    for t = 1:size(es, 3)
        [U, S, V] = svd(es{:,:,t});
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

        r1 = cat(3, r1, R1);
        r2 = cat(3, r2, R2);
        t1 = [t1,T1];
        t2 = [t2,T2];

    end
end

[rg, tg] = getVeridicalGroupsOfPoses(r1, r2, t1, t2);

N = size(tg, 2);
combinations = nchoosek(1:N, 2);
numCombinations = size(combinations, 1);

deltaR = zeros(numCombinations, 1);
deltaT = zeros(numCombinations, 1);

for i = 1:numCombinations
    idx1 = combinations(i, 1);
    idx2 = combinations(i, 2);

    deltaR(i) = real(acos(0.5 * (trace(rg( :, :,idx1) * transpose(rg( :, :, idx2)) - eye(3)))));
    deltaT(i) = abs(dot(tg(:, idx1), tg(:, idx2)) - 1);
end

figure;
subplot(2, 1, 1);
histogram(deltaR, 'Normalization', 'probability', 'BinWidth', 0.1);
title('Rotation Differences');


subplot(2, 1, 2);
histogram(deltaT, 'Normalization', 'probability', 'BinWidth', 0.1);
title('Translation Differences');
