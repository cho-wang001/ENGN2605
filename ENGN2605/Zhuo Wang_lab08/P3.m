clc;
clear;
close all;

K = load("Question2\CalibrationMatrix.mat").K;

R = load("Question2\RelativePose.mat").rotation;
T = load("Question2\RelativePose.mat").translation;

% R21 = R12';T21=-R12'*T12
R21 = R';
T21 = -R'*T;

C1 = load("Question2\Two_View_Correspondences.mat").img1_points;
C2 = load("Question2\Two_View_Correspondences.mat").img2_points;

g1 = zeros(3, size(C1, 1));
g2 = zeros(3, size(C1, 1));

p1 = zeros(3, size(C1, 1));
p11 = zeros(3, size(C1, 1));
p2 = zeros(3, size(C1, 1));
p21 = zeros(3, size(C1, 1));
P1 = zeros(3, size(C1, 1));
P2 = zeros(3, size(C1, 1));

gg1 = zeros(3, size(C1, 1));
gg2 = zeros(3, size(C1, 1));

ik = inv(K);

for i = 1:size(C1, 1)
    g1(:, i) = ik*[C1(i, :), 1]';
    g2(:, i) = ik*[C2(i, :), 1]';

    p = pinv([-R*g1(:, i), g2(:, i)]) * T;

    p1(:, i) = p(1).* g1(:, i);
    p2(:, i) = R*p1(:, i) + T;
    
    p111 = pinv([-R21*g2(:, i), g1(:, i)]) * T21;

    p21(:, i) = p111(1).* g2(:, i);
    p11(:, i) = R21*p21(:, i) + T21;

    P1(:, i) = (p1(:, i)+p11(:, i))./2;
    P2(:, i) = (p2(:, i)+p21(:, i))./2;
end


for i = 1:size(P1, 2)
    GG1 = K*P1(:, i);
    gg1(:, i) = GG1./GG1(3);
    GG2 = K*P2(:, i);
    gg2(:, i) = GG2./GG2(3);
end

gg1 = gg1(1:2,:)';
gg2 = gg2(1:2,:)';

skew_T = [0, -T(3), T(2); T(3), 0, -T(1); -T(2), T(1), 0];

E = skew_T * R;

X = [gg1(1,1), gg1(1,2), gg2(1,1), gg2(1,2), 1];

Energy = @(x) (evaluateEnergyFunction(x, C1(1,:), C2(1,:), E));
options = optimoptions('lsqnonlin', 'Display', 'iter');
options.Algorithm = 'levenberg-marquardt';

optimal_points = lsqnonlin(Energy, X, [], [], options);


oc1 = optimal_points(1:2);
oc2 = optimal_points(3:4);


%mid-point triangulation method
MP1 = K*p11(:,1);
tmp1 = MP1./MP1(3);
tmp1 = tmp1(1:2)';

MP2 = K*p2(:,1);
tmp2 = MP2./MP2(3);
tmp2 = tmp2(1:2)';



errors_m = norm(tmp1 - gg1(1,:))^2 + norm(tmp2 - gg2(1,:))^2;
disp("Reprojection errors (the mid-point triangulation method):"+errors_m);
errors_o = norm(oc1 - gg1(1,:))^2 + norm(oc2 - gg2(1,:))^2;
disp("Reprojection errors (the non-linear optimization method):"+errors_o);



function Energy = evaluateEnergyFunction(x, gamma, gamma_bar, E)
    gamma_hat = x(1:2);
    gamma_bar_hat = x(3:4);
    lambda = x(5);

    term1 = norm(gamma_hat - gamma)^2;
    term2 = norm(gamma_bar_hat - gamma_bar)^2;

    term3 = lambda.*([gamma_bar_hat, 1]*E* [gamma_hat, 1]');

    Energy = term1 + term2 + term3;
end


