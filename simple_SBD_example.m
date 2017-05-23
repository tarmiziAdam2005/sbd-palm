clc; clear;


%% I. SIMULATE DATA FOR SBD:
%  =========================
% For the simple example, we will generate a single-slice, random kernel
% with small problem sizes for faster solution.

%% 1. Random kernel and observation
k = [8 8];           	% kernel size

% GENERATE
A0 = proj2oblique(randn(k));

%% 2. Activation map
m = [64 64];            % image size for each slice / observation grid

%   Each pixel has probability theta of being a kernel location
theta = 3e-3;           % activation concentration
eta = 0; 1e-3;             % additive noise variance

% GENERATE
X0_good = false;
while ~X0_good
    X0 = double(rand(m) <= theta);              % activations are on / off
    X0_good = sum(X0(:) ~= 0) > 0;
end

Y = convfft2(A0, X0) + sqrt(eta)*randn(m);     % observation

%% II. Sparse Blind Deconvolution:
%  ===============================
[Aout, Xout, costs] = AXsolve_PALM( Y, k, 0.1, 100 );

figure(1); plot(0:size(costs,1)-1, costs(:,1), '.'); hold on;
plot((0:size(costs,1)-1)+0.5, costs(:,2), 'r.'); hold off;
xlim([0 100]);

figure(2); subplot(121); imagesc(abs(A0)); subplot(122); imagesc(abs(Aout));
figure(3); subplot(121); imagesc(abs(X0)); subplot(122); imagesc(abs(Xout));