clc; clear all;
obj=VideoReader('video.mp4'); 
vidFrames = read(obj); 
numFrames = get(obj,'numberOfFrames');  
dt = get(obj,'FrameRate');
%%
X = [];
for k = 1 : numFrames
gray(:, :) = double(rgb2gray(vidFrames(:,:,:,k)));
[width, length] = size(gray);
X = [X, reshape(gray,[width*length, 1])];
end

%% 
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
[U, S, V] = svd(X1, 'econ');
r = size(U, 2); % rank truncation
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);
%% Build Atilde and DMD Modes
Atilde = Ur'*X2*Vr/Sr;
[W, D] = eig(Atilde);
Phi = X2*Vr/Sr*W; % DMD Modes
%% DMD Spectra
t = 1:numFrames;
dt = t(2) - t(1);
% t = linspace(0, numFrames*dt, numFrames);
%%
lambda = diag(D);
omega = log(lambda)/dt;

%% Seperation
[omega_bg, min_index] = min(abs(omega));
Phi_bg = Phi(:,min_index);
%% Compute DMD Solution
%% X_bg
b = Phi_bg\X(:, 1);

time_dynamics_bg = [];
for iter = 1: numFrames
    time_dynamics_bg= [time_dynamics_bg, b*10*exp(omega_bg*t(iter))];
end
X_bg = Phi_bg * time_dynamics_bg;

%% X_fg
X_bg(X_bg < 0) = 0;
X_fg = X-X_bg;
X_fg(X_fg < 0) = 0;
%% Play Video
for k = 1: numFrames
    colormap gray;
    imagesc(reshape(abs(X_bg(:, 2)), [width, length]));
end
%%
for k = 1:numFrames
    colormap gray;
    [w, l] = size(time_dynamics_bg);
    imagesc(reshape(abs(time_dynamics_bg(:, k)), [w, l]));
end


