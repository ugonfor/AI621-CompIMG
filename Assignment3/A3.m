
%% Toy Problem
clear
% Load the source image
s = imread('./data/toy_problem.png');
s = im2double(s);
[N, M] = size(s);

index = zeros(N, M);
index(1:N*M) = 1: N*M;

A = zeros( (N-1)*M + N*(M-1) + 1, N*M);
b = zeros((N-1)*M + N*(M-1) + 1, 1);

row = 1;
for n=1:N
    for m=1:M-1
        A(row, index(n,m+1)) = 1;
        A(row, index(n,m)) = -1;
        b(row) = s(n, m+1) - s(n,m);
        row = row + 1;
    end
end
for n=1:N-1
    for m=1:M
        A(row, index(n+1,m)) = 1;
        A(row, index(n,m)) = -1;
        b(row) = s(n+1, m) - s(n,m);
        row = row + 1;
    end
end
A(row, index(1,1)) = 1;
b(row) = s(1,1);

A = sparse(A);

v = A\b;
v_img = reshape(v, [N,M]);
figure; imshow(v_img); title("TOY PROBLEM");

disp(['Error: ' num2str(sqrt(sum(s(:)-v_img(:))))])

%% Poisson Blending

% starter script for project 3
DO_BLEND = true;
DO_MIXED  = false;
DO_COLOR2GRAY = false;

if DO_BLEND
    im_background = imresize(im2double(imread('./data/hiking.jpg')), 0.5, 'bilinear');
    im_object = imresize(im2double(imread('./data/penguin-chick.jpeg')), 0.5, 'bilinear');

    % get source region mask from the user
    objmask = getMask(im_object);
    % align im_s and mask_s with im_background
    [im_s, mask_s] = alignSource(im_object, objmask, im_background);

    % blend
    im_blend = poissonBlend(im_s, mask_s, im_background);
    figure(3), hold off, imshow(im_blend), title("RESULT");
end

if DO_MIXED
    % read images
    %...
    
    % blend
    im_blend = mixedBlend(im_s, mask_s, im_bg);
    figure(3), hold off, imshow(im_blend);
end

if DO_COLOR2GRAY
    im_rgb = im2double(imread('./samples/colorBlindTest35.png'));
    im_gr = color2gray(im_rgb);
    figure(4), hold off, imagesc(im_gr), axis image, colormap gray
end


