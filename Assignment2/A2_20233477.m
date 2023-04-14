%% 1
disp("INITIALS AND COLOR TRANSFORMATION (5 PTS)");

% read and convert to double [0,1]
baby2 = VideoReader("./data/baby2.mp4");
H = baby2.Height;
W = baby2.Width;
frames = baby2.numframes;

baby2_frames = zeros(H, W, 3, frames);
for i = 1:frames
    baby2_frames(:,:,:,i) = double(readFrame(baby2)) / 255;
end

% convert to YIQ color space
baby2_yiq = zeros(H, W, 3, frames);
for i = 1:frames
    baby2_yiq(:,:,:,i) = rgb2ntsc(baby2_frames(:,:,:,i));
end

% plot
imwrite(baby2_frames(:,:,:,1), "./figure/initial_face.jpg");
imwrite(baby2_yiq(:,:,:,1), "./figure/initial_YIQ.jpg");


%% 2
disp("LAPLACIAN PYRAMID (20PTS)");

% for each frames
for i = 1:frames
    
    % Gaussian pyramid
    % https://kr.mathworks.com/help/images/ref/impyramid.html
    G0 = baby2_yiq(:,:,:,i); % original frame
    G1 = impyramid(G0, 'reduce');
    G2 = impyramid(G1, 'reduce');
    G3 = impyramid(G2, 'reduce');
    G4 = impyramid(G3, 'reduce');
    
    % upsample and residual
    L4 = G4;
    L3 = imresize(G4, 2) - G3;
    L2 = imresize(G3, 2) - G2;
    L1 = imresize(G2, 2) - G1;
    L0 = imresize(G1, 2) - G0;
    
    % for interpreability, I use mat2gray
    if i == 1
        imwrite(mat2gray(G4), './figure/G4.png');
        imwrite(mat2gray(G3), './figure/G3.png');
        imwrite(mat2gray(G2), './figure/G2.png');
        imwrite(mat2gray(G1), './figure/G1.png');
        imwrite(mat2gray(G0), './figure/G0.png');
        
        imwrite(mat2gray(L4), './figure/L4.png');
        imwrite(mat2gray(L3), './figure/L3.png');
        imwrite(mat2gray(L2), './figure/L2.png');
        imwrite(mat2gray(L1), './figure/L1.png');
        imwrite(mat2gray(L0), './figure/L0.png');
    end
end

%% 3
disp("TEMPORAL FILTERING (30PTS)");








