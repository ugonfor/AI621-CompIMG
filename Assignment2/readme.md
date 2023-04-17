# Assignment#2 Eulerian Video Magnification

> Late Day: 3

Author: Hyogon Ryu (20233477)



All figures images files can find on here

## INITIALS AND COLOR TRANSFORMATION (5 PTS)

### Description

I load the face.mp4 and baby2.mp4

### Code

```matlab
%% 1
disp("INITIALS AND COLOR TRANSFORMATION (5 PTS)");

% read and convert to double [0,1]
baby2 = VideoReader("./data/face.mp4");
H = baby2.Height;
W = baby2.Width;
num_frames = baby2.numframes;
Fs = baby2.FrameRate;

baby2_frames = zeros(H, W, 3, num_frames);
for i = 1:num_frames
    baby2_frames(:,:,:,i) = double(readFrame(baby2)) / 255;
end

% convert to YIQ color space
baby2_yiq = zeros(H, W, 3, num_frames);
for i = 1:num_frames
    baby2_yiq(:,:,:,i) = rgb2ntsc(baby2_frames(:,:,:,i));
end

% plot
imwrite(baby2_frames(:,:,:,1), "./figure/initial_face.jpg");
imwrite(baby2_yiq(:,:,:,1), "./figure/initial_YIQ.jpg");

```

### images

![initial_face](.\figure\initial_face.jpg)

`./figure/initial_face.jpg`



![initial_face](.\figure\initial_yiq.jpg)

`.\figure\initial_yiq.jpg`





## LAPLACIAN PYRAMID (20PTS)

### Description

By using Gaussian Filter, I reduce the size of each frame. (use `impyramid` function)

Then, get the residual frame(Laplacian frame) by substitute up-sampled Gaussian Images from the original frame. (`Laplacian_k = Gaussian_k - imresize(Gaussian_k+1, 2)` )

### Code

```matlab

disp("LAPLACIAN PYRAMID (20PTS)");

L4_frames = zeros(H/16,W/16,3,num_frames);
L3_frames = zeros(H/8,W/8,3,num_frames);
L2_frames = zeros(H/4,W/4,3,num_frames);
L1_frames = zeros(H/2,W/2,3,num_frames);
L0_frames = zeros(H,W,3,num_frames);

% for each frames
for i = 1:num_frames
    
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
    
    L4_frames(:,:,:,i) = L4;
    L3_frames(:,:,:,i) = L3;
    L2_frames(:,:,:,i) = L2;
    L1_frames(:,:,:,i) = L1;
    L0_frames(:,:,:,i) = L0;
    
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

clear('L0');clear('L1');clear('L2');clear('L3');clear('L4');
clear('G0');clear('G1');clear('G2');clear('G3');clear('G4');

```

### Images

![](./figure/initial_YIQ.jpg)

`./figure/Initial_YIQ.jpg`

![L0](.\figure\L0.png)

`.\figure\L0.png`



![L0](.\figure\L1.png)

`.\figure\L1.png`



![L0](.\figure\L2.png)

`.\figure\L2.png`



![L0](.\figure\L3.png)

`.\figure\L3.png`



![L0](.\figure\L4.png)

`.\figure\L4.png`



## TEMPORAL FILTERING (30PTS)

### Description

By using FFT, I check the frequency domain of the video.

I filter the temporal component of the video, and use the FFT.

after that, I compare the frequent domain and butterworthBandpassFilter. (To determine the hyperparameter of the bandpass filter. I follow the paper's setting and It seems good.)

### Code

```matlab
%% 3
disp("TEMPORAL FILTERING (30PTS)");

[H, W, C, N] = size(L0_frames);

temporal_tmp = zeros(N, 1);
temporal_fft = zeros(fix(N/2) + 1 , 1);
for h = 1:H
    for w = 1:W
        for c = 1:C
            temporal_tmp(:,1) = L0_frames(h,w,c,:);
            fftx = fft(temporal_tmp);

            tmp1 = abs(fftx/N);
            tmp2 = tmp1(1:fix(N/2)+1);
            tmp2(2:end-1) = 2*tmp2(2:end-1);
            
            % average
            temporal_fft(:,1) = temporal_fft(:,1) + tmp2(:,1);
        end
    end
end

figure;
subplot(2, 1, 1);
freq = Fs * (0:fix(N/2)) / N;
plot(freq, temporal_fft);

subplot(2, 1, 2);
N = 256;
Fc1 = 0.83;
Fc2 = 1;
% Fc1 = 2.33;
% Fc2 = 2.67;
Hd = butterworthBandpassFilter(Fs, N, Fc1, Fc2);
fftHd = freqz(Hd, size(freq,2));
plot(freq, abs(fftHd));
```

### Images

![](C:\Users\ryuhy\Desktop\coursework\AI621 CompIMG\Assignment #2\figure\TEMPORAL_FILTERING.jpg)



## EXTRACTING THE FREQUENCY BAND OF INTEREST (30PTS)

### Description

I decided the filter just before. So I extract the frequency domain of the video by using that filter.

### code

```matlab
%% 4
disp("EXTRACTING THE FREQUENCY BAND OF INTEREST (30PTS)");


fftHd = freqz(Hd, num_frames);
L0_filtered = Temporal_Filtering(L0_frames, fftHd);
L1_filtered = Temporal_Filtering(L1_frames, fftHd);
L2_filtered = Temporal_Filtering(L2_frames, fftHd);
L3_filtered = Temporal_Filtering(L3_frames, fftHd);
L4_filtered = Temporal_Filtering(L4_frames, fftHd);
```

```matlab

function [ res ] = Temporal_Filtering(frames, fftHd)
    [H,W,C,N] = size(frames);
    
    res = zeros(H,W,C,N);
    
    fftHd_pad = zeros(2*N, 1);
    fftHd_pad(1:N) = fftHd;
    fftHd_pad(N+1:end) = fftHd(N:-1:1);
    
    for h = 1:H
        for w = 1:W
            pixel = frames(h,w,1,:);
            fft_p = fft(pixel, 2*N);
            tmp = real(ifft(fft_p .* reshape(fftHd_pad, [1,1,1,2*N])));
            res(h,w,1,:) = tmp(:,:,:,1:N);
            
            pixel = frames(h,w,2,:);
            fft_p = fft(pixel, 2*N);
            tmp = real(ifft(fft_p .* reshape(fftHd_pad, [1,1,1,2*N])));
            res(h,w,2,:) = tmp(:,:,:,1:N);
            
            pixel = frames(h,w,3,:);
            fft_p = fft(pixel, 2*N);
            tmp = real(ifft(fft_p .* reshape(fftHd_pad, [1,1,1,2*N])));
            res(h,w,3,:) = tmp(:,:,:,1:N);
        end
    end
end
```

### Images

No image here.



## IMAGE RECONSTRUCTION (20PTS, WHICH INCLUDES EVALUATION OF YOUR RESULTS)

### Description

Image reconstruction by adding the original video and filtered Laplacian frames. 

To magnitude the some signal, there are coefficient of the filter. (in code, `a` is the coefficient)

I select the channel 3 (in code, `c` is the selected channel), since changing third-component of the frames brings the color change. 



And then, I made the video with reconstructed frames.

### Code

```matlab
%% 5
disp("IMAGE RECONSTRUCTION");

a = [80 80 80 80 80];
c = 3;

RECONSTRUCTED = zeros(H, W, 3, num_frames);
for n = 1:num_frames
    RECONSTRUCTED_frame = baby2_yiq(:,:,c,n) + a(1) * L0_filtered(:,:,c,n) + a(2) * imresize(L1_filtered(:,:,c,n), [H,W])  ...
       + a(3) * imresize(L2_filtered(:,:,c,n), [H,W]) + a(4) * imresize(L3_filtered(:,:,c,n), [H,W]) + a(5) * imresize(L4_filtered(:,:,c,n), [H,W]);

    RECONSTRUCTED(:,:,c,n) = RECONSTRUCTED_frame(:,:,1);
end
RECONSTRUCTED(:,:,1,:) = baby2_yiq(:,:,1,:); % change as c
RECONSTRUCTED(:,:,2,:) = baby2_yiq(:,:,2,:); % change as c
% RECONSTRUCTED(:,:,3,:) = baby2_yiq(:,:,3,:); % change as c


%% 6
disp("MAKE VIDEO");

baby2_rgb_RECON = zeros(H, W, 3, num_frames);
for i=1:N
    frame = ntsc2rgb(RECONSTRUCTED(:,:,:,i));
    frame(frame > 1) = 1;
    frame(frame < 0 ) = 0;
    baby2_rgb_RECON(:,:,:,i) = frame;
end

vw1 = VideoWriter('./figure/face_RECON.avi');
open(vw1);

for i=1:N
    writeVideo(vw1, squeeze(baby2_rgb_RECON(:, :, :, i)));
end

close(vw1);

disp("DONE!");
```

### Video

