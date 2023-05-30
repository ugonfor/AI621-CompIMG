%%% HDR Imaging
clear;

%% Load
image = imread("./data/door.tiff");
image = im2double(image);

% ratio
ratio = 0.01;

jpg_img = imread("./data/exposure_stack/exposure1.jpg");
jpg_img = imresize(jpg_img, ratio);
jpg_imgs = zeros([size(jpg_img), 16], 'uint8');

k = 16; % 16 images
for i = 1:k
    img_path = sprintf("./data/exposure_stack/exposure%d.jpg", i);
    jpg_imgs(:,:,:,i) = imresize(imread(img_path), ratio);
end

[H, W, C] = size(jpg_img);

%% Linearize Rendered Images

Z_r = zeros(H*W,k);
Z_g = zeros(H*W,k);
Z_b = zeros(H*W,k);

for i=1:k
    tmp = jpg_imgs(:,:,1,i);
    Z_r(:,i) = tmp(:);
    tmp = jpg_imgs(:,:,2,i);
    Z_g(:,i) = tmp(:);
    tmp = jpg_imgs(:,:,3,i);
    Z_b(:,i) = tmp(:);
end

% exposure
t = zeros(16,1);
for i=1:k
    t(i) = log(1/2048 * pow2(i-1));
end

w_uni = zeros(256,1);
for i=1:256
    w_uni(i) = 1;
end

w_tent = zeros(256,1);
for i=1:128
    w_tent(i) = ( 2*(i-1) )/256;
end
for i=1:128
    w_tent(i + 128) = 1 - ( 2*(i-1) )/256;
end

% w_gau = zeros(256,1);
%for i=1:256
%    w_gau(i) = % gaussian
%end

[g_r_uni, l_r_uni] = minimize(Z_r, t, w_uni);
[g_g_uni, l_g_uni] = minimize(Z_g, t, w_uni);
[g_b_uni, l_b_uni] = minimize(Z_b, t, w_uni);

[g_r_tent, l_r_tent] = minimize(Z_r, t, w_tent);
[g_g_tent, l_g_tent] = minimize(Z_g, t, w_tent);
[g_b_tent, l_b_tent] = minimize(Z_b, t, w_tent);

%[g_r_gau, l_r_gau] = minimize(Z_r, t, w_gau);
%[g_g_gau, l_g_gau] = minimize(Z_g, t, w_gau);
%[g_b_gau, l_b_gau] = minimize(Z_b, t, w_gau);

figure(1);
plot(g_r_uni, 'r');
hold on;
plot(g_g_uni, 'g');
hold on;
plot(g_b_uni, 'b');
hold on;
title('uniform');

figure(2);
plot(g_r_tent, 'r');
hold on;
plot(g_g_tent, 'g');
hold on;
plot(g_b_tent, 'b');
hold on;
title('tent');

% linearize
linearize_img_uni = zeros(size(jpg_imgs));
linearize_img_tent = zeros(size(jpg_imgs));
for i=1:k
    for c=1:C
        for h=1:H
            for w=1:W
                if c==1
                    linearize_img_uni(h,w,c,i) = exp(g_r_uni(jpg_imgs(h,w,c,i)+1));
                    linearize_img_tent(h,w,c,i) = exp(g_r_tent(jpg_imgs(h,w,c,i)+1));
                end
                if c==2
                    linearize_img_uni(h,w,c,i) = exp(g_g_uni(jpg_imgs(h,w,c,i)+1));
                    linearize_img_tent(h,w,c,i) = exp(g_g_tent(jpg_imgs(h,w,c,i)+1));
                end
                if c==3
                    linearize_img_uni(h,w,c,i) = exp(g_b_uni(jpg_imgs(h,w,c,i)+1));
                    linearize_img_tent(h,w,c,i) = exp(g_b_tent(jpg_imgs(h,w,c,i)+1));
                end
            end
        end
    end
end

%save('./result/linearize_tent.mat', "linearize_img_tent");
%save('./result/linearize_uni.mat', "linearize_img_uni");

%%
clear
                
%% Merge Exposure Stack Into HDR Image
ratio = 0.03;
k = 16;

t = zeros(16,1);
for i=1:k
    t(i) = 1/2048 * pow2(i-1);
end

% weight
w_uni = zeros(256,1);
for i=1:256
    w_uni(i) = 1;
end

w_tent = zeros(256,1);
for i=1:128
    w_tent(i) = ( 2*(i-1) )/256;
end
for i=1:128
    w_tent(i + 128) = 1 - ( 2*(i-1) )/256;
end

tiff_tmp = imresize(imread("./data/exposure_stack/exposure1.tiff"), ratio);
[raw_h, raw_w, raw_c] = size(tiff_tmp);

raw_linear_uni_img = zeros([raw_h, raw_w, raw_c]);
raw_linear_tent_img = zeros([raw_h, raw_w, raw_c]);
raw_log_uni_img = zeros([raw_h, raw_w, raw_c]);
raw_log_tent_img = zeros([raw_h, raw_w, raw_c]);

jpg_tmp = imresize(imread("./data/exposure_stack/exposure1.jpg"), ratio);
[rendered_h, rendered_w, rendered_c] = size(jpg_tmp);

rendered_linear_uni_img = zeros([rendered_h, rendered_w, rendered_c]);
rendered_linear_tent_img = zeros([rendered_h, rendered_w, rendered_c]);
rendered_log_uni_img = zeros([rendered_h, rendered_w, rendered_c]);
rendered_log_tent_img = zeros([rendered_h, rendered_w, rendered_c]);


% init raw data
imgs_raw = zeros([size(tiff_tmp), 16]);
imgs_raw_ldr = zeros([size(tiff_tmp), 16], 'uint8');
for i = 1:k
    img_path = sprintf("./data/exposure_stack/exposure%d.tiff", i);
    imgs_raw(:,:,:,i) = im2double(imresize(imread(img_path), ratio));
    imgs_raw_ldr(:,:,:,i) = uint8(imresize(imread(img_path), ratio));
end


% init rendered data
imgs_rendered = zeros([size(jpg_tmp), 16]);
imgs_rendered_ldr = zeros([size(jpg_tmp), 16], 'uint8');
for i = 1:k
    img_path = sprintf("./data/exposure_stack/exposure%d.jpg", i);
    imgs_rendered(:,:,:,i) = im2double(imresize(imread(img_path), ratio));
    imgs_rendered_ldr(:,:,:,i) = uint8(imresize(imread(img_path), ratio));
end
imgs_rendered_ldr_uni = uint8( load("./result/linearize_uni.mat").linearize_img_uni * 256 );
imgs_rendered_ldr_tent = uint8( load("./result/linearize_tent.mat").linearize_img_tent * 256 );

% raw 
for c=1:raw_c
    for h=1:raw_h
        for w=1:raw_w
            
            % 분자
            linear_uni_top = 0;
            linear_tent_top = 0;
            log_uni_top = 0;
            log_tent_top = 0;
            
            % 분모
            linear_uni_bottom = 0;
            linear_tent_bottom = 0;
            log_uni_bottom = 0;
            log_tent_bottom = 0;
            
            for i=1:k
               linear_uni_top = linear_uni_top + w_uni( imgs_raw_ldr(h,w,c,i) + 1) * imgs_raw(h,w,c,i) / t(i);
               linear_tent_top = linear_tent_top + w_tent( imgs_raw_ldr(h,w,c,i) + 1) * imgs_raw(h,w,c,i) / t(i);
               log_uni_top = log_uni_top + w_uni( imgs_raw_ldr(h,w,c,i) + 1) * ( log(imgs_raw(h,w,c,i)) - log(t(i)) );
               log_tent_top = log_tent_top + w_tent( imgs_raw_ldr(h,w,c,i) + 1) * ( log(imgs_raw(h,w,c,i)) - log(t(i)) );              
               
               linear_uni_bottom = linear_uni_bottom + w_uni( imgs_raw_ldr(h,w,c,i) + 1);
               linear_tent_bottom = linear_tent_bottom + w_tent( imgs_raw_ldr(h,w,c,i) + 1);
               log_uni_bottom = log_uni_bottom + w_uni( imgs_raw_ldr(h,w,c,i) + 1);
               log_tent_bottom = log_tent_bottom + w_tent( imgs_raw_ldr(h,w,c,i) + 1);
            end
            
            raw_linear_uni_img(h,w,c) = linear_uni_top / linear_uni_bottom;
            raw_linear_tent_img(h,w,c) = linear_tent_top / linear_tent_bottom;
            raw_log_uni_img(h,w,c) = exp(log_uni_top / log_uni_bottom);
            raw_log_tent_img(h,w,c) = exp(log_tent_top / log_tent_bottom);
            
        end
    end
end


% render
for c=1:rendered_c
    for h=1:rendered_h
        for w=1:rendered_w
            
            % 분자
            linear_uni_top = 0;
            linear_tent_top = 0;
            log_uni_top = 0;
            log_tent_top = 0;
            
            % 분모
            linear_uni_bottom = 0;
            linear_tent_bottom = 0;
            log_uni_bottom = 0;
            log_tent_bottom = 0;
            
            for i=1:k
               linear_uni_top = linear_uni_top + w_uni( imgs_rendered_ldr_uni(h,w,c,i) + 1) * imgs_rendered(h,w,c,i) / t(i);
               linear_tent_top = linear_tent_top + w_tent( imgs_rendered_ldr_tent(h,w,c,i) + 1) * imgs_rendered(h,w,c,i) / t(i);
               log_uni_top = log_uni_top + w_uni( imgs_rendered_ldr_uni(h,w,c,i) + 1) * imgs_rendered(h,w,c,i) - log(t(i));
               log_tent_top = log_tent_top + w_tent( imgs_rendered_ldr_tent(h,w,c,i) + 1) * imgs_rendered(h,w,c,i) - log(t(i));              
               
               linear_uni_bottom = linear_uni_bottom + w_uni( imgs_rendered_ldr_uni(h,w,c,i) + 1);
               linear_tent_bottom = linear_tent_bottom + w_tent( imgs_rendered_ldr_tent(h,w,c,i) + 1);
               log_uni_bottom = log_uni_bottom + w_uni( imgs_rendered_ldr_uni(h,w,c,i) + 1);
               log_tent_bottom = log_tent_bottom + w_tent( imgs_rendered_ldr_tent(h,w,c,i) + 1);
            end
            
            rendered_linear_uni_img(h,w,c) = linear_uni_top / linear_uni_bottom;
            rendered_linear_tent_img(h,w,c) = linear_tent_top / linear_tent_bottom;
            rendered_log_uni_img(h,w,c) = exp(log_uni_top / log_uni_bottom);
            rendered_log_tent_img(h,w,c) = exp(log_tent_top / log_tent_bottom);
            
        end
    end
end

imwrite(tonemap(raw_linear_uni_img), "./result/raw_linear_uni_img.png");
imwrite(tonemap(raw_linear_tent_img), "./result/raw_linear_tent_img.png");
imwrite(tonemap(raw_log_uni_img), "./result/raw_log_uni_img.png");
imwrite(tonemap(raw_log_tent_img), "./result/raw_log_tent_img.png");
imwrite(tonemap(rendered_linear_uni_img), "./result/rendered_linear_uni_img.png");
imwrite(tonemap(rendered_linear_tent_img), "./result/rendered_linear_tent_img.png");
imwrite(tonemap(rendered_log_uni_img), "./result/rendered_log_uni_img.png");
imwrite(tonemap(rendered_log_tent_img), "./result/rendered_log_tent_img.png");

hdrwrite(raw_linear_uni_img, "./result/raw_linear_uni_img.hdr");
hdrwrite(raw_linear_tent_img, "./result/raw_linear_tent_img.hdr");
hdrwrite(raw_log_uni_img, "./result/raw_log_uni_img.hdr");
hdrwrite(raw_log_tent_img, "./result/raw_log_tent_img.hdr");
hdrwrite(rendered_linear_uni_img, "./result/rendered_linear_uni_img.hdr");
hdrwrite(rendered_linear_tent_img, "./result/rendered_linear_tent_img.hdr");
hdrwrite(rendered_log_uni_img, "./result/rendered_log_uni_img.hdr");
hdrwrite(rendered_log_tent_img, "./result/rendered_log_tent_img.hdr");

%% Evaluation

% HDR 이미지 로드
hdrImage1 = hdrread('raw_linear_uni_img.hdr');
hdrImage2 = hdrread('raw_linear_tent_img.hdr');
% 색 체커 좌표 로드 (ginput으로 얻은 좌표)
load('color_checker_coordinates.mat');

% 평균 조도 값 저장을 위한 배열 초기화
avgLuminance1 = zeros(6, 1);
avgLuminance2 = zeros(6, 1);

% HDR 이미지 1 처리
xyzImage1 = rgb2xyz(hdrImage1, 'Colorspace', 'linear-rgb');
luminance1 = xyzImage1(:, :, 2); % Y 채널 추출
for i = 1:6
    % 각 중립 패치에 대해 사각형 영역 자르기
    patchRegion = imcrop(luminance1, colorCheckerCoordinates{i});
    % 평균 조도 계산
    avgLuminance1(i) = mean(patchRegion(:));
end

% HDR 이미지 2 처리 (위와 유사한 단계)
xyzImage2 = rgb2xyz(hdrImage2, 'Colorspace', 'linear-rgb');
luminance2 = xyzImage2(:, :, 2);
for i = 1:6
    patchRegion = imcrop(luminance2, colorCheckerCoordinates{i});
    avgLuminance2(i) = mean(patchRegion(:));
end

% 평균 조도 값의 로그 계산
logAvgLuminance1 = log(avgLuminance1);
logAvgLuminance2 = log(avgLuminance2);

% 선형 회귀 수행
coeffs1 = polyfit(logAvgLuminance1, 1:6, 1);
lineFit1 = polyval(coeffs1, logAvgLuminance1);
coeffs2 = polyfit(logAvgLuminance2, 1:6, 1);
lineFit2 = polyval(coeffs2, logAvgLuminance2);

% 최소 제곱 오차 계산
error1 = sum((lineFit1 - (1:6)).^2);
error2 = sum((lineFit2 - (1:6)).^2);

% 로그 플롯 생성
figure;
plot(logAvgLuminance1, 1:6, 'ro', 'MarkerSize', 8);
hold on;
plot(logAvgLuminance2, 1:6, 'bo', 'MarkerSize', 8);
plot(lineFit1, 1:6, 'r--');
plot(lineFit2, 1:6, 'b--');
hold off;
grid on;
xlabel('Log Average Luminance');
ylabel('Patch Number');
legend('HDR Image 1', 'HDR Image 2', 'Linear Fit 1', 'Linear Fit 2');
title('Comparison of Average Luminances');

% 오차 출력
disp('Error for HDR Image 1:');
disp(error1);
disp('Error for HDR Image 2:');
disp(error2);