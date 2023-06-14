%% LIGHT FIELD RENDERING, FOCAL STACKS, AND DEPTH FROM DEFOCUS (100 POINTS)
clear;

%% Initials (5 points)

img = imread("./data/chessboard_lightfield.png");
%img = imresize(img, 0.4);
img = im2double(img);
[h, w, c] = size(img);


u = 16;
v = 16;
s = h/u;
t = w/v;
c = c;
LF = zeros(u,v,s,t,c);
disp(size(LF));
for i=1:h
    for j=1:w
        for k=1:c
            LF( mod(i-1, u)+1, mod(j-1, v)+1, ...
                floor((i-1)/u)+1, floor((j-1)/v)+1, ...
                k) = img(i,j,k);
        end
    end
end

%% Sub aperture view (20 points)
SB_view = zeros(s*u,t*v,c);
for i=1:u
    for j=1:v
        SB_view( (s*(i-1))+1:(s*(i)), ...
                 (t*(j-1))+1:(t*(j)),:) = LF(i,j,:,:,:);
    end
end

imshow(SB_view);
imwrite(SB_view,"Sub aperture.png");

%% Refocusing and focal-stack generation (40 points)

average_img = zeros(s,t,c);
for i=1:u
    for j=1:v
        %average_img(:,:,:) = average_img(:,:,:) + SB_view( (s*(i-1))+1:(s*(i)), ...
        %                                                   (t*(j-1))+1:(t*(j)),:) / (u*v);
        tmp = zeros(s,t,c);
        tmp(:) = LF(i,j,:,:,:);
        average_img(:,:,:) = average_img(:,:,:) + tmp / (u*v);
    end
end
% imshow(average_img);

depth = [-0.5, 0, 0.5, 1, 1.5, 2];
f_img = zeros(s,t,c);
Focal_Stack = zeros(s,t,c,size(depth,2));

count = 0;
for d=1:size(depth,2)
    for i=1:s
        for j=1:t
            for k=1:c
                for x=1:u
                    for y=1:v
                        tmp_i = i + round((x-u/2)*depth(d));
                        tmp_j = j - round((y-v/2)*depth(d));
                        if tmp_i < 1 || s < tmp_i
                            continue
                        end
                        if tmp_j < 1 || t < tmp_j
                            continue
                        end
                        f_img(i,j,k) = f_img(i,j,k) + LF(x,y, tmp_i, tmp_j,k);
                        count = count+1;
                    end
                end
                f_img(i,j,k) = f_img(i,j,k)/count;
                count = 0;
            end
        end
    end
    Focal_Stack(:,:,:,d) = f_img;
end

imshow(f_img);
imwrite(Focal_Stack(:,:,:,1),"Focal_Stack1.png");
imwrite(Focal_Stack(:,:,:,2),"Focal_Stack2.png");
imwrite(Focal_Stack(:,:,:,3),"Focal_Stack3.png");
imwrite(Focal_Stack(:,:,:,4),"Focal_Stack4.png");
imwrite(Focal_Stack(:,:,:,5),"Focal_Stack5.png");
imwrite(Focal_Stack(:,:,:,6),"Focal_Stack6.png");

%% All-focus image and depth from defocus (35 points). 

sigma1=7;
sigma2=5;
lambda=0.2;

luminance = zeros(s,t,size(depth,2));
low = zeros(s,t,size(depth,2));
high = zeros(s,t,size(depth,2));
sharp = zeros(s,t,size(depth,2));

for d=1:size(depth,2)
    tmp = rgb2xyz(Focal_Stack(:,:,:,d));
    luminance(:,:,d) = tmp(:,:,2);
    
    low(:,:,d) = imgaussfilt(luminance(:,:,d), sigma1);
    high(:,:,d) = luminance(:,:,d) - low(:,:,d);
    sharp(:,:,d) = imgaussfilt(high(:,:,d) .^ 2 , sigma2);
end

all_focus = zeros(s,t,c);
f_depth = zeros(s,t);

for i=1:s
    for j=1:t
        d_num = 0;
        d_den = 0;
        for d=1:size(depth,2)
           d_num = d_num + sharp(i,j,d)*d*lambda;
           d_den = d_den + sharp(i,j,d);
        end
        f_depth(i,j) = d_num/d_den;
        
        for k=1:c
            af_num = 0;
            af_den = 0;
            for d=1:size(depth,2)
               af_num = af_num + sharp(i,j,d)*Focal_Stack(i,j,k,d);
               af_den = af_den + sharp(i,j,d);
            end
            all_focus(i,j,k)=af_num/af_den;
        end
    end
end

imshow(all_focus);
imwrite(all_focus,"all_focus.png");
imshow(1-f_depth);
imwrite(all_focus,"f_depth.png");

