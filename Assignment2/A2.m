%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initials and Color Transformation(5pts)
fprintf("INITIALS AND COLOR TRANSFORMATION (5 PTS)\n");

% Open the video file using VideoReader
vidObj_baby2 = VideoReader('./data/face.mp4');

% Get the video properties
numFrames = vidObj_baby2.NumFrames;
frameRate = vidObj_baby2.FrameRate;

% Preallocate memory for the video frames
video = zeros(vidObj_baby2.Height, vidObj_baby2.Width, 3, numFrames, 'double');

% Loop through each frame and read it into the video matrix
for i = 1:numFrames
    video(:,:,:,i) = double(readFrame(vidObj_baby2)) / 255;
end

% Preallocate memory for the video frames
video_yiq = zeros(vidObj_baby2.Height, vidObj_baby2.Width, 3, numFrames, 'double');

% Convert each frame to the YIQ color space
for i = 1:numFrames
    yiq = rgb2ntsc(video(:,:,:,i));
    video_yiq(:,:,:,i) = yiq;
end

% imshow for test
% figure(1); imshow(video(:,:,:,1));
imwrite(video_yiq(:,:,:,1), "./figure/Initials.jpg");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAPLACIAN PYRAMID (20PTS)
fprintf("LAPLACIAN PYRAMID (20PTS)\n");

pyramid_level = 5;
face_L_pyramid = cell(size(video ,4),pyramid_level + 1);
for i = 1:size(video,4)
    face_L_pyramid(i,:) = genpyramid(video(:,:,:,i),pyramid_level,"laplace");
end
imwrite(face_L_pyramid{1, 1}, "./figure/Laplacian_Pyramid1.jpg");
imwrite(face_L_pyramid{1, 2}, "./figure/Laplacian_Pyramid2.jpg");
imwrite(face_L_pyramid{1, 3}, "./figure/Laplacian_Pyramid3.jpg");
imwrite(face_L_pyramid{1, 4}, "./figure/Laplacian_Pyramid4.jpg");
imwrite(face_L_pyramid{1, 5}, "./figure/Laplacian_Pyramid5.jpg");
imwrite(face_L_pyramid{1, 6}, "./figure/Laplacian_Pyramid6.jpg");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMPORAL FILTERING (30PTS)
fprintf("TEMPORAL FILTERING (30PTS)\n");

face_1 = cell(size(face_L_pyramid{1,1}));
face_2 = cell(size(face_L_pyramid{1,2}));
face_3 = cell(size(face_L_pyramid{1,3}));
face_4 = cell(size(face_L_pyramid{1,4}));
face_5 = cell(size(face_L_pyramid{1,5}));
face_6 = cell(size(face_L_pyramid{1,6}));

h = fdesign.bandpass('N,F3dB1,F3dB2',256,0.83,1,30);
Hd = design(h, 'butter');
fftHd = freqz(Hd,301);

% nasty code because of lack of memory during code 

for i = 1:size(face_1,1)
    for j = 1:size(face_1,2)
        for k = 1:size(face_1,3)
            face_1{i,j,k} = Filtering(face_L_pyramid,1,i,j,k,fftHd);
        end
    end
end

for i = 1:size(face_2,1)
    for j = 1:size(face_2,2)
        for k = 1:size(face_2,3)
            face_2{i,j,k} = Filtering(face_L_pyramid,2,i,j,k,fftHd);
        end
    end
end

for i = 1:size(face_3,1)
    for j = 1:size(face_3,2)
        for k = 1:size(face_3,3)
            face_3{i,j,k} = Filtering(face_L_pyramid,3,i,j,k,fftHd);
        end
    end
end

for i = 1:size(face_4,1)
    for j = 1:size(face_4,2)
        for k = 1:size(face_4,3)
            face_4{i,j,k} = Filtering(face_L_pyramid,4,i,j,k,fftHd);
        end
    end
end

face_1 = convert(face_1);
face_2 = convert(face_2);
face_3 = convert(face_3);
face_4 = convert(face_4);


fprintf("DONE!\n");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions

function [pyramid] = genpyramid(img,level,type)
% level <= 5 

G_pyramid = cell(1,level+1);
G_pyramid{1} = img;

for i = 2:level+1
     G_pyramid{i} = G_filter(G_pyramid{i-1});
end

L_pyramid = cell(1,level + 1);

for i = 1:level
    L_pyramid{i} = L_filter(G_pyramid{i},G_pyramid{i+1});
end

L_pyramid{level+1} = G_pyramid{level + 1};

if strcmp(type,"gauss")
    pyramid = G_pyramid{5};
elseif strcmp(type,"laplace")
    pyramid = L_pyramid;
end

end

function [filter_img] = G_filter(img)

g_filter = [1 2 1; 2 4 2; 1 2 1]/16;
temp_img = imfilter(img,g_filter,"replicate","same");
filter_img = zeros(ceil(size(temp_img,1)/2),ceil(size(temp_img,2)/2),size(temp_img,3));
for j = 1:3
    ttemp_img = temp_img(:,:,j);
    filter_img(:,:,j) = ttemp_img(1:2:size(temp_img,1),1:2:size(temp_img,2));
end

end

function [Lap_img] = L_filter(ori,blur)
    Lap_img = zeros(size(ori));
    for j = 1:3
        Lap_img(:,:,j) = ori(:,:,j) - imresize(blur(:,:,j),size(ori(:,:,j)));
    end
end

function [Signal] = Filtering(pyramid,level,i,j,k,filter)

Signal = zeros(size(pyramid,1),1);
    for t = 1:size(pyramid,1)
        Signal(t) = pyramid{t,level}(i,j,k);
    end      
    
Signal = fftshift(fft(Signal));
Signal = real(ifft(ifftshift(Signal.*filter)));

end

function [result] = convert(input)
    result = cell(length(input{1,1,1}),1);
    for t = 1:length(input{1,1,1})
        result{t} = zeros(size(input));
        for i = 1:size(input,1)
            for j = 1:size(input,2)
                for k = 1:size(input,3)
                    result{t}(i,j,k) = input{i,j,k}(t);
                end
            end
        end
    end
end
