face = VideoReader('face.mp4');
baby2 = VideoReader('baby2.mp4');

frame_face = read(face,[1 inf]);
frame_baby2 = read(baby2,[1 inf]);

%% initials and Transformation
fprintf("====");

frame_face = convertframedouble(frame_face);
YIQframe_face = convertframeYIQ(frame_face);

frame_baby2 = convertframedouble(frame_baby2);
YIQframe_baby2 = convertframeYIQ(frame_baby2);

%% Laplacian pyramid Generation
fprintf("====");
pyramid_level = 4;
face_L_pyramid = cell(size(YIQframe_face,4),pyramid_level);
for i = 1:size(YIQframe_face,4)
    face_L_pyramid(i,:) = genpyramid(YIQframe_face(:,:,:,i),pyramid_level,"laplace");
end
baby2_L_pyramid = cell(size(YIQframe_baby2,4),pyramid_level);
for i = 1:size(YIQframe_baby2,4)
    baby2_L_pyramid(i,:) = genpyramid(YIQframe_baby2(:,:,:,i),pyramid_level,"laplace");
end

%% Temporal Filtering
fprintf("====");

face_1 = cell(size(face_L_pyramid{1,1}));
face_2 = cell(size(face_L_pyramid{1,2}));
face_3 = cell(size(face_L_pyramid{1,3}));
face_4 = cell(size(face_L_pyramid{1,4}));

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

baby2_1 = cell(size(baby2_L_pyramid{1,1}));
baby2_2 = cell(size(baby2_L_pyramid{1,2}));
baby2_3 = cell(size(baby2_L_pyramid{1,3}));
baby2_4 = cell(size(baby2_L_pyramid{1,4}));

h = fdesign.bandpass('N,F3dB1,F3dB2',256,2.33,2.67,30);
Hd = design(h, 'butter');
fftHd = freqz(Hd,900);

% nasty code because of lack of memory during code 

for i = 1:size(baby2_1,1)
    for j = 1:size(baby2_1,2)
        for k = 1:size(baby2_1,3)
            baby2_1{i,j,k} = Filtering(baby2_L_pyramid,1,i,j,k,fftHd);
        end
    end
end

for i = 1:size(baby2_2,1)
    for j = 1:size(baby2_2,2)
        for k = 1:size(baby2_2,3)
            baby2_2{i,j,k} = Filtering(baby2_L_pyramid,2,i,j,k,fftHd);
        end
    end
end

for i = 1:size(baby2_3,1)
    for j = 1:size(baby2_3,2)
        for k = 1:size(baby2_3,3)
            baby2_3{i,j,k} = Filtering(baby2_L_pyramid,3,i,j,k,fftHd);
        end
    end
end

for i = 1:size(baby2_4,1)
    for j = 1:size(baby2_4,2)
        for k = 1:size(baby2_4,3)
            baby2_4{i,j,k} = Filtering(baby2_L_pyramid,4,i,j,k,fftHd);
        end
    end
end

baby2_1 = convert(baby2_1);
baby2_2 = convert(baby2_2);
baby2_3 = convert(baby2_3);
baby2_4 = convert(baby2_4);


%% IMAGE RECONSTRUCTION
fprintf("====");

face_G_ori = cell(size(YIQframe_face,4),1);
for i = 1:size(YIQframe_face,4)
    face_G_ori{i} = genpyramid(YIQframe_face(:,:,:,i),pyramid_level,"gauss");
end

recon = zeros(size(frame_face));
g_filter = [1 2 1; 2 4 2; 1 2 1]/16;
for i = 1:size(frame_face,4)
    temp = imresize(face_G_ori{i},2);
    temp = imfilter(temp,g_filter,"replicate","same")+face_4{i};
    temp = imresize(temp,2);
    temp = imfilter(temp,g_filter,"replicate","same")+face_3{i};
    temp = imresize(temp,2);
    temp = imfilter(temp,g_filter,"replicate","same")+face_2{i};
    temp = imresize(temp,2);
    temp = imfilter(temp,g_filter,"replicate","same")+face_1{i};
end

%% functions 

function [double] = convertframedouble(frame)

double = zeros(size(frame));
for i = 1:size(frame,4)
    double(:,:,:,i) = im2double(frame(:,:,:,i));
end

end


function [YIQ] = convertframeYIQ(frame)

YIQ = zeros(size(frame));
for i = 1:size(frame,4)
    YIQ(:,:,:,i) = rgb2ntsc(frame(:,:,:,i));
end

end

function [pyramid] = genpyramid(img,level,type)
% level <= 5 

G_pyramid = cell(1,level+1);
G_pyramid{1} = img;

for i = 2:level+1
     G_pyramid{i} = G_filter(G_pyramid{i-1});
end

L_pyramid = cell(1,level);

for i = 1:level
    L_pyramid{i} = L_filter(G_pyramid{i},G_pyramid{i+1});
end

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
