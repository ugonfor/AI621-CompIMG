function [im_blend] = poissonBlend(im_s, mask_s, im_background)

disp(["source size:" size(im_s)]);
disp(["mask size:" size(mask_s)]);
disp(["background size:" size(im_background)]);

figure(10);
subplot(1,3,1); imshow(im_s); title("source"); 
subplot(1,3,2); imshow(mask_s); title("mask");
subplot(1,3,3); imshow(im_background); title("background");

[N,M,C] = size(im_s);

index = zeros(N, M);
index(1:N*M) = 1: N*M;
v_img = zeros(N, M, C);

masks = sum(sum(mask_s));

for c=1:C
    
    A = sparse(8*masks, N*M);
    b = zeros(8*masks, 1);

    row = 1;
    for n=1:N
        for m=1:M
            if mask_s(n,m) == 1
                A(row, index(n,m)) = 1;
                A(row, index(n,m+1)) = -1;
                b(row) = im_s(n,m,c) - im_s(n,m+1,c);
                row = row + 1;
                A(row, index(n,m)) = 1;
                A(row, index(n,m-1)) = -1;
                b(row) = im_s(n,m,c) - im_s(n,m-1,c);
                row = row + 1;
                A(row, index(n,m)) = 1;
                A(row, index(n+1,m)) = -1;
                b(row) = im_s(n,m,c) - im_s(n+1,m,c);
                row = row + 1;
                A(row, index(n,m)) = 1;
                A(row, index(n-1,m)) = -1;
                b(row) = im_s(n,m,c) - im_s(n-1,m,c);
                row = row + 1;
            end
        end
    end

    for n=1:N
        for m=1:M
            if mask_s(n,m) == 1
                if mask_s(n,m + 1) == 0
                    A(row, index(n,m)) = 1;
                    b(row) = im_s(n,m,c) - im_s(n,m+1,c) + im_background(n, m+1,c);
                    row = row + 1;
                end
                if mask_s(n,m - 1) == 0
                    A(row, index(n,m)) = 1;
                    b(row) = im_s(n,m,c) - im_s(n,m-1,c) + im_background(n, m-1,c);
                    row = row + 1;
                end
                if mask_s(n+1,m) == 0
                    A(row, index(n,m)) = 1;
                    b(row) = im_s(n,m,c) - im_s(n+1,m,c) + im_background(n+1, m,c);
                    row = row + 1;
                end
                if mask_s(n-1, m) == 0
                    A(row, index(n,m)) = 1;
                    b(row) = im_s(n,m,c) - im_s(n-1,m,c) + im_background(n-1, m,c);
                    row = row + 1;
                end
            end
        end
    end
    disp(row);   
    disp(size(A));
    disp(size(b));
    v = A\b;
    v_img_c = reshape(v, [N,M]);
    v_img(:,:,c) = v_img_c;
end

figure(11); imshow(v_img); title("POISSON Blend");
im_blend = v_img;
im_blend(im_blend == 0) = im_background(im_blend==0);