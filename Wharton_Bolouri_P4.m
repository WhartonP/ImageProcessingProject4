%Main script to read the image in and filter out the periodic signal from
%the image and then correct the non-uniform illumination 
clear; close all;

% fileName = input("Enter the name of the image with extention ex. Image.tif -> ", 's');
imOrig = imread("Proj4.tif");
figure; imshow(imOrig);

FT = fftshift(fft2(double(imOrig)));
FTView = (1 + abs(FT));
figure; plot(FTView);

%for a butterworth filter
% Initialize filter.
% filter1 = ones(2*nx-1,2*ny-1);
% filter2 = ones(2*nx-1,2*ny-1);
% filter3 = ones(2*nx-1,2*ny-1);
[nx, ny] = size(imOrig);
d1 = 400;
d0 = 2^8;
n = 2;

for i = 1:nx
    for j = 1:ny
        dist = ((i-(nx/2))^2 + (j-(ny/2))^2)^.5;
        % Create Butterworth filter.
        filter1(i,j)= 1/(1 + (dist/d1)^(2*n));
        filter2(i,j) = 1/(1 + (dist/d0)^(2*n));
        filter3(i,j)= 1.0 - filter2(i,j);
        filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end

imSharpen = filter3 .* FT;
imSharpenView = log(1 + abs(imSharpen));
figure; plot(imSharpenView);

g = real(ifft2(fftshift(imSharpen)));
figure; imshow(g);



