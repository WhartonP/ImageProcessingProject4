%Main script to read the image in and filter out the periodic signal from
%the image and then correct the non-uniform illumination 
clear; close all;

% fileName = input("Enter the name of the image with extention ex. Image.tif -> ", 's');
imOrig = imread("Proj4.tif");

[nx, ny] = size(imOrig);

FT = log(1 + abs(fftshift(fft2(double(imOrig)))));

FT2 = fft2(double(imOrig));

%for a butterworth filter
% Initialize filter.
% filter1 = ones(2*nx-1,2*ny-1);
% filter2 = ones(2*nx-1,2*ny-1);
% filter3 = ones(2*nx-1,2*ny-1);
d1 = 6;
d0 = 4;
n = 30;

for i = 1:nx
    for j = 1:ny
        dist = ((i-(nx/2))^2 + (j-(ny/2))^2)^.5;
        % Create Butterworth filter.
        filter1(i,j)= 1/(1 + (dist/d1)^(2*n));
%         filter2(i,j) = 1/(1 + (dist/d0)^(2*n));
        filter3(i,j)= 1.0 - filter1(i,j);
%         filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end

imSharpen = FT2 .* filter3;
g = real(ifft(imSharpen));
imshow(g);
% 
% finalImage = log(1 + abs(ifftshift(ifft2(double(imSharpen)))));
% 
% imshow(finalImage);