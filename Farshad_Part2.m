clear
close all

figure
I = imread('Proj4.tif');
imshow(I)

se = strel('disk',15);
background = imopen(I,se);
%imshow(background)

I2 = I - background;
%imshow(I2)

figure
I3 = (I2 + uint8(100));  %Lightenning up the image
imshow(I3)

