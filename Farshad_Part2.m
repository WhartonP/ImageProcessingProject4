clear
close all

figure
I = imread('Proj4.tif');
% imshow(I)

figure; imshow(imread("Proj4_uniform.tif"));
pause()
lightValue = mean(mean(I));
% for i = 1:5:120
% se = strel('disk',15);
% se = strel('sphere',i);
% background = imopen(I,se);
% I2 = I - background;
% I3 = (I2 + uint8(lightValue));  %Lightenning up the image

% figure; imshow(I3);

se2 = offsetstrel('ball',120,16,8);
background2 = imopen(I, se2);
I4 = I - background2;
I5 = I4 + uint8(lightValue);
figure; imshow(I5);
%imshow(background)


%imshow(I2)

% figure
% I3 = (I2 + uint8(lightValue));  %Lightenning up the image
% imshow(I3)
% end

