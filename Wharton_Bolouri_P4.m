%Main script to read the image in and filter out the periodic signal from
%the image and then correct the non-uniform illumination 
%By Peter Wharton & Farshad Bolouri
clear; close all;

fileName = input("Enter the name of the image with extention ex. Image.tif -> ", 's');
imOrig = imread(fileName);
% figure; imshow(imOrig);


%For extracting the pattern
FT = (fftshift(fft2((imOrig))));
ave = mean(mean(FT));
% imagesc(log(1+abs(FT)))
% colormap("gray");
% colorbar()


% Finding Spikes
temp = FT;
spikes = FT(isoutlier(max(log(1+abs(FT)),[],2)),:);
spikes = (max(log(1+abs(spikes)),[],2));

indices = zeros(length(spikes),2);
i = 1;
while i < length(spikes) + 1
    [row, col] = find(temp == max(max(temp)));
    if i == 1
        maxRow= row;
        maxCol = col;
    end
    if isempty((find(row == indices(:,1))) == 0) &&...
            isempty((find(col == indices(:,2))) == 0)
        indices(i,:) = [row, col];
        i = i + 1;
    end
    temp(row, col) = 0;
end

% Extracting Pattern
flag = true;
for i = 1:size(FT,1)
    for j = 1:size(FT,2)
        for k = 1:size(indices,1)
            if FT(i,j) == FT(indices(k,1), indices(k,2))
                flag =false;
                if FT(i,j) ~= FT(maxRow, maxCol)
                    FT(i,j) = FT(i,j) * 10;
                end
                break;
            end
        end
        if flag == true
            FT(i,j) = 0;
        end
        flag = true;
    end
end

% figure
% imagesc(log(1+abs(FT)))
% colorbar()
% colormap("gray")

finalPatternImg = fftshift(ifft2(FT));
finalPatternImg = imadjust(uint8(abs(finalPatternImg)));

%For correcting the non-uniform lighting

lightValue = mean(mean(imOrig));
se = offsetstrel('ball',120,16,8);
background = imopen(imOrig, se);
imMinusBack = imOrig - background;
finalBrightImg = imsharpen(imMinusBack + uint8(lightValue), 'Amount', .8);

subplot(2, 2, 1); imshow(imOrig);
title('Original Input Image');
subplot(2, 2, 2); imshow(finalPatternImg);
title('Extracted Pattern');
subplot(2 ,2, 3); imshow(finalBrightImg);
title('Corrected non-uniform Lighting');

