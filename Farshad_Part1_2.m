%Main script to read the image in and filter out the periodic signal from
%the image and then correct the non-uniform illumination
clear; close all;

% fileName = input("Enter the name of the image with extention ex. Image.tif -> ", 's');
imOrig = (imread("Proj4.tif"));

[nx, ny] = size(imOrig);

FT = (fftshift(fft2((imOrig))));

%FT2 = fft2(double(imOrig));
%FT2 =  (1 + abs(fftshift(fft2((imOrig)))));

imagesc(log(1+abs(FT)))

colormap("gray");
colorbar()
%colormap("gray")
%k= [];
%m = 1;
ave = mean(mean(FT));
%% Finding Spikes
temp = FT;
spikes = 7;
indices = zeros(spikes,2);
i = 1;
while i < 8
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
%% 

flag = true;
for i = 1:size(FT,1)
    for j = 1:size(FT,2)
        for k = 1:spikes
            if FT(i,j) == FT(indices(k,1), indices(k,2))
                flag =false;
                if FT(i,j) ~= FT(maxRow, maxCol)
                    FT(i,j) = FT(i,j) * 10;
                end
                break;
            end
        end
        if flag == true
            FT(i,j) = ave;
        end
        flag = true;
    end
end


figure
imagesc(log(1+abs(FT)))
colorbar()
colormap("gray")

%

%g = real(ifft2(ifftshift((FT))));
g = (fftshift(ifft2(FT)));
figure;
imshow(imadjust(uint8(abs(g))));