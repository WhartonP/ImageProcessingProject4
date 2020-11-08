
%Main script to read the image in and filter out the periodic signal from
%the image and then correct the non-uniform illumination
clear; close all;

% fileName = input("Enter the name of the image with extention ex. Image.tif -> ", 's');
imOrig = (imread("Proj4.tif"));

[nx, ny] = size(imOrig);

FT = (fftshift(fft2((imOrig))));

%FT2 = fft2(double(imOrig));
%FT2 =  (1 + abs(fftshift(fft2((imOrig)))));
figure;
imagesc(log(1+abs(FT)))
colorbar()
colormap("gray")%colorbar()
figure;
plot(abs((FT)))
%% Filter Design
fc = 182;
fp = 228;
fs = 1E12;

% [b,a] = butter(6,[fc/(fs/2) fp/((fs/2))],'bandpass');
% %[b,a] = butter(6,fc/(fs/2),'high');
% FT = filter2([b,a],FT);
% figure;
% imagesc(log(1+abs(FT)))
% colorbar()
% colormap("gray")

hi = bandpassfilter([3,3],fc/(fs/2),fp/((fs/2)),6);
 
figure;
plot(log(1+abs(FT)))
%

%g = real(ifft2(ifftshift((FT))));
g = (fftshift(ifft2(FT)));
figure;
imshowpair((imOrig-uint8(abs(g))),imOrig,'montage');
figure;
imshowpair((uint8(abs(g))),imread("Proj4_pattern.tif"),'montage');
figure;
imshowpair(histeq(uint8(abs(g))),imread("Proj4_pattern.tif"),'montage');

%% Functions
function Y = imfft(X)
Y = fftshift(fft2(X));
end

function Y = imifft(X)
Y = abs(ifft2(ifftshift(X)));
end

function immagphase(X)
subplot(121)
A = abs(X);
if min(min(A)) == 0
    Y = mat2gray(log(A + 1e-7));
else
    Y = mat2gray(log(A));
end
subimage([-1 1],[-1 1],Y)
title 'Magnitude'
subplot(122)
Y = mat2gray(angle(X));
subimage([-1 1],[-1 1],Y)
title 'Phase'
end

% LOWPASSFILTER - Constructs a low-pass butterworth filter.
%
% usage: f = lowpassfilter(sze, cutoff, n)
%
% where: sze    is a two element vector specifying the size of filter
%               to construct.
%        cutoff is the cutoff frequency of the filter 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%
% The frequency origin of the returned filter is at the corners.
%
% See also: HIGHPASSFILTER, HIGHBOOSTFILTER, BANDPASSFILTER
%

% Peter Kovesi   pk@cs.uwa.edu.au
% Department of Computer Science & Software Engineering
% The University of Western Australia
%
% October 1999
%
% Modified Rob Gaddi   gaddi@rice.edu
% ELEC 301
% Rice University
%
% December 2001

function f = lowpassfilter(sze, cutoff, n)

if cutoff < 0 | cutoff > 0.5
    error('cutoff frequency must be between 0 and 0.5');
end

if rem(n,1) ~= 0 | n < 1
    error('n must be an integer >= 1');
end

%  Modification ELEC 301 Project Group, Dec 2001
%  Original code [rows, cols] = sze was not accepted by Matlab
rows = sze(1);
cols = sze(2);
%  End Alteration

% X and Y matrices with ranges normalised to +/- 0.5
x =  (ones(rows,1) * [1:cols]  - (fix(cols/2)+1))/cols;
y =  ([1:rows]' * ones(1,cols) - (fix(rows/2)+1))/rows;

radius = sqrt(x.^2 + y.^2);        % A matrix with every pixel = radius relative to centre.

%  Alteration, ELEC 301 Project Group, Dec 2001
%  Original code fftshifted the filter before output.  Since
%  imFFT and imIFFT already shift, the output should remain low-centered.
f = 1 ./ (1.0 + (radius ./ cutoff).^(2*n));   % The filter
%  End Alteration
end

% BANDPASSFILTER - Constructs a band-pass butterworth filter
%
% usage: f = bandpassfilter(sze, cutin, cutoff, n)
%
% where: sze    is a two element vector specifying the size of filter
%               to construct.
%        cutin and cutoff are the frequencies defining the band pass 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%
% The frequency origin of the returned filter is at the corners.
%
% See also: LOWPASSFILTER, HIGHPASSFILTER, HIGHBOOSTFILTER
%

% Peter Kovesi   pk@cs.uwa.edu.au
% Department of Computer Science & Software Engineering
% The University of Western Australia
%
% October 1999

function f = bandpassfilter(sze, cutin, cutoff, n)

if cutin < 0 | cutin > 0.5 | cutoff < 0 | cutoff > 0.5
    error('frequencies must be between 0 and 0.5');
end

if rem(n,1) ~= 0 | n < 1
    error('n must be an integer >= 1');
end

f = lowpassfilter(sze, cutoff, n) - lowpassfilter(sze, cutin, n);
end

% HIGHPASSFILTER  - Constructs a high-pass butterworth filter.
%
% usage: f = highpassfilter(sze, cutoff, n)
%
% where: sze    is a two element vector specifying the size of filter
%               to construct.
%        cutoff is the cutoff frequency of the filter 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%
% The frequency origin of the returned filter is at the corners.
%
% See also: LOWPASSFILTER, HIGHBOOSTFILTER, BANDPASSFILTER
%

% Peter Kovesi   pk@cs.uwa.edu.au
% Department of Computer Science & Software Engineering
% The University of Western Australia
%
% October 1999

function f = highpassfilter(sze, cutoff, n)

if cutoff < 0 | cutoff > 0.5
    error('cutoff frequency must be between 0 and 0.5');
end

if rem(n,1) ~= 0 | n < 1
    error('n must be an integer >= 1');
end

f = 1.0 - lowpassfilter(sze, cutoff, n);
end

%EFILTER constructs an elliptical lowpass filter Butterworth filter
%   E = EFILTER(Size, cutoffM, cutoffm, n) designs an Nth order elliptical
%   lowpass digital Butterworth filter where Size is a two element
%   [rows, cols] vector specifying the size of the filter to construct,
%   cutoffM and cutoffm, the cutoff freqency on the major and minor
%   axes are 0 < cutoff <= 1.
%
%   If E = EFilter(Size, cutoffM, cutoffm, n, alpha), where alpha is an angle
%   in radians, it will return and plot an elliptical filter rotated
%   counter-clockwise through alpha.
%
%   If E = EFilter(Size, cutoffM, cutoffm, n, alpha, xoff, yoff), where xoff
%   and yoff are offsets in the x and y direction, it will return and
%   plot an eliptical filter which is offset by the specified amount.
%   An offset of 0 corresponds to the center and an offset of 1
%   corresponds to the edge of the filter. A positive offset shifts the
%   filter in the positive direction.
%
%   Calling EFilter(...) without assigning the output variable
%   plots the 3D surface described by the function.

% Katie Streit   kstreit@rice.edu
% ELEC 301
% Rice University
%
% December 2001

% Much of this code was based on Peter Kovesi's  (pk@cs.uwa.edu.au)
% Matlab function for a lowpass Butterworth filter.


function varargout = EFilter(sze, cutoffM, cutoffm, n, varargin);

if nargin == 4
    alpha = 0;
    offx = 0;
    offy = 0;
elseif nargin == 5
    offx = 0;
    offy = 0;
    alpha = varargin{1};
elseif nargin == 7
    alpha = varargin{1};
    offx = varargin{2};
    offy = varargin{3};
else
    error('Invalid number of input arguments');
end

if nargout > 1
    error('Invalid number of output arguments');
end

if cutoffM < 0 | cutoffM > 1
    error('cutoffM frequency must be between 0 and 1');
end

if cutoffm < 0 | cutoffm > 1
    error('cutoffm frequency must be between 0 and 1');
end

if rem(n,1) ~= 0 | n < 1
    error('n must be an integer >= 1');
end


%%extracts the sizes from sze
rows = sze(1);
cols = sze(2);

%x and y matrices normalized to +/-.5 and an offset of offx or offy
x =  ((((ones(rows,1) * [1:cols])-offx*rows/2)  - (fix(cols/2)+1))/cols);
y =  ((([1:rows]' * ones(1,cols))-offy*rows/2) - (fix(rows/2)+1))/rows;

%applies a linear transformation to rotate through alpha. Note that it takes
% uses negative alpha, which is caused by x and y being independent matrices.
x2 = (x*cos(alpha) - y*sin(-alpha));
y2 = (x*sin(-alpha) + y*cos(alpha));

%constructs an elliptical cone (defined by a and b) of height r on each at
%each elliptical ring. (r is effectively the "radius")
%r = sqrt(((x2/a).^2 + (y2/b).^2));

%Designs the filter
%f = 1./(1.0 + (r./cutoff).^(2*n));

a = cutoffM/2;
b = cutoffm/2;

f = 1./(1+((x2/(a)).^2 + (y2/(b)).^2).^n);

if nargout > 0
    varargout{1} = f;
else
    %Plots a normalized (+/- 1), interpolated 3D image of the filter
    surf([-1:2/(cols-1):1],[-1:2/(rows-1):1], f);
    shading interp;
    title('Elliptical Butterworth filter');
    xlabel('x');
    ylabel('y');
    zlabel('intensity');
    grid on;
end
end