%Main script to read the image in and filter out the periodic signal from
%the image and then correct the non-uniform illumination 
clear; close all;

fileName = input("Enter the name of the image with extention ex. Image.tif -> ", 's');
imOrig = imread(fileName);