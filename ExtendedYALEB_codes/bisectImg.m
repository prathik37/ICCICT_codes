%This function splits the image into two parts, computes DCT coefficient of
%each part and fuses the two DCT's obtained.

function [fusedDct]=bisectImg(imgIn,dctRow,dctCol)


% clc
% clear all
% close all

%imgIn=imread('s1\1.pgm');

% PreProcessed 
% imgIn=imresize(imgIn,0.125);
% imgIn=imadjust(imgIn,[],[],0.5);
 imgIn=histeq(imgIn);


[numRow numCol]=size(imgIn);

%Edge detection
 myFilt=fspecial('sobel');
%imgIn=imfilter(imgIn,myFilt,'same');

%Disect Image
imgPart1=imgIn(1:numRow,1:(numCol/2));
imgPart2=imgIn(1:numRow,(numCol/2)+1:numCol);


% Verify
% figure,
% subplot(121)
% imshow(imgPart1);
% subplot(122)
% imshow(imgPart2);

imgPart1=imfilter(imgPart1,myFilt,'same');
imgPart2=imfilter(imgPart2,myFilt,'same');

%To compute Dct
img1Dct=dct2(imgPart1);
img2Dct=dct2(imgPart2);

%Verify
% figure,
% subplot(121)
% image(img1Dct);
% subplot(122)
% image(img2Dct);


%Choose Dct size
% dctRow=5;
% dctCol=ceil(2*dctRow/3);

ImgDct1=img1Dct(1:dctRow,1:dctCol);
imgDct2=img2Dct(1:dctRow,1:dctCol);

%Verify
% r1=idct2(ImgDct1);
%  r2=idct2(imgDct2);
%  figure,
% subplot(121)
% imshow(mat2gray(r1));
% subplot(122)
%  imshow(mat2gray(r2));


%Fuse Dcts
temp=reshape(ImgDct1.',1,[]);
fusedDct(1:(dctRow*dctCol))=temp;
temp=reshape(imgDct2.',1,[]);
fusedDct((dctRow*dctCol)+1:2*(dctRow*dctCol))=temp;

%No Fusion
% temp=dct2(imgIn);
% temp=reshape(temp.',1,[]);
% fusedDct(1:2*(dctRow*dctCol))=temp(1:2*(dctRow*dctCol));

%Verify
% figure
% imshow(fusedDct);


end