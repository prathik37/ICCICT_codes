%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To create Scale Normalised database using skin detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all
close all

% Define number of classes and images in a subject
classes = 80;
total = 2;

for i=1:classes
    for j=1:total
        
        img = imread(strcat(pwd,'\s',num2str(i),'\',num2str(j),'.ppm'));
        img = imresize(img,0.25);
        [out edgeImg]=generate_skinmap(img);

        grayImg=rgb2gray(img);

        [rows columns]=size(grayImg);
        sizeImg=size(edgeImg);
        
        %To find rmin,rmax,cmin,cmax        
        rowSums=sum(edgeImg.');
        
        for rowCounter=1:sizeImg(1);
            if (rowSums(rowCounter)~=0)
                rowMin=rowCounter;
                break;
            end
        end
        
        
        for rowCounter=sizeImg(1):-1:1;
            if (rowSums(rowCounter)~=0)
                rowMax=rowCounter;
                break;
            end
        end
        
        colSums=sum(edgeImg);
        
        for colCounter=1:sizeImg(2);
            if (colSums(colCounter)~=0)
                colMin=colCounter;
                break;
            end
        end
        
        
        for colCounter=sizeImg(2):-1:1;
            if (colSums(colCounter)~=0)
                colMax=colCounter;
                break;
            end
        end
        
        %Finding normalization factor        
        rowNormalFact=rowMax-rowMin;
        colNormalFact=colMax-colMin;
        
        scaledImg=grayImg(rowMin:rowMax,colMin:colMax);
        [x y]=size(scaledImg);
       
        
        %Writing the image to a desired location
        img=['H:\ICCICT 2012 K Manikantan, Prathik P, Rahul Ajay Nafde\databases\fa_fb_scaled\s',num2str(i),'\',num2str(j),'.ppm'];
        imwrite(scaledImg,img,'ppm');
    end
end