%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To scale all the normalised images into a standard size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

% Define number of classes and images in a subject
classes = 80;
total = 2;

%define minimum value of rows and columns which will be further altered
rmin=100;
cmin=100;

%to find the least number of rows and columns
for i=1:classes
    for j=1:total
        
        img = rgb2gray(imread(strcat(pwd,'\s',num2str(i),'\',num2str(j),'.ppm')));
        [r c]= size(img);
        
        if(r<rmin)
            rmin=r;
        end
        
        if(c<cmin)
            cmin=c;
        end
        
    end
end


%writing the image into a desired location by scaling appropriately
for i=1:classes
    for j=1:total
        img = rgb2gray(imread(strcat(pwd,'\s',num2str(i),'\',num2str(j),'.ppm')));
        scaledImg=imresize(img,[rmin cmin],'bicubic');
        img=['H:\ICCICT 2012 K Manikantan, Prathik P, Rahul Ajay Nafde\databases\fa_fb_scaled\s',num2str(i),'\',num2str(j),'.ppm'];
        imwrite(scaledImg,img,'ppm');
        
    end
end