function [histeqImg]=histeqLocal(imgIn)

Img=imgIn; %input image equated to function's local variable Img.
%Computation performed on Img.

%Select window size.
M=10;   %row size.
N=20;   %column size.

mid_val=round((M*N)/2); %position of middle value of the region defined by window size.

%To find the number of rows and columns to be padded with zero.
in=0;
for i=1:M
    for j=1:N
        in=in+1;
        if(in==mid_val)
            PadM=i-1;
            PadN=j-1;
            break;
        end
    end
end

%Padding the image with zero on all sides.
B=padarray(imgIn,[PadM,PadN]);

for i= 1:size(B,1)-((PadM*2)+1)
    for j=1:size(B,2)-((PadN*2)+1)
        cdf=zeros(256,1);    %variable to store cumulative density function (CDF)
                             %for each region.
                             
        inc=1;               %inc -> counter to track middle element in every region.
        
        %Region defined by window size.
        for x=1:M
            for y=1:N
                %To find the middle element in the region.
                if(inc==mid_val)
                    ele=B(i+x-1,j+y-1)+1;
                end
                pos=B(i+x-1,j+y-1)+1;
                cdf(pos)=cdf(pos)+1;    %storing the probability density of
                                        %each value for the region under
                                        %consideration.
                                        
                inc=inc+1;              %update counter.
            end
        end
        
        %To compute the cumulative density function (CDF) for the values in the region.
        for l=2:256
            cdf(l)=cdf(l)+cdf(l-1);
        end
        Img(i,j)=round(cdf(ele)/(M*N)*255);    % Maximum Range is 255 for gray image.
    end
end

histeqImg=Img;

%Verify the histograms.
%figure,imshow(Img);
% figure,
% subplot(2,1,1);title('Before Local Histogram Equalization'); imhist(inImg);
% subplot(2,1,2);title('After Local Histogram Equalization'); imhist(histeqImg);

end