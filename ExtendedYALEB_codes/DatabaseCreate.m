%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is used to pre-process all the images of Extended Yale B and 
% store them in the desired location for training and testing purpose.
% 
% Pre-processing steps carried out are: 
% Raw image -> Resolution reduction -> gamma intensity correction -> 
% -> Local Histogram Equalization 
%
% NOTE: Change current folder to "ExtendedYALEB_subset5" folder and its 
% path. Before executing the program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               INITIALISATION OF PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of standard resized image
% numberRows=60;
% numberColumns=40;

% Number of subjects and images of each subject
totalSubjects=28;
totalImgSubject=19;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      DATABASE CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for subjectNum=1:totalSubjects;
    
    %Specify the directory path where the images are to be stored
    parentDir=strcat('C:\Users\Rahul Nafde\Desktop\ICCICT2012\databases\PreprocessedYaleB_subset5');
    childDir=strcat('s',num2str(subjectNum));
    mkdir(parentDir,childDir);
    
    
    for subjectImgNum=1:totalImgSubject;
        subjectImg=strcat('s',num2str(subjectNum),'/',num2str(subjectImgNum),'.pgm');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       PROCESSING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imgIn=imread(subjectImg);
        
        imgIn=imresize(imgIn,0.125);
        imgIn=imadjust(imgIn,[],[],0.5);
        imgLocalHist=histeqLocal(imgIn);
         
        %%%%%%%%%%%%%%% END OF PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        imgOut=['H:\ICCICT 2012 K Manikantan, Prathik P, Rahul Ajay Nafde\databases\PreprocessedYaleB_subset5\s',num2str(subjectNum),'\',num2str(subjectImgNum),'.pgm'];
        imwrite(imgLocalHist,imgOut,'pgm');
    end
end
%%%%%%%%%%%% DATABASE CREATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Database Created');
timeElapsed=toc/60;
disp('Time Elapsed');
disp(timeElapsed);
