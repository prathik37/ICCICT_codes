%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Extraction using DCT Fusion based on
% Facial Symmetry for Enhanced Face Recognition
%
% K. Manikantan, Prathik P., Rahul Ajay Nafde
%
% Database : Extended Yale B
% Training : Testing ratio per subject : 3:16
% Mentioned result : 99.04%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

%% Initialisation of variables

cimg=19;                 %no. of images in a class
training=3;              %no. of images chosen for testing
g=1;                     %parameter to store RR for each trial
                         
trials=10;                %no. of trials

%dimensions for fusion
dctRow=18;               %choose DCT size. Specify one only either dctRow
                         %dctCol. Other calculated based on aspect ratio.
dctCol=ceil(2*dctRow/3);

p=2*dctRow*dctCol;       %dimension of 1-Dimensional fused DCT matrix.

global Imil Imol classes %Imil is mean of each of the 28 classes,
                         %Imol is the grand mean.
classes=28;


%% Training of images

for x=1:trials
    disp('Trial no.');
    disp(x);
    a1=[1:cimg];
    
    %random number generation
    for i=1:classes
        b1(i,1:training)=randperm(cimg,training);
        c1(i,1:(cimg-training))=setdiff(a1,b1(i,1:training));
    end
    
    %reading an image, splitting into two, taking dct, dct subset selection
    %conversion into row vectors and storage in feature face gallery
    for i=1:classes
        for j=1:training
            Img= imread(strcat(pwd,'\s',num2str(i),'\',num2str(b1(i,j)),'.pgm'));
            [imgDct]=bisectImg(Img,dctRow,dctCol);
            Idctcomp{i,j}=imgDct;
        end
    end
    
    %to find mean of 28 classes
    for i=1:classes
        Imi{i}=0;
        for j=1:training;
            Imi{i}=Imi{i}+(Idctcomp{i,j})./training;
        end
        Imil{i}=Imi{i};        %to convert mean of all classes to
                               %1-Dimensional arrays
    end
    
    %to find grand mean of all 28 classes
    Imo=zeros(1,p);
    for i=1:classes
        Imo=Imo+Imi{1,i};
    end
    Imo=Imo./classes;
    Imol=Imo;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Feature selection using BPSO
    
    %BPSO parameters
    dim = p;
    iterations = 100;
    inertia = 0.6;
    correction_factor = 2;
    N = 30;
    
    %BPSO initialisation
    for i = 1:N
        for d = 1:dim
            x(i,d) = (rand > 0.5);
            vel(i,d) = 0;
            p_best(i,d) = x(i,d);
        end
        
        fit_pbest(i) = bpso_fit(x(i,:));
    end
    [fit_gbest,index] = max(fit_pbest);
    g_best = p_best(index,:);
    
    
    %To find optimum p_best and g_best
    for iter = 1:iterations
        for i=1:N
            f = bpso_fit(x(i,:));
            
            if (f > fit_pbest(i))
                fit_pbest(i) = f;
                p_best(i,:) = x(i,:);
            end
            
            [fit_gbest,index] = max(fit_pbest);
            g_best = p_best(index,:);
            
            for d=1:dim
                vel(i,d) = (inertia)*(vel(i,d)) + correction_factor * (p_best(i,d) - x(i,d))* rand + correction_factor * (g_best(1,d) - x(i,d))* rand;
                x(i,d) = (rand) < (1/(1 + exp(-vel(i,d))));
                
            end
        end
    end
    
    
    %feature gallery construction (choosing selected features)
    for i=1:classes
        for j=1:training
            zidct{i,j}=Idctcomp{i,j}.*g_best;
        end
    end
    
    %% Testing of images
    
    %reading an image, splitting into two, taking dct, dct subset selection
    %conversion into row vectors and storage in testing gallery
    for i=1:classes
        for j=1:(cimg-training)
            Img = imread(strcat(pwd,'\s',num2str(i),'\',num2str(c1(i,j)),'.pgm'));
            [imgDct]=bisectImg(Img,dctRow,dctCol);
            Tdctcomp{i,j}=imgDct;
        end
    end
    
    
    %test gallery construction (choosing selected features)
    for i=1:classes
        for j=1:(cimg-training)
            Tidct{i,j}=Tdctcomp{i,j}.*g_best;
        end
    end
    
    %% To find Recognition Rate
    
    %classifier (euclidean distance)
    s=1;
    for i=1:classes
        for j=1:(cimg-training)
            
            
            a=Tidct{i,j}; %testing images
            
            for k=1:classes
                for l=1:training
                    b=zidct{k,l}; %trained images
                    
                    f=(a-b).^2;
                    yy(s)=sqrt(sum(f(:)));
                    s=s+1;
                    
                    
                end
            end
        end
    end
    
    final_values = vec2mat(yy,training*classes);
    for i=1:((cimg-training)*classes)
        [value(i),index(i)] = min(final_values(i,:));
        norm_index(i)=ceil(index(i)/training);
    end
    
    norm_mat=vec2mat(norm_index,(cimg-training));
    count=0;
    for i=1:classes
        for j=1:(cimg-training)
            if (norm_mat(i,j)==i)
                count=count+1;
            end
        end
    end
    
    recognition(g) = count/((cimg-training)*classes);
    disp('recognition rate');
    disp(recognition(g)*100);              %recognition rate for each case
    g=g+1;
end


% to compute final recognition rate
v=0;
for i=1:trials
    v= v+ recognition(i);
end

recognition_rate=(v/trials)*100;
disp('Average Recognition rate for the below no of trials'); disp(trials);
disp(recognition_rate);








