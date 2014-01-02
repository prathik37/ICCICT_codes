%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Extraction using DCT Fusion based on
% Facial Symmetry for Enhanced Face Recognition
%
% K. Manikantan, Prathik P., Rahul Ajay Nafde
%
% Database : Color FERET
% Training : Testing ratio per subject : 1:1
% Mentioned result : 99.60%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% Initialisation of variables

%Global variables, classes is no. of classes in ORL database,
%Imil is mean of each of the 40 classes, Imol is the grand mean
global classes Imil Imol

%DCT subset matrix selection
m_row=12;
m_col=6;

%No. of classes in Color FERET database
classes=80;

%No. of training images per subject
training=1;

%Total no. of images in each subject
total=2;

%No. of trials
trials=10;

%Parameter in logarithm transform
c=1.5;

%Resize factor ( 1 since images have been resized while creating
%                altered database)
rsz=1;

%% Creation of feature face gallery and feature test gallery

for tr=1:trials
    a1=[1:total];
    b1=zeros(classes,training);
    
    %random number generation 
    for i=1:classes
        b1(i,1:training)=randperm(total,training);
        c1(i,1:(total-training))=setdiff(a1,b1(i,1:training));        
    end
    
    %reading an image, splitting into two, applying log transform, using 
    %Laplacian of Gaussian (LoG) edge detection, taking dct, dct subset 
    %selection,conversion into row vectors and storage in feature face gallery
    
    for i=1:classes
        for j=1:training
            img = imread(strcat(pwd,'\s',num2str(i),'\',num2str(b1(i,j)),'.ppm'));
            img = imresize(img,rsz);
            img = rgb2gray(img);
            
            [s_x s_y]= size(img);
            
            img1 = img(1:s_x,1:ceil(s_y/2));
            img2 = img(1:s_x,ceil((s_y/2))+1:s_y);
            
            
            img1=c*log( 1+ double(img1));
            img1=im2uint8(mat2gray(img1));
            img1 = edge(img1,'log');
            
            img1 = dct2(img1);
            img1 = img1(1:m_row,1:m_col);
            img1 = reshape(img1.',1,[]);
            
            img2=c*log( 1+ double(img2));
            img2=im2uint8(mat2gray(img2));
            img2 = edge(img2,'log');
            
            img2 = dct2(img2);
            img2 = img2(1:m_row,1:m_col);
            img2 = reshape(img2.',1,[]);
            
            idctcomp{i,j}=[img1 img2];
                 
        end
    end
    
    %reading an image, splitting into two, applying log transform, using 
    %Laplacian of Gaussian (LoG) edge detection, taking dct, dct subset 
    %selection,conversion into row vectors and storage in feature test gallery
    for i=1:classes
        for j=1:total-training
            img = imread(strcat(pwd,'\s',num2str(i),'\',num2str(c1(i,j)),'.ppm'));
            img = imresize(img,rsz);
            img = rgb2gray(img);
            
            [s_x s_y]= size(img);
            
            img1 = img(1:s_x,1:ceil(s_y/2));
            img2 = img(1:s_x,ceil((s_y/2))+1:s_y);
            
            
            img1=c*log( 1+ double(img1));
            img1=im2uint8(mat2gray(img1));
            img1 = edge(img1,'log');
            
            img1 = dct2(img1);
            img1 = img1(1:m_row,1:m_col);
            img1 = reshape(img1.',1,[]);
            
            img2=c*log( 1+ double(img2));
            img2=im2uint8(mat2gray(img2));
            img2 = edge(img2,'log');
            
            img2 = dct2(img2);
            img2 = img2(1:m_row,1:m_col);
            img2 = reshape(img2.',1,[]);
            
            tdctcomp{i,j}=[img1 img2];
             
        end
    end
    
    %To find means of number of classes
    for i=1:classes
        Imi{i}=0;
        for j=1:training;
            Imi{i}= (Imi{i} + idctcomp{i,j})./training;
        end
        Imil{i}=reshape(Imi{i}.',1,[]);                                    
    end
    
    %to find grand mean of all the classes
    Imo=0;
    for i=1:classes
        Imo=Imo+Imi{1,i};
    end
    Imo=Imo./classes;
    Imol=reshape(Imo.',1,[]);                                              
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Feature selection using BPSO
    
    %BPSO parameters
    dim = 2*m_row*m_col;
    iterations = 30;
    inertia = 0.6;
    correction_factor = 2;
    N = 30;
    
    
    %BPSO initialisation
    for i = 1:N
        for d = 1:dim
            x(i,d)=randi([0,1]);
            vel(i,d) = 0;
            p_best(i,d) = x(i,d);
        end
        
        fit_pbest(i) = bpso_fit(x(i,:));                                      
    end
    
    [fit_gbest,index] = max(fit_pbest);                                      
    g_best = p_best(index,:);
    
    
    %Procedure to find optimum p_best and g_best
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
    
    %face feature gallery construction (choosing selected features)
    for i=1:classes
        for j=1:training
            idctcomp{i,j}= idctcomp{i,j}.*g_best;
        end
    end
    
     %test feature gallery construction (choosing selected features)
    for i=1:classes
        for j=1:total-training
            tdctcomp{i,j}= tdctcomp{i,j}.*g_best;
        end
    end
    
     %% To find recognition rate
    
    %classifier (euclidean distance)
    s=1;
    for i=1:classes
        for j=1:total-training
            
            a=tdctcomp{i,j};
            
            for k=1:classes
                for l=1:training
                    b=idctcomp{k,l};
                    
                    f=(a-b).^2;
                    yy(s)=sqrt(sum(f(:)));
                    s=s+1;
                end
            end
        end
    end
    
    final_values = vec2mat(yy,classes*training);                                          
    for i=1:classes*(total-training)
        [value(i),index(i)] = min(final_values(i,:));
        norm_index(i)=ceil(index(i)/training);                                            
    end
    
    norm_mat=vec2mat(norm_index,total-training);                                
    count=0;
    for i=1:classes
        for j=1:total-training
            if (norm_mat(i,j)==i)                                              
                count=count+1;
            end
        end
    end
    
    
    
    recognition(tr) = count/(classes*(total-training));
    disp('recognition rate');
    disp(recognition(tr)*100);                                                 
    
    
    
    
end

%% To compute final recognition rate
disp('Average recognition rate');
disp((sum(recognition)/trials)*100);
