%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Extraction using DCT Fusion based on
% Facial Symmetry for Enhanced Face Recognition
%
% K. Manikantan, Prathik P., Rahul Ajay Nafde
%
% Database : ORL
% Training : Testing ratio per subject : 9:1
% Mentioned result : 98.50%
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

%No. of classes in ORL database
classes=40;  

%No. of training images per subject
training=9;

%Total no. of images in each subject
total=10;

%parameter to store recognition rates for each trial
g=1;   

%No. of trials 
trials=10;                                                                                                                              

%% Training of images

for x=1:trials                                                             
    disp('Trial no.');
    disp(x);
    a1=[1:total];
    b1=zeros(classes,training);
    
    
    %random number generation    
    for i=1:classes
        b1(i,1:training)=randperm(total,training);
        c1(i,1:(total-training))=setdiff(a1,b1(i,1:training));
        
    end
    
    
    %reading an image, splitting into two, taking dct, dct subset selection
    %conversion into row vectors and storage in feature face gallery
    for i=1:classes                                                            
        for j=1:training
            
            img = imread(strcat(pwd,'\s',num2str(i),'\',num2str(b1(i,j)),'.pgm'));  
            img = imresize(img,0.5,'bicubic');
            
            [rr rc]= size(img);
            img1 = img(1:rr,1:rc/2);
            img2 = img(1:rr,(rc/2)+1:rc);
            temp =dct2(img1);
            temp = temp(1:m_row,1:m_col);
            Idct1=reshape(temp.',1,[]);
            temp= dct2(img2);
            temp = temp(1:m_row,1:m_col);
            Idct2=reshape(temp.',1,[]);
            Idctcomp{i,j}= [Idct1 Idct2];
                       
        end
    end
    
    %To find means of 40 classes
    for i=1:classes
        Imi{i}=0;
        for j=1:training;
            Imi{i}= (Imi{i} + Idctcomp{i,j})./training;
        end
        Imil{i}=reshape(Imi{i}.',1,[]);                                    
    end
    
    %To find grand mean of all 40 classes
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
        for j=1:total-training
            
            Timg = imread(strcat(pwd,'\s',num2str(i),'\',num2str(c1(i,j)),'.pgm'));
            Timg = imresize(Timg,0.5,'bicubic');
            
            [rt rc] = size(Timg);
            Timg1 = Timg(1:rt,1:rc/2);
            Timg2 =Timg(1:rr,(rc/2)+1:rc);
            temp= dct2(Timg1);
            temp= temp(1:m_row,1:m_col);
            Tdct1=reshape(temp.',1,[]);
            temp= dct2(Timg2);
            temp= temp(1:m_row,1:m_col);
            Tdct2=reshape(temp.',1,[]);
            Tdctcomp{i,j}= [Tdct1 Tdct2];
            
        end
    end
    
    
   %test gallery construction (choosing selected features)
    for i=1:classes
        for j=1:total-training
            Tidct{i,j}=Tdctcomp{i,j}.*g_best;
        end
    end
    
    %% To find recognition rate
    
    %classifier (euclidean distance)
    s=1;
    for i=1:classes
        for j=1:total-training
            
            a=Tidct{i,j};
            
            for k=1:classes
                for l=1:training
                    b=zidct{k,l};
                    
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
    
    recognition(g) = count/(classes*(total-training));
    disp('recognition rate');
    disp(recognition(g)*100);                                                  
    g=g+1;
end


%% to compute final recognition rate
v=0;
for i=1:trials
    v= v+ recognition(i);
end

recognition_rate=(v/trials)*100;
disp('Average Recognition rate');
disp(recognition_rate);









