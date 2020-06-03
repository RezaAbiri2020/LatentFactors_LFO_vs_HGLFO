
% performing Dimensionality Reduction using PCA 

close all;
clear all;

%% General parameters that should be specified at the begining

% choosing brain area to do analysis on that part and do referencing
% options:
%Selected_Chs;  case 1 all electrodes
%Selected_PMChs; case 2 
%Selected_SMChs; case 3
%Selected_PMSMChs; case 4
%Selected_HandChs; case 5
Brain_part=1;

%type of referencing the ECoG
Reference=1; % for median
%Reference=2; % for mean 

% which band to do analysis: 6 bands
FilterBand=6;

% which feature to do analysis
%Feature=1; % pure filtered bands
Feature=2;   % abs hilbert or envelope


%% loading and breaking raw ECoG data into trials
load('..\ECoGData\ECoG_data.mat');

% Final list of bad channels
BadChs=[65 66 128 128+1 128+2 128+16 128+18 128+20 128+24 128+30 ...
    128+31 128+32 128+44 128+62 128+64];

% make channels Ready for median or mean
All_Index=ones(size(lfp,1),1);
All_Index(BadChs')=0;

% only 256 Channels were used for real recording
Nch_Record=256;
Selected_Chs=logical(All_Index(1:Nch_Record,1));

% time markers in ecog
[Nch,N] = size(lfp);
ecog_time = (0:N-1)/Fs;
anin = anin(1,1:N);

% from Daniel function: start and end of movements
[trial_start_time,trial_end_time] = get_trial_times(anin,Fs);

trial_start=trial_start_time*Fs;
trial_stop=trial_end_time*Fs;

% F1_ECoGdata=lfp(1:Nch_Record,trial_start(1):trial_stop(1));
% F2_ECoGdata=lfp(1:Nch_Record,trial_start(2):trial_stop(2));
% F3_ECoGdata=lfp(1:Nch_Record,trial_start(3):trial_stop(3));
% F4_ECoGdata=lfp(1:Nch_Record,trial_start(4):trial_stop(4));
% F5_ECoGdata=lfp(1:Nch_Record,trial_start(5):trial_stop(5));
% G1_ECoGdata=lfp(1:Nch_Record,trial_start(6):trial_stop(6));
% G2_ECoGdata=lfp(1:Nch_Record,trial_start(7):trial_stop(7));
% G3_ECoGdata=lfp(1:Nch_Record,trial_start(8):trial_stop(8));
% Alldata={F1_ECoGdata;F2_ECoGdata;F3_ECoGdata;F4_ECoGdata;F5_ECoGdata;G1_ECoGdata;G2_ECoGdata;G3_ECoGdata};


%% Finding the precentral (for Primary Motor=PM) and postcentral (somatosensory=SM) electrodes/channels

 load('..\ImagingData\TDT_elecs_all.mat')
 
 PM_chs=[];
 SM_chs=[];
 PMSM_chs=[];
 
 for i=65:320
     
     % removing/passing bad channels 
     if ismember(i-64,BadChs)
         continue 
         
      % finding PM channels   
     elseif strcmp(anatomy{i,4},'precentral')
         PM_chs=[PM_chs;(i-64)];
         PMSM_chs=[PMSM_chs;(i-64)];
         
      % finding SM channels   
     elseif strcmp(anatomy{i,4},'postcentral')
         SM_chs=[SM_chs;(i-64)];
         PMSM_chs=[PMSM_chs;(i-64)];
         
     end   
     
     
 end
 

% logic PM channels
All_Index_PM=zeros(size(lfp,1),1);
All_Index_PM(PM_chs)=1;
Selected_PMChs=logical(All_Index_PM(1:Nch_Record,1));
 

% logic SM channels
All_Index_SM=zeros(size(lfp,1),1);
All_Index_SM(SM_chs)=1;
Selected_SMChs=logical(All_Index_SM(1:Nch_Record,1));

% logic PMSM channels
All_Index_PMSM=zeros(size(lfp,1),1);
All_Index_PMSM(PMSM_chs)=1;
Selected_PMSMChs=logical(All_Index_PMSM(1:Nch_Record,1));


%% Finding the hand area electrodes/channels with observation

% hand area channels
for i=0:4
    HandChs(((i*7)+1):((i*7)+7))=(122:128)+i*16;
end 

 All_Index_Hand=zeros(Nch_Record,1);
 for i=1:Nch_Record
     
     % removing/passing bad channels 
     if ismember(i,BadChs)
         continue 
         
      % finding hand channels   
     elseif ismember(i,HandChs)
         All_Index_Hand(i,1)=1;
         

         
     end   
     
     
 end
 
 % logic Hand channels
Selected_HandChs=logical(All_Index_Hand);

%% Determining brain section for analysis

switch Brain_part
    
    case 1
        Brain_Area=Selected_Chs;
    case 2
        Brain_Area=Selected_PMChs;
    case 3
        Brain_Area=Selected_SMChs;
    case 4
        Brain_Area=Selected_PMSMChs;
    case 5
        Brain_Area=Selected_HandChs;
end


 %% All Kinematics data for all movements: loading/preparation 
num_trials =8;

for trial=1:num_trials
    % load trial data
    filename = '../LeapMotionData/20180301/10h02m38s/EC171_20180301_10h02m38s_trial00';
    load([filename,num2str(trial),'.mat']);
   
    % grab kinematics
    [time_trial,hand_trial,fingers_trial] = leap_motion_kinematics(data.leap_frame);
   
    % time stamps for kinematics data
    time_AllTrials{trial} = time_trial; 
    % hand pos
    hand{trial}.palm_pos = hand_trial.palm_pos;
    hand{trial}.palm_vel = hand_trial.palm_vel;
    hand{trial}.palm_normal = hand_trial.palm_normal;
    % finger positions
    for f=1:5
        fingers{trial}(f).pos = fingers_trial(f).pos;
        fingers{trial}(f).vel = fingers_trial(f).vel;
        fingers{trial}(f).dir = fingers_trial(f).dir;
    end
end

% smooth all the trajectories and calc velocities
for trial=1:num_trials
    for f=1:5
        fingers{trial}(f).spos = smooth_kinematics(fingers{trial}(f).pos);
        fingers{trial}(f).svel = smooth_kinematics(fingers{trial}(f).vel);
        fingers{trial}(f).tvel = sqrt(sum(fingers{trial}(f).svel.^2,2)); % tangential
    end
end

%% Fingers movement analysis: upsampling & syncing

% F1 movement analysis: upsampling & syncing
Xq_F1=0:1/(Fs):time_AllTrials{1,1}(end);
Posq_F1=interp1(time_AllTrials{1,1},fingers{1,1}(1).spos,Xq_F1);

% for finger 1 
[F1_n,F1_m]=size(Xq_F1');

%plot(Posq_F1(:,1))

Velq_F1=interp1(time_AllTrials{1,1},fingers{1,1}(1).svel,Xq_F1);
AbsVel_F1=sqrt(Velq_F1(:,1).^2+Velq_F1(:,2).^2+Velq_F1(:,3).^2);

%plot(AbsVel_F1)

% F2 movement analysis: upsampling & syncing
Xq_F2=0:1/(Fs):time_AllTrials{1,2}(end);
Posq_F2=interp1(time_AllTrials{1,2},fingers{1,2}(2).spos,Xq_F2);

% for finger 1 
[F2_n,F2_m]=size(Xq_F2');

% plot(Posq_F2(:,3))

Velq_F2=interp1(time_AllTrials{1,2},fingers{1,2}(2).svel,Xq_F2);
AbsVel_F2=sqrt(Velq_F2(:,1).^2+Velq_F2(:,2).^2+Velq_F2(:,3).^2);

% plot(AbsVel_F2)

% F3 movement analysis: upsampling & syncing

Xq_F3=0:1/(Fs):time_AllTrials{1,3}(end);
Posq_F3=interp1(time_AllTrials{1,3},fingers{1,3}(3).spos,Xq_F3);

% for finger 1 
[F3_n,F3_m]=size(Xq_F3');

% plot(Posq_F3(:,3))
 
Velq_F3=interp1(time_AllTrials{1,3},fingers{1,3}(3).svel,Xq_F3);
AbsVel_F3=sqrt(Velq_F3(:,1).^2+Velq_F3(:,2).^2+Velq_F3(:,3).^2);

% plot(AbsVel_F3)

% F4 movement analysis: upsampling & syncing
Xq_F4=0:1/(Fs):time_AllTrials{1,4}(end);
Posq_F4=interp1(time_AllTrials{1,4},fingers{1,4}(4).spos,Xq_F4);

% for finger 1 
[F4_n,F4_m]=size(Xq_F4');
% plot(Posq_F4(:,3))

Velq_F4=interp1(time_AllTrials{1,4},fingers{1,4}(4).svel,Xq_F4);
AbsVel_F4=sqrt(Velq_F4(:,1).^2+Velq_F4(:,2).^2+Velq_F4(:,3).^2);

% plot(AbsVel_F4)

% F5 movement analysis: upsampling & syncing
Xq_F5=0:1/(Fs):time_AllTrials{1,5}(end);
Posq_F5=interp1(time_AllTrials{1,5},fingers{1,5}(5).spos,Xq_F5);

% for finger 1 
[F5_n,F5_m]=size(Xq_F5');

% plot(Posq_F5(:,3))

Velq_F5=interp1(time_AllTrials{1,5},fingers{1,5}(5).svel,Xq_F5);
AbsVel_F5=sqrt(Velq_F5(:,1).^2+Velq_F5(:,2).^2+Velq_F5(:,3).^2);

%plot(AbsVel_F5)

%% Optional spatial filtering before extracting th features
SpatialFilteringOne=0;

if SpatialFilteringOne
    % grid layout
    Ch_num_1=1:Nch_Record;
    Ch_num_2=reshape(Ch_num_1,[16,16]);
    % Ch_num_3 is the grid layout
    Ch_num_3=rot90(rot90(Ch_num_2));
    
    % modification to sync
    lfp_1=lfp(1:Nch_Record,:);
    %lfp_1(BadChs,:)=nan;
    % z score
    lfp_2=zscore(lfp_1')';
    %lfp_2=lfp_1;
    
    % re-arranging 
    [R,C]=size(Ch_num_3);
    lfp_Spatial_1=zeros(R,C,length(lfp_2));
    for r=1:R
        for c=1:C
            lfp_Spatial_1(r,c,:)=lfp_2(Ch_num_3(r,c),:);
        end
    end
    
    % performing spatial filtering for each sample
    lfp_Spatial_2=zeros(R,C,length(lfp_2));
    for sample=1:length(lfp_2)
        Values_1=squeeze(lfp_Spatial_1(:,:,sample));
        Values_2=medfilt2(Values_1,[3,3],'symmetric');
        lfp_Spatial_2(:,:,sample)= Values_2;
        
    end
    
    % return to original setup
    lfp_3=zeros(size(lfp_2,1),size(lfp_2,2));
    for ch=1:Nch_Record
        [r c]=find(Ch_num_3==ch);
        lfp_3(ch,:)=squeeze(lfp_Spatial_2(r,c,:));
        
    end
       
end

%% Referencing & Filtering options for generating raw bands for all Fingers

Modification=[length(Xq_F1),length(Xq_F2),length(Xq_F3),length(Xq_F4),length(Xq_F5)];

for Fi=1:5
    
    if SpatialFilteringOne
        
        ECoGdata_1=lfp_3(1:Nch_Record,trial_start(Fi):(trial_start(Fi)+Modification(Fi)-1));
        
        %ECoGdata_Modified=lfp_3(1:Nch_Record,trial_start(Fi):(trial_start(Fi)+Modification(Fi)-1));
        %ECoGdata_1=zscore(ECoGdata_Modified')';
        
    else
        % modification to sync
        ECoGdata_Modified=lfp(1:Nch_Record,trial_start(Fi):(trial_start(Fi)+Modification(Fi)-1));
        % z score
        ECoGdata_1=zscore(ECoGdata_Modified')';
        
    end
    
    % option median or mean:
    switch Reference
        case 1
            %if SpatialFilteringOne
             %   ECoGdata_2= ECoGdata_1;
            %else
                F_Median =median(ECoGdata_1(Brain_Area,:),1);
                ECoGdata_2= ECoGdata_1-repmat(F_Median,Nch_Record,1);
            %end
            
        case 2
            %if SpatialFilteringOne
             %   ECoGdata_2= ECoGdata_1;
            %else
                F_Mean =mean(ECoGdata_1(Brain_Area,:),1);
                ECoGdata_2= ECoGdata_1-repmat(F_Mean,Nch_Record,1);
            %end     
    end
    
    % choose the band freq for analysis
    switch FilterBand
        
        case 1 % delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            
        case 2 % theta
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            
        case 3 % alpha
            [b,a]=butter(3,[8, 13]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            
        case 4 % beta
            F_FilteredAll=[];
            F_HilbertAll=[];
            
            Betabands={
                [13,19]
                [19,30]};
            
            for band=1:2
                [b,a]=butter(3,Betabands{band}/(Fs/2));
                F_FilteredAll(:,:,band)=filtfilt(b,a,ECoGdata_2');
                
                % Envelope & Hilbert
                [F_Envelope,Lower]=envelope(F_FilteredAll(:,:,band));
                % Hilbert & Filter
                [bb,aa]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertAll(:,:,band)=filtfilt(bb,aa,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                
            end
            
            F_Filtered=mean(F_FilteredAll,3);
            F_Hilbert=mean(F_HilbertAll,3);
            
            
        case 5 % low gamma
            F_FilteredAll=[];
            F_HilbertAll=[];
            
            LGbands={
                [30,36]
                [36,42]
                [42,50]};
            
            for band=1:3
                [b,a]=butter(3,LGbands{band}/(Fs/2));
                F_FilteredAll(:,:,band)=filtfilt(b,a,ECoGdata_2');
                
                % Envelope & Hilbert
                [F_Envelope,Lower]=envelope(F_FilteredAll(:,:,band));
                % Hilbert & Filter
                [bb,aa]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertAll(:,:,band)=filtfilt(bb,aa,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            end
            
            F_Filtered=mean(F_FilteredAll,3);
            F_Hilbert=mean(F_HilbertAll,3);
            
            
        case 6
            F_FilteredAll=[];
            F_HilbertAll=[];
            
            HGbands={
                [70,77]
                [77,85]
                [85,93]
                [93,102]
                [102,113]
                [113,124]
                [124,136]
                [136,150]};
            
            for band=1:8
                [b,a]=butter(3,HGbands{band}/(Fs/2));
                F_FilteredAll(:,:,band)=filtfilt(b,a,ECoGdata_2');
                
                % Envelope & Hilbert
                [F_Envelope,Lower]=envelope(F_FilteredAll(:,:,band));
                % Hilbert & Filter
                [bb,aa]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertAll(:,:,band)=filtfilt(bb,aa,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                
            end
            
            F_Filtered=mean(F_FilteredAll,3);
            F_Hilbert=mean(F_HilbertAll,3);
                  
  
    end
       
    % pulling out Feature:Filtered data or Feature:Hilberet
    switch Feature
        
        case 1
            FingersFeatures(Fi).Finger=F_Filtered;
            
        case 2
            FingersFeatures(Fi).Finger=F_Hilbert;       
    end 
     
end

%% Optional spatial filtering after extracting the features
SpatialFilteringTwo=0;

if SpatialFilteringTwo
    
        % grid layout
    Ch_num_1=1:Nch_Record;
    Ch_num_2=reshape(Ch_num_1,[16,16]);
    % Ch_num_3 is the grid layout
    Ch_num_3=rot90(rot90(Ch_num_2));
    
    for Fi=1:5
        Data_Finger_1=FingersFeatures(Fi).Finger;
        %Data_Finger_1(:,BadChs)=nan;
        % re-arranging
        [R,C]=size(Ch_num_3);
        Feature_Spatial_1=zeros(R,C,length(Data_Finger_1));
        for r=1:R
            for c=1:C
                Feature_Spatial_1(r,c,:)=Data_Finger_1(:,Ch_num_3(r,c));
            end
        end
        
        % performing spatial filtering for each sample
        Feature_Spatial_2=zeros(R,C,length(Data_Finger_1));
        for sample=1:length(Data_Finger_1)
            Values_1=squeeze(Feature_Spatial_1(:,:,sample));
            Values_2=medfilt2(Values_1);
            Feature_Spatial_2(:,:,sample)= Values_2;
            
        end
        % return to original setup
        Data_Finger_2=zeros(size(Data_Finger_1,1),size(Data_Finger_1,2));
        for ch=1:Nch_Record
            [r c]=find(Ch_num_3==ch);
            Data_Finger_2(:,ch)=squeeze(Feature_Spatial_2(r,c,:));
            
        end
        
      FingersFeatures_Spatial(Fi).Finger=Data_Finger_2; 
        
    end
    
end

%% generating the moving PCA analysis of features
% number of sample point
PCA_Win=4000;
Jump_Win=200;
Num_PCs=10;

Fi_Sizes=[F1_n F2_n F3_n F4_n F5_n];
% PCA for features
for Fi=1:5
    
    if SpatialFilteringTwo
        
        Data_Finger=FingersFeatures_Spatial(Fi).Finger(:,Brain_Area);
    else
        Data_Finger=FingersFeatures(Fi).Finger(:,Brain_Area);
    end
    
        
    F_score_all=[];
    F_explained_all=[];
    
    Data_Finger=[zeros(PCA_Win-Jump_Win+1,size(Data_Finger,2));Data_Finger];
    k=ceil(((size(Data_Finger,1))-PCA_Win)/Jump_Win);
    
    for i=1:k 
        Data_PCA=Data_Finger(Jump_Win*(i-1)+1:(Jump_Win*(i-1))+PCA_Win,:);
        [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data_PCA);
        F_score_all=[F_score_all; F_score(PCA_Win-Jump_Win+1:end,1:Num_PCs)];
        F_explained_all=[F_explained_all, F_explained(1:Num_PCs,:)];   
    end 
    
    Fingers_Score(Fi)={F_score_all};
    Fingers_Explained(Fi)={F_explained_all};
          
end


%% Position prediction & Vel prediction 

Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

 %prediction the position
for Fi=1:5
  
    Projected_PCs=cell2mat(Fingers_Score(Fi));
    Input_dynamics=Projected_PCs;
    
    Output_Dynamics=Kin_Pos{1,Fi}(1:length(Projected_PCs));
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end

% prediction the velocities
for Fi=1:5
  
    Projected_PCs=cell2mat(Fingers_Score(Fi));
    Input_dynamics=Projected_PCs;
    Output_Dynamics=Kin_Vel{1,Fi}(1:length(Projected_PCs));
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end


%% generating one PCA for all duration of finger movement

if 1
Num_PCs=10;

Fi_Sizes=[F1_n F2_n F3_n F4_n F5_n];
% PCA for features
for Fi=1:5
    
    if SpatialFilteringTwo
        
        Data_Finger=FingersFeatures_Spatial(Fi).Finger(:,Brain_Area);
    else
        Data_Finger=FingersFeatures(Fi).Finger(:,Brain_Area);
    end
    
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data_Finger);
    %Fingers_ScoreLG(Fi).Finger=F_score;
      
end

%save('ScoresLG.mat','Fingers_ScoreLG')

% Position prediction & Vel prediction 

Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

 %prediction the position
for Fi=1:5
  
    Projected_PCs=F_score(:,1:Num_PCs);
    Input_dynamics=Projected_PCs;
    Output_Dynamics=Kin_Pos{1,Fi}(1:length(Projected_PCs));
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end

% prediction the velocities
for Fi=1:5
  
    Projected_PCs=F_score(:,1:Num_PCs);
    Input_dynamics=Projected_PCs;
    Output_Dynamics=Kin_Vel{1,Fi}(1:length(Projected_PCs));
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end

end

%% Without any PCA; using filtered data or hilbert features; Position prediction & Vel prediction 

if 0
    
Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

%prediction the position
for Fi=1:5
    
    Input_dynamics=FingersFeatures(Fi).Finger;
    Output_Dynamics=Kin_Pos{1,Fi}(1:end);
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end

% prediction the velocities
for Fi=1:5
  
    Input_dynamics=FingersFeatures(Fi).Finger;
    Output_Dynamics=Kin_Vel{1,Fi}(1:end);
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end

end 

%% combining PCs from different bands to estimate the kinematics

if 0
% Position prediction & Vel prediction 
Common_Dim=40;

Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

 %prediction the position
for Fi=1:5
  
    Projected_PCs=[Fingers_ScoreHG(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreDelta(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreTheta(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreMu(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreBeta(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreLG(Fi).Finger(:,1:Common_Dim)]; % 
                
    Input_dynamics=Projected_PCs;
    Output_Dynamics=Kin_Pos{1,Fi}(1:length(Projected_PCs));
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
      
end

% prediction the velocities
for Fi=1:5
    
    Projected_PCs=[Fingers_ScoreHG(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreDelta(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreTheta(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreMu(Fi).Finger(:,1:Common_Dim),...
                    Fingers_ScoreBeta(Fi).Finger(:,1:Common_Dim)...
                    Fingers_ScoreLG(Fi).Finger(:,1:Common_Dim)]; % 
    
    Input_dynamics=Projected_PCs;
    Output_Dynamics=Kin_Vel{1,Fi}(1:length(Projected_PCs));
    Figures=1;
    Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
    
end

end 

%% performing analysis for phase plane
if 0
    
load('ScoresHG.mat')    

% Preeya method: 

Fi=1;
% CurrentSamples=A*previousSamples
Common_Dim=2;
Samples=Fingers_ScoreDelta(Fi).Finger(:,[1 2]);
Samples=Samples';
Current_Sample=Samples(:,2:end);
Previous_Sample=Samples(:,1:(end-1));
A=Current_Sample*Previous_Sample'*pinv(Previous_Sample*Previous_Sample'); 

% State space model estimation:
% Xdot=A*X
Fi=1;
Samples=Fingers_ScoreDelta(Fi).Finger(:,[1 2]);
% magnify the values
X=1*Samples;
Xdot=diff(X);

X=X(1:(end-1),:);
X=X';
Xdot=Xdot';
A=Xdot*X'*pinv(X*X'); 

% State space model with real values:
figure;
set(gcf, 'Position', [100, 100, 900, 600]);
suptitle(['Feature: PC1 PC2:Trajectory of first-2PC']);
X_PCs=X';
Traj_PC1=X_PCs(1:1000,1);
Traj_PC2=X_PCs(1:1000,2);
u=[0;diff(Traj_PC1)];
v=[0;diff(Traj_PC2)];
quiver(Traj_PC1(1:10:end),Traj_PC2(1:10:end),u(1:10:end),v(1:10:end),'MaxHeadSize',0.05,...
    'LineWidth',1.5,'AutoScale','on','color','b')
xlabel('Real Neural State 1 for PC1')
ylabel('Real Neural State 2 for PC2')




end









    