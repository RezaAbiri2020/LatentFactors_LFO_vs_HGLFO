
%  performing Dimensionality Reduction using ERP and PCA 

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

%Analysis of ECoG dynamics based on which kinematic parameter
ECoGDynamics=1; % based on MaxVel
%ECoGDynamics=2; % based on MaxPos

% which assembly type 
%AssembledData=1; % Using cancatanated raw bands
AssembledData=2;  % Using ERP based for that band 

% The time window for maxvel and maxpos?
window=250;
%window=750;

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

%% F1 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F1=0:1/(Fs):time_AllTrials{1,1}(end);
Posq_F1=interp1(time_AllTrials{1,1},fingers{1,1}(1).spos,Xq_F1);

half_window=200; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F1_n,F1_m]=size(Xq_F1');

% finding max and min time stamps
j=1;
jj=1;
for i=15*half_window:F1_n-(half_window+1)
    
    logic1=(Posq_F1(i-half_window:i+half_window,1)< Posq_F1(i,1));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(1).F2E_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F1(i-half_window:i+half_window,1)> Posq_F1(i,1));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(1).E2F_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F1(:,1))
hold on; plot(FingersKinInfo.Finger(1).F2E_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(1).E2F_MaxPos,-5,'or')

%% F1 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F1=interp1(time_AllTrials{1,1},fingers{1,1}(1).svel,Xq_F1);
AbsVel_F1=sqrt(Velq_F1(:,1).^2+Velq_F1(:,2).^2+Velq_F1(:,3).^2);

Index_F1=sort([1; FingersKinInfo.Finger(1).F2E_MaxPos;...
    FingersKinInfo.Finger(1).E2F_MaxPos; F1_n]);


j=1;
for i=1:2:(length(Index_F1)-1)
    Section=AbsVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(1).F2E_MaxVel(1,j)=(n_S+Index_F1(i)-1);  
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F1)-1)
    Section=AbsVel_F1(Index_F1(i):(Index_F1(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(1).E2F_MaxVel(1,j)=(n_S+Index_F1(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F1)
hold on; plot(FingersKinInfo.Finger(1).F2E_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(1).E2F_MaxVel,150,'ob')

figure;
plot(Posq_F1(:,1))
hold on; plot(FingersKinInfo.Finger(1).F2E_MaxPos,-30,'og')
hold on; plot(FingersKinInfo.Finger(1).E2F_MaxPos,-5,'or')
hold on; vline(FingersKinInfo.Finger(1).F2E_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(1).E2F_MaxVel,'r--')

FingersKinData.Finger(1).PurePos=Posq_F1(:,1);
FingersKinData.Finger(1).PureVel=AbsVel_F1;


%% F2 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F2=0:1/(Fs):time_AllTrials{1,2}(end);
Posq_F2=interp1(time_AllTrials{1,2},fingers{1,2}(2).spos,Xq_F2);

half_window=200; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F2_n,F2_m]=size(Xq_F2');



% finding max and min time stamps
j=1;
jj=1;

for i=10*half_window:F2_n-(half_window+1)
    
    logic1=(Posq_F2(i-half_window:i+half_window,3)< Posq_F2(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(2).F2E_MaxPos(j,1)=i;
        j=j+1;   
    end
    
    logic2=(Posq_F2(i-half_window:i+half_window,3)> Posq_F2(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(2).E2F_MaxPos(jj,1)=i;
        jj=jj+1;  
    end
    
    
end

% check on the plot
figure;
plot(Posq_F2(:,3))
hold on; plot(FingersKinInfo.Finger(2).F2E_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(2).E2F_MaxPos,-50,'or')

%% F2 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F2=interp1(time_AllTrials{1,2},fingers{1,2}(2).svel,Xq_F2);
AbsVel_F2=sqrt(Velq_F2(:,1).^2+Velq_F2(:,2).^2+Velq_F2(:,3).^2);

Index_F2=sort([1; FingersKinInfo.Finger(2).F2E_MaxPos;...
    FingersKinInfo.Finger(2).E2F_MaxPos; F2_n]);

j=1;
for i=1:2:(length(Index_F2)-1)
    Section=AbsVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(2).F2E_MaxVel(1,j)=(n_S+Index_F2(i)-1);
    j=j+1;
    
end

j=1;
for i=2:2:(length(Index_F2)-1)
    Section=AbsVel_F2(Index_F2(i):(Index_F2(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(2).E2F_MaxVel(1,j)=(n_S+Index_F2(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F2)
hold on; plot(FingersKinInfo.Finger(2).F2E_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(2).E2F_MaxVel,150,'ob')

figure;
plot(Posq_F2(:,3))
hold on; plot(FingersKinInfo.Finger(2).F2E_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(2).E2F_MaxPos,-50,'or')
hold on; vline(FingersKinInfo.Finger(2).F2E_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(2).E2F_MaxVel,'r--')

FingersKinData.Finger(2).PurePos=Posq_F2(:,3);
FingersKinData.Finger(2).PureVel=AbsVel_F2;

%% F3 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F3=0:1/(Fs):time_AllTrials{1,3}(end);
Posq_F3=interp1(time_AllTrials{1,3},fingers{1,3}(3).spos,Xq_F3);

half_window=250; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F3_n,F3_m]=size(Xq_F3');

% finding max and min time stamps
j=1;
jj=1;
for i=10*half_window:F3_n-(half_window+1)
    
    logic1=(Posq_F3(i-half_window:i+half_window,3)< Posq_F3(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(3).F2E_MaxPos(j,1)=i;
        j=j+1;  
    end
    
    logic2=(Posq_F3(i-half_window:i+half_window,3)> Posq_F3(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(3).E2F_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

FingersKinInfo.Finger(3).F2E_MaxPos=FingersKinInfo.Finger(3).F2E_MaxPos([1:3,5:end]);
% check on the plot
figure;
plot(Posq_F3(:,3))
hold on; plot(FingersKinInfo.Finger(3).F2E_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(3).E2F_MaxPos,-50,'or')

%% F3 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F3=interp1(time_AllTrials{1,3},fingers{1,3}(3).svel,Xq_F3);
AbsVel_F3=sqrt(Velq_F3(:,1).^2+Velq_F3(:,2).^2+Velq_F3(:,3).^2);

Index_F3=sort([1; FingersKinInfo.Finger(3).F2E_MaxPos;...
    FingersKinInfo.Finger(3).E2F_MaxPos; F3_n]);

j=1;
for i=1:2:(length(Index_F3)-1)
    Section=AbsVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(3).F2E_MaxVel(1,j)=(n_S+Index_F3(i)-1);
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F3)-1)
    Section=AbsVel_F3(Index_F3(i):(Index_F3(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(3).E2F_MaxVel(1,j)=(n_S+Index_F3(i)-1);
    j=j+1;
end

% check on the plot
figure;
plot(AbsVel_F3)
hold on; plot(FingersKinInfo.Finger(3).F2E_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(3).E2F_MaxVel,150,'ob')

figure;
plot(Posq_F3(:,3))
hold on; plot(FingersKinInfo.Finger(3).F2E_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(3).E2F_MaxPos,-50,'or')
hold on; vline(FingersKinInfo.Finger(3).F2E_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(3).E2F_MaxVel,'r--')

FingersKinData.Finger(3).PurePos=Posq_F3(:,3);
FingersKinData.Finger(3).PureVel=AbsVel_F3;

%% F4 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F4=0:1/(Fs):time_AllTrials{1,4}(end);
Posq_F4=interp1(time_AllTrials{1,4},fingers{1,4}(4).spos,Xq_F4);

half_window=200; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F4_n,F4_m]=size(Xq_F4');

% finding max and min time stamps
j=1;
jj=1;
for i=10*half_window:F4_n-(half_window+1)
    
    logic1=(Posq_F4(i-half_window:i+half_window,3)< Posq_F4(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(4).F2E_MaxPos(j,1)=i;
        j=j+1;
        
    end
    
    logic2=(Posq_F4(i-half_window:i+half_window,3)> Posq_F4(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(4).E2F_MaxPos(jj,1)=i;
        jj=jj+1;
    end
    
    
end

% check on the plot
figure;
plot(Posq_F4(:,3))
hold on; plot(FingersKinInfo.Finger(4).F2E_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(4).E2F_MaxPos,-50,'or')

%% F4 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F4=interp1(time_AllTrials{1,4},fingers{1,4}(4).svel,Xq_F4);
AbsVel_F4=sqrt(Velq_F4(:,1).^2+Velq_F4(:,2).^2+Velq_F4(:,3).^2);

Index_F4=sort([1; FingersKinInfo.Finger(4).F2E_MaxPos;...
    FingersKinInfo.Finger(4).E2F_MaxPos; F4_n]);

j=1;
for i=1:2:(length(Index_F4)-1)
    Section=AbsVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(4).F2E_MaxVel(1,j)=(n_S+Index_F4(i)-1);
    j=j+1;
    
end

j=1;
for i=2:2:(length(Index_F4)-1)
    Section=AbsVel_F4(Index_F4(i):(Index_F4(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(4).E2F_MaxVel(1,j)=(n_S+Index_F4(i)-1);
    j=j+1;
end


% check on the plot
figure;
plot(AbsVel_F4)
hold on; plot(FingersKinInfo.Finger(4).F2E_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(4).E2F_MaxVel,150,'ob')

figure;
plot(Posq_F4(:,3))
hold on; plot(FingersKinInfo.Finger(4).F2E_MaxPos,-80,'og')
hold on; plot(FingersKinInfo.Finger(4).E2F_MaxPos,-50,'or')
hold on; vline(FingersKinInfo.Finger(4).F2E_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(4).E2F_MaxVel,'r--')

FingersKinData.Finger(4).PurePos=Posq_F4(:,3);
FingersKinData.Finger(4).PureVel=AbsVel_F4;

%% F5 movement analysis: Finding the peak points for transition from flextion to extention or extension to flextion

Xq_F5=0:1/(Fs):time_AllTrials{1,5}(end);
Posq_F5=interp1(time_AllTrials{1,5},fingers{1,5}(5).spos,Xq_F5);

half_window=350; %points; about 1 sec: for finding local max and local min in position data

% for finger 1 
[F5_n,F5_m]=size(Xq_F5');

% finding max and min time stamps
j=1;
jj=1;
for i=5*half_window:F5_n-(half_window+1)
    
    logic1=(Posq_F5(i-half_window:i+half_window,3)< Posq_F5(i,3));
    
    if sum(logic1)==2*half_window
        FingersKinInfo.Finger(5).F2E_MaxPos(j,1)=i;
        j=j+1;
    end
    
    logic2=(Posq_F5(i-half_window:i+half_window,3)> Posq_F5(i,3));
    
    if sum(logic2)==2*half_window
        FingersKinInfo.Finger(5).E2F_MaxPos(jj,1)=i;
        jj=jj+1;
        
    end
     
end

FingersKinInfo.Finger(5).E2F_MaxPos=FingersKinInfo.Finger(5).E2F_MaxPos([1:16,18:end]);

% check on the plot
figure;
plot(Posq_F5(:,3))
hold on; plot(FingersKinInfo.Finger(5).F2E_MaxPos,-50,'og')
hold on; plot(FingersKinInfo.Finger(5).E2F_MaxPos,-20,'or')

%% F5 movement analysis: Finding the max velocity during transition from flextion to extention or extension to flextion

Velq_F5=interp1(time_AllTrials{1,5},fingers{1,5}(5).svel,Xq_F5);
AbsVel_F5=sqrt(Velq_F5(:,1).^2+Velq_F5(:,2).^2+Velq_F5(:,3).^2);

Index_F5=sort([1; FingersKinInfo.Finger(5).F2E_MaxPos;...
    FingersKinInfo.Finger(5).E2F_MaxPos; F5_n-(half_window+1)]);

j=1;
for i=1:2:(length(Index_F5)-1)
    Section=AbsVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(5).F2E_MaxVel(1,j)=(n_S+Index_F5(i)-1);
    j=j+1;
end

j=1;
for i=2:2:(length(Index_F5)-1)
    Section=AbsVel_F5(Index_F5(i):(Index_F5(i+1)));
    [m_S n_S]=max(Section);
    FingersKinInfo.Finger(5).E2F_MaxVel(1,j)=(n_S+Index_F5(i)-1);
    j=j+1;
end

% check on the plot
figure;
plot(AbsVel_F5)
hold on; plot(FingersKinInfo.Finger(5).F2E_MaxVel,100,'ok')
hold on; plot(FingersKinInfo.Finger(5).E2F_MaxVel,150,'ob')

figure;
plot(Posq_F5(:,3))
hold on; plot(FingersKinInfo.Finger(5).F2E_MaxPos,-50,'og')
hold on; plot(FingersKinInfo.Finger(5).E2F_MaxPos,-20,'or')
hold on; vline(FingersKinInfo.Finger(5).F2E_MaxVel,'g--')
hold on; vline(FingersKinInfo.Finger(5).E2F_MaxVel,'r--')

FingersKinData.Finger(5).PurePos=Posq_F5(:,3);
FingersKinData.Finger(5).PureVel=AbsVel_F5;

%% Referencing & Filtering options for generating raw bands for all Fingers

Modification=[length(Xq_F1),length(Xq_F2),length(Xq_F3),length(Xq_F4),length(Xq_F5)];

for Fi=1:5
    % modification to sync
    ECoGdata_Modified=lfp(1:Nch_Record,trial_start(Fi):(trial_start(Fi)+Modification(Fi)-1));
    % z score
    ECoGdata_1=zscore(ECoGdata_Modified')';
    
    % option median or mean:
    switch Reference
        
        case 1
            F_Median =median(ECoGdata_1(Brain_Area,:),1);
            ECoGdata_2= ECoGdata_1-repmat(F_Median,Nch_Record,1);
            
        case 2
            F_Mean =mean(ECoGdata_1(Brain_Area,:),1);
            ECoGdata_2= ECoGdata_1-repmat(F_Mean,Nch_Record,1);
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
            F_EnvelopeAll=[];
            
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
                F_EnvelopeAll(:,:,band)=F_Envelope;
                
            end
            F_EnvMean=mean(F_EnvelopeAll,3);
            % Hilbert & Filter
            [bb,aa]=butter(3,[.5, 4]/(Fs/2));
            F_Filtered=filtfilt(bb,aa,F_EnvMean);      
            F_Hilbert=abs(hilbert(F_Filtered));
    end
       
    % pulling out Feature:Filtered data or Feature:Hilberet
    switch Feature
        
        case 1
            FingersFeatures(Fi).Finger=F_Filtered;
            
        case 2
            FingersFeatures(Fi).Finger=F_Hilbert;       
    end
     
    % choose the option for analysis E2F or F2E or both
    % Pulling out features: position-based pulling out or Max-Vel-based pulling out
    Features=FingersFeatures(Fi).Finger;
    
    switch ECoGDynamics
        
        case 1 % maxvel-based analysis
            % F2E
            for i=1:length(FingersKinInfo.Finger(Fi).F2E_MaxVel)
                FingersECoG(Fi).FingerF2E{i,1}=Features(FingersKinInfo.Finger(Fi).F2E_MaxVel(i)-window:...
                    FingersKinInfo.Finger(Fi).F2E_MaxVel(i)+window,:);
            end
            %E2F
            for i=1:length(FingersKinInfo.Finger(Fi).E2F_MaxVel)
                FingersECoG(Fi).FingerE2F{i,1}=Features(FingersKinInfo.Finger(Fi).E2F_MaxVel(i)-window:...
                    FingersKinInfo.Finger(Fi).E2F_MaxVel(i)+window,:);
            end
            
        case 2 % maxpos-based analysis
            %F2E
            for i=1:length(FingersKinInfo.Finger(Fi).F2E_MaxPos)
                FingersECoG(Fi).FingerF2E{i,1}=Features(FingersKinInfo.Finger(Fi).F2E_MaxPos(i)-window:...
                    FingersKinInfo.Finger(Fi).F2E_MaxPos(i)+window,:);
            end
            %F2E
            for i=1:length(FingersKinInfo.Finger(Fi).E2F_MaxPos)
                FingersECoG(Fi).FingerE2F{i,1}=Features(FingersKinInfo.Finger(Fi).E2F_MaxPos(i)-window:...
                    FingersKinInfo.Finger(Fi).E2F_MaxPos(i)+window,:);
            end
                   
    end  
    
    % cancatanteting or ERP 
    
    switch AssembledData
        
        case 1 % cancatanteting
            %F2E
            ECoG_Final(Fi).FingerF2E=cell2mat(FingersECoG(Fi).FingerF2E);
            %E2F
            ECoG_Final(Fi).FingerE2F=cell2mat(FingersECoG(Fi).FingerE2F);
            
        case 2  % ERP
            %F2E
            for i=1:length(FingersECoG(Fi).FingerF2E)
                ECoG_portionF2E(:,:,i)=FingersECoG(Fi).FingerF2E{i,1};
            end
            %E2F
            for i=1:length(FingersECoG(Fi).FingerE2F)
                ECoG_portionE2F(:,:,i)=FingersECoG(Fi).FingerE2F{i,1};
            end
            
            ECoG_Final(Fi).FingerF2E=mean(ECoG_portionF2E,3);
            ECoG_Final(Fi).FingerE2F=mean(ECoG_portionE2F,3);        
    end
     
end


%% 1- Return to Modify Kinematics data for Pos-based analysis; 2- Add single/concatenated Kinematics data for max-vel-based analysis
% Fingers: Finding and including kinematics data
% generating the same sampling points for kinematics

for Fi=1:5     
            
    % maxvel-based windows
    %F2E
    for i=1:length(FingersKinInfo.Finger(Fi).F2E_MaxVel)
        FingersKinData.FingerF2E(Fi).MaxVelWins{i,1}=FingersKinData.Finger(Fi).PureVel(...
            FingersKinInfo.Finger(Fi).F2E_MaxVel(i)-window:FingersKinInfo.Finger(Fi).F2E_MaxVel(i)+window);
    end
    %E2F
    for i=1:length(FingersKinInfo.Finger(Fi).E2F_MaxVel)
        FingersKinData.FingerE2F(Fi).MaxVelWins{i,1}=FingersKinData.Finger(Fi).PureVel(...
            FingersKinInfo.Finger(Fi).E2F_MaxVel(i)-window:FingersKinInfo.Finger(Fi).E2F_MaxVel(i)+window);
    end
    
    % maxpos-based windows
    %F2E
    for i=1:length(FingersKinInfo.Finger(Fi).F2E_MaxPos)
        FingersKinData.FingerF2E(Fi).MaxPosWins{i,1}=FingersKinData.Finger(Fi).PurePos(...
            FingersKinInfo.Finger(Fi).F2E_MaxPos(i)-window:FingersKinInfo.Finger(Fi).F2E_MaxPos(i)+window);
    end
    %E2F
    for i=1:length(FingersKinInfo.Finger(Fi).E2F_MaxPos)
        FingersKinData.FingerE2F(Fi).MaxPosWins{i,1}=FingersKinData.Finger(Fi).PurePos(...
            FingersKinInfo.Finger(Fi).E2F_MaxPos(i)-window:FingersKinInfo.Finger(Fi).E2F_MaxPos(i)+window);
    end
    
    % concatanating and making equal size for later analysis
    FingersKinData.FingerF2E(Fi).CatMaxVels=cell2mat(FingersKinData.FingerF2E(Fi).MaxVelWins);
    FingersKinData.FingerE2F(Fi).CatMaxVels=cell2mat(FingersKinData.FingerE2F(Fi).MaxVelWins);
    FingersKinData.FingerF2E(Fi).CatMaxPos=cell2mat(FingersKinData.FingerF2E(Fi).MaxPosWins);
    FingersKinData.FingerE2F(Fi).CatMaxPos=cell2mat(FingersKinData.FingerE2F(Fi).MaxPosWins);
    
    
    Cat_Size=min([size(FingersKinData.FingerF2E(Fi).CatMaxVels,1),size(FingersKinData.FingerF2E(Fi).CatMaxPos,1),...
        size(FingersKinData.FingerE2F(Fi).CatMaxVels,1),size(FingersKinData.FingerE2F(Fi).CatMaxPos,1)]);
    
    FingersKinData.FingerF2E(Fi).CatMaxVels=FingersKinData.FingerF2E(Fi).CatMaxVels(1:Cat_Size,1);
    FingersKinData.FingerE2F(Fi).CatMaxVels=FingersKinData.FingerE2F(Fi).CatMaxVels(1:Cat_Size,1);
    FingersKinData.FingerF2E(Fi).CatMaxPos=FingersKinData.FingerF2E(Fi).CatMaxPos(1:Cat_Size,1); 
    FingersKinData.FingerE2F(Fi).CatMaxPos=FingersKinData.FingerE2F(Fi).CatMaxPos(1:Cat_Size,1); 
      
end


%% PCA analysis of cancatanated or ERP results


% PCA for ECoG F2E
for Fi=1:5
    Data=ECoG_Final(Fi).FingerF2E(:,Brain_Area);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.FingerF2E(Fi).Coeff=F_coeff;
    PCA.FingerF2E(Fi).Score=F_score;
    PCA.FingerF2E(Fi).Explained=F_explained;
    PCA.FingerF2E(Fi).Variability=F_variability;
    
    NumDim=50;
    figure;
    subplot(1,2,1)
    bar(F_explained)
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('PC #')
    
    subplot(1,2,2)
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
     
    % Putting zero for bad channels for performing brain plots
    %Modify the number of rows:
    j=0;
    for i=1:Nch_Record
        if Brain_Area(i)==1
            PCA.FingerF2E(Fi).CoeffModified(i,:)=(F_coeff(i-j,:));  
        else
            PCA.FingerF2E(Fi).CoeffModified(i,:)=zeros(1,size(F_coeff,1));
            j=j+1;
        end    
    end
    %Modify the number of columns:
    F_coeff_Modified=[];
    F_coeff_Modified=PCA.FingerF2E(Fi).CoeffModified;
    j=0;
    for i=1:Nch_Record
        if Brain_Area(i)==1
            PCA.FingerF2E(Fi).CoeffModified(:,i)=(F_coeff_Modified(:,i-j));
        else
            PCA.FingerF2E(Fi).CoeffModified(:,i)=zeros(Nch_Record,1);
            j=j+1;
        end
    end
    
       
end

% PCA for ECoG E2F
for Fi=1:5
    Data=ECoG_Final(Fi).FingerE2F(:,Brain_Area);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.FingerE2F(Fi).Coeff=F_coeff;
    PCA.FingerE2F(Fi).Score=F_score;
    PCA.FingerE2F(Fi).Explained=F_explained;
    PCA.FingerE2F(Fi).Variability=F_variability;
    
    NumDim=50;
    figure;
    subplot(1,2,1)
    bar(F_explained)
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('PC #')
    
    subplot(1,2,2)
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
     
    % Putting zero for bad channels for performing brain plots
    %Modify the number of rows:
    j=0;
    for i=1:Nch_Record
        if Brain_Area(i)==1
            PCA.FingerE2F(Fi).CoeffModified(i,:)=(F_coeff(i-j,:));  
        else
            PCA.FingerE2F(Fi).CoeffModified(i,:)=zeros(1,size(F_coeff,1));
            j=j+1;
        end    
    end
    %Modify the number of columns:
    F_coeff_Modified=[];
    F_coeff_Modified=PCA.FingerE2F(Fi).CoeffModified;
    j=0;
    for i=1:Nch_Record
        if Brain_Area(i)==1
            PCA.FingerE2F(Fi).CoeffModified(:,i)=(F_coeff_Modified(:,i-j));
        else
            PCA.FingerE2F(Fi).CoeffModified(:,i)=zeros(Nch_Record,1);
            j=j+1;
        end
    end
    
       
end

%% Position prediction & Vel prediction using cancatanteting or ERP results for windows
titles={'PCA; Finger 1: Thumb', 'PCA; Finger 2: Index','PCA; Finger 3: Middle','PCA; Finger 4:Ring', 'PCA; Finger 5: Pinkie'};

% projection for observation
Num_PCs=10;

% predictions based on max-vel or max-pos analysis of ECoG
% F2E ECoG dynamics
for prediction=1:4
    
    for Fi=1:5
        % Finding the first nonzeros Coeffs
        All_PCs=PCA.FingerF2E(Fi).CoeffModified;
        Targeted_PCs=All_PCs(:,Brain_Area);
        Employed_PCs=Targeted_PCs(:,1:Num_PCs);
        
        if prediction==1
            %prediction the maxVelsWins F2E using F2E dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerF2E)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerF2E(Fi).CatMaxVels;
            
        elseif prediction==2
            %prediction the maxposWins F2E using F2E dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerF2E)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerF2E(Fi).CatMaxPos;
            
        elseif prediction==3
            %prediction the maxVelWins E2F using F2E dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerF2E)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerE2F(Fi).CatMaxVels;
            
        elseif prediction==4
            %prediction the maxposWins E2F using F2E dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerF2E)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerE2F(Fi).CatMaxPos;
            
        end
        
        % making consistency in size for regression
        Reg_Size=min(size(Input_dynamics,1),size(Output_Dynamics,1));
        Input_dynamics=Input_dynamics(1:Reg_Size,:);
        Output_Dynamics=Output_Dynamics(1:Reg_Size,:);
        
        Figures=1;
        Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
    end
    
end

% E2F ECoG dynamics
for prediction=1:4
    
    for Fi=1:5
        % Finding the first nonzeros Coeffs
        All_PCs=PCA.FingerE2F(Fi).CoeffModified;
        Targeted_PCs=All_PCs(:,Brain_Area);
        Employed_PCs=Targeted_PCs(:,1:Num_PCs);
        
        if prediction==1
            %prediction the maxVelsWins E2F using E2F dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerE2F)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerE2F(Fi).CatMaxVels;
            
        elseif prediction==2
            %prediction the maxposWins E2F using E2F dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerE2F)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerE2F(Fi).CatMaxPos;
            
        elseif prediction==3
            %prediction the maxVelWins F2E using E2F dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerE2F)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerF2E(Fi).CatMaxVels;
            
        elseif prediction==4
            %prediction the maxposWins F2E using E2F dynamics of ECoG
            Input_dynamics=cell2mat(FingersECoG(Fi).FingerE2F)*Employed_PCs;
            Output_Dynamics=FingersKinData.FingerF2E(Fi).CatMaxPos;
            
        end
        
        % making consistency in size for regression
        Reg_Size=min(size(Input_dynamics,1),size(Output_Dynamics,1));
        Input_dynamics=Input_dynamics(1:Reg_Size,:);
        Output_Dynamics=Output_Dynamics(1:Reg_Size,:);
        
        Figures=1;
        Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
    end
    
end
%% Position prediction & Vel prediction using cancatanteting or ERP results for all movements

titles={'PCA; Finger 1: Thumb', 'PCA; Finger 2: Index','PCA; Finger 3: Middle','PCA; Finger 4:Ring', 'PCA; Finger 5: Pinkie'};

% projection for observation
Num_PCs=10;

% predictions based on max-vel or max-pos analysis of ECoG
% F2E ECoG dynamics
for prediction=2%1:2
    
    for Fi=1:5
        % Finding the first nonzeros Coeffs
        All_PCs=PCA.FingerF2E(Fi).CoeffModified;
        Targeted_PCs=All_PCs(:,Brain_Area);
        Employed_PCs=Targeted_PCs(:,1:Num_PCs);
        
        if prediction==1
            %prediction the whole Vels using F2E dynamics of ECoG
            Input_dynamics=FingersFeatures(Fi).Finger*Employed_PCs;
            Output_Dynamics=FingersKinData.Finger(Fi).PureVel;
            
        elseif prediction==2
            %prediction the whole pos using F2E dynamics of ECoG
            Input_dynamics=FingersFeatures(Fi).Finger*Employed_PCs;
            Output_Dynamics=FingersKinData.Finger(Fi).PurePos;
            
        end
        
        Figures=1;
        Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
    end
    
end

% E2F ECoG dynamics
for prediction=1%:2
    
    for Fi=1:5
        % Finding the first nonzeros Coeffs
        All_PCs=PCA.FingerE2F(Fi).CoeffModified;
        Targeted_PCs=All_PCs(:,Brain_Area);
        Employed_PCs=Targeted_PCs(:,1:Num_PCs);
        
        if prediction==1
            %prediction the whole Vels using E2F dynamics of ECoG
            Input_dynamics=FingersFeatures(Fi).Finger*Employed_PCs;
            Output_Dynamics=FingersKinData.Finger(Fi).PureVel;
            
        elseif prediction==2
            %prediction the whole pos using E2F dynamics of ECoG
            Input_dynamics=FingersFeatures(Fi).Finger*Employed_PCs;
            Output_Dynamics=FingersKinData.Finger(Fi).PurePos;
 
        end

        Figures=1;
        Rsq_Pre_output=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
    end
    
end


%% phase plan analysis

% State space model estimation:
% Xdot=A*X
Fi=5;
Samples=PCA.FingerE2F(Fi).Score(:,[1 2]);
% magnify the values
X=1*Samples;
Xdot=diff(X);
X=X(1:(end-1),:);
X=X';
Xdot=Xdot';
A=Xdot*X'*pinv(X*X'); 
% solving for phase plane
save('AMatrix.mat','A')
tspan=[0,4000];
icond={[1, 0]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',0.01);


% the real trajectories 
figure;
set(gcf, 'Position', [100, 100, 900, 600]);
suptitle(['Feature: PC1 PC2:Trajectory of first-2PC']);

Traj_PC1=Input_dynamics(503:1003,1);
Traj_PC2=Input_dynamics(503:1003,2);
u=[0;diff(Traj_PC1)];
v=[0;diff(Traj_PC2)];
quiver(Traj_PC1(1:10:end),Traj_PC2(1:10:end),u(1:10:end),v(1:10:end),'MaxHeadSize',0.05,...
    'LineWidth',1.5,'AutoScale','on','color','b')
xlabel('Neural State 1 for PC1')
ylabel('Neural State 2 for PC2')
zlabel('Neural State 3 for PC3')

%for high gamma
% eig(A)
% ans =
% 
%    0.0039 + 0.0063i
%    0.0039 - 0.0063i

% for delta
% eig(A)
% 
% ans =
% 
%    0.0001 + 0.0095i
%    0.0001 - 0.0095i

        
    