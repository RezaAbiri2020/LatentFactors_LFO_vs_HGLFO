% Analyzing the MLR for subbands of high gamma; is there any patterns in
% lags predicttion

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
load('E:\ECoGLeapMotion\DataPatientTwo\ECoGData\ECoG_data.mat');

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

%% Finding the precentral (for Primary Motor=PM) and postcentral (somatosensory=SM) electrodes/channels

 load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')
 
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
    filename = 'E:\ECoGLeapMotion\DataPatientTwo/LeapMotionData/20180301/10h02m38s/EC171_20180301_10h02m38s_trial00';
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
    fprintf(['Finger:',num2str(Fi),'\n'])
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
            Delta_Fingers(Fi).Filtered=F_Filtered;
            F_Hilbert=abs(hilbert(F_Filtered));
            Delta_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 2 % theta
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            Theta_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 3 % alpha
            [b,a]=butter(3,[8, 13]/(Fs/2));
            F_Filtered=filtfilt(b,a,ECoGdata_2');
            F_Hilbert=abs(hilbert(F_Filtered));
            Alpha_Fingers(Fi).Hilbert=F_Hilbert;
            
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
                [b,a]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertAll(:,:,band)=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                
            end
            
            F_Filtered=mean(F_FilteredAll,3);
            F_Hilbert=mean(F_HilbertAll,3);
            Beta_Fingers(Fi).Hilbert=F_Hilbert;
            
        case 5 % low gamma
            LGbands={
                [30,36]
                [36,42]
                [42,50]};
            
            % for whole band
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[30 50]/(Fs/2));
            F_FilteredAll=filtfilt(b,a,ECoGdata_2');
           
            % Envelope & Hilbert
            [F_Envelope,Lower]=envelope(F_FilteredAll);
            LG_FingersWholeBand(Fi).PureEnv=F_Envelope;
            % Hilbert & Filter in delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            LG_FingersWholeBand(Fi).EnvDeltaWOutMean=F_HilbertAll;
            
            % Hilbert & Filter in theta
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            LG_FingersWholeBand(Fi).EnvThetaWOutMean=F_HilbertAll;
            
            % for different sub bands
            AllSubBands=[];
            for band=1:3
                [b,a]=butter(3,LGbands{band}/(Fs/2));
                F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                               
                % Envelope & Hilbert
                [F_Envelope,Lower]=envelope(F_FilteredAll);
                LG_Fingers(Fi).SubBands(band).PureEnv=F_Envelope;
                AllSubBands(:,:,band)=F_Envelope;
                
                % Hilbert & Filter for delta
                [b,a]=butter(3,[.5, 4]/(Fs/2));
                F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                LG_Fingers(Fi).SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
               
                % Hilbert & Filter for theta
                [b,a]=butter(3,[4, 8]/(Fs/2));
                F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                LG_Fingers(Fi).SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                
            end
            
            % for the averege of bands
            LG_FingersAvgSubBands(Fi).PureEnv=mean(AllSubBands,3);
            % Hilbert & Filter for delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_HilbertDelta1=filtfilt(b,a,mean(AllSubBands,3));
            LG_FingersAvgSubBands(Fi).EnvDeltaWOutMean=F_HilbertDelta1;
            
            % Hilbert & Filter for theta
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_HilbertTheta1=filtfilt(b,a,mean(AllSubBands,3));
            LG_FingersAvgSubBands(Fi).EnvThetaWOutMean=F_HilbertTheta1;
            
        case 6
            
            HGbands={
                [70,77]
                [77,85]
                [85,93]
                [93,102]
                [102,113]
                [113,124]
                [124,136]
                [136,150]};
            
            % for whole band
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[70 150]/(Fs/2));
            F_FilteredAll=filtfilt(b,a,ECoGdata_2');
            HG_FingersWholeBand(Fi).Filtered=F_FilteredAll;
            % Envelope & Hilbert
            [F_Envelope,Lower]=envelope(F_FilteredAll);
            HG_FingersWholeBand(Fi).PureEnv=F_Envelope;
            % Hilbert & Filter in delta
            [b,a]=butter(3,[.5, 4]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvDeltaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvDeltaWithMean=F_HilbertAll;
            
            % Hilbert & Filter in theta
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[4, 8]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvThetaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvThetaWithMean=F_HilbertAll;
            
            % Hilbert & Filter in alpha
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[8, 13]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvAlphaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvAlphaWithMean=F_HilbertAll;
            
            % Hilbert & Filter in beta
            F_FilteredAll=[];
            F_HilbertAll=[];
            [b,a]=butter(3,[13, 30]/(Fs/2));
            F_HilbertAll=filtfilt(b,a,F_Envelope);
            HG_FingersWholeBand(Fi).EnvBetaWOutMean=F_HilbertAll;
            F_HilbertAll=[];
            F_HilbertAll=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
            HG_FingersWholeBand(Fi).EnvBetaWithMean=F_HilbertAll;
            
            
            % for different sub bands
            if Fi==1
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger1SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger1SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger1SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger1SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger1SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger1SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger1SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==2
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger2SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger2SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger2SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger2SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger2SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger2SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger2SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==3
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger3SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger3SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger3SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger3SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger3SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger3SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger3SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==4
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger4SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger4SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger4SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger4SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger4SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger4SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger4SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            elseif Fi==5
                for band=1:8
                    [b,a]=butter(3,HGbands{band}/(Fs/2));
                    F_FilteredAll=filtfilt(b,a,ECoGdata_2');
                    HG_Finger5SubBands(band).Filtered=F_FilteredAll;
                    
                    % Envelope & Hilbert
                    [F_Envelope,Lower]=envelope(F_FilteredAll);
                    HG_Finger5SubBands(band).PureEnv=F_Envelope;
                    % Hilbert & Filter for delta
                    [b,a]=butter(3,[.5, 4]/(Fs/2));
                    F_HilbertDelta1=filtfilt(b,a,F_Envelope);
                    F_HilbertDelta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvDeltaWOutMean=F_HilbertDelta1;
                    HG_Finger5SubBands(band).EnvDeltaWithMean=F_HilbertDelta2;
                    
                    % Hilbert & Filter for theta
                    [b,a]=butter(3,[4, 8]/(Fs/2));
                    F_HilbertTheta1=filtfilt(b,a,F_Envelope);
                    F_HilbertTheta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvThetaWOutMean=F_HilbertTheta1;
                    HG_Finger5SubBands(band).EnvThetaWithMean=F_HilbertTheta2;
                    
                    % Hilbert & Filter for alpha
                    [b,a]=butter(3,[8, 13]/(Fs/2));
                    F_HilbertAlpha1=filtfilt(b,a,F_Envelope);
                    F_HilbertAlpha2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvAlphaWOutMean=F_HilbertAlpha1;
                    HG_Finger5SubBands(band).EnvAlphaWithMean=F_HilbertAlpha2;
                    
                    % Hilbert & Filter for beta
                    [b,a]=butter(3,[13, 30]/(Fs/2));
                    F_HilbertBeta1=filtfilt(b,a,F_Envelope);
                    F_HilbertBeta2=filtfilt(b,a,F_Envelope)+repmat(mean(F_Envelope),size(F_Envelope,1),1);
                    HG_Finger5SubBands(band).EnvBetaWOutMean=F_HilbertBeta1;
                    HG_Finger5SubBands(band).EnvBetaWithMean=F_HilbertBeta2;
                    
                end
            end      
    end
      
     
end

%% Choose Channels: In case of choosing top predictive channels

load('HG_AllWeights.mat');
load('Delta_AllWeights.mat');

% observaing and comparing the weights of R2 on grid for all channels
Nch_Record=256;

% grid layout
Ch_num_1=1:Nch_Record;
Ch_num_2=reshape(Ch_num_1,[16,16]);
% Ch_num_3 is the grid layout
Ch_num_3=rot90(rot90(Ch_num_2));

[R,C]=size(Ch_num_3);

% for vel for high gamma 
for Fi=5
    Weights=HG_WeightsVel(Fi).Finger;
    Weights=Weights(2:end);
    Layout=zeros(R,C);
    
    for r=1:R
        for c=1:C
            Layout(r,c)=Weights(Ch_num_3(r,c),1);
        end
    end
    figure;
    imagesc(Layout)
    title(['Finger: ',num2str(Fi)])
    colorbar
    
end
[Values, ChNo]=sort(abs(Weights),'descend');

% For vel: by observation in Finger 4, 5; see & check Ch:187 (highly negative) Ch:186 (highly positive)


% for pos for high gamma 
for Fi=1:5
    Weights=HG_WeightsPos(Fi).Finger;
    Weights=Weights(2:end);
    Layout=zeros(R,C);
    
    for r=1:R
        for c=1:C
            Layout(r,c)=Weights(Ch_num_3(r,c),1);
        end
    end
    figure;
    imagesc(Layout)
    title(['Finger: ',num2str(Fi)])
    colorbar
    
end

% For Pos: by observation in Finger 1; see & check Ch:49 (highly negative) Ch:50 (highly positive)

%% MLR using lags in subbands & record the predicted weights for each subband 

Kin_Pos={Posq_F1(:,1) Posq_F2(:,3) Posq_F3(:,3) Posq_F4(:,3) Posq_F5(:,3)};
Kin_Vel={AbsVel_F1 AbsVel_F2 AbsVel_F3 AbsVel_F4 AbsVel_F5};

%ch=logical(All_Index_Hand);
ch=[90:94,106:110];
%ch=[95];
Subbands=[];
for band=[1,3,6]
    
    Hilbert_band=abs(hilbert(HG_Finger5SubBands(band).EnvDeltaWOutMean(:,ch)));
    %Hilbert_band=HG_Finger5SubBands(band).PureEnv(:,ch);
    Hilbert_band_Resample=Hilbert_band;
    % for downsampling
    %Hilbert_band_Resample=resample(Hilbert_band,80,508); % resample at rate of (50/Fs)*Fs
    
    Subbands=[Subbands,Hilbert_band_Resample]; 
end

% generate the lags
LagN=2;
Input_dynamics=Subbands;
for i=1:LagN
    Bands_Lag=[zeros(i,size(Subbands,2)); Subbands(1:(end-i),:)];
    Input_dynamics=[Input_dynamics,Bands_Lag]; 
end 

Fi=5;
Figures=1;
WeightVel=[];
WeightPos=[];
Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
Output_DynamicsVel_Resample=Output_DynamicsVel;
%Output_DynamicsVel_Resample=resample(Output_DynamicsVel,80,508);
Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
Output_DynamicsPos_Resample=Output_DynamicsPos;
%Output_DynamicsPos_Resample=resample(Output_DynamicsPos,80,508);
WeightVel=MultipleRegFuncCalWeight(Fi,Input_dynamics,Output_DynamicsVel_Resample,Figures);
WeightPos=MultipleRegFuncCalWeight(Fi,Input_dynamics,Output_DynamicsPos_Resample,Figures);

for i=1:10%sum(ch)
    figure;
    set(gcf, 'Position', [300, 300, 1200, 600]);
    Weight_Lag=reshape(WeightVel(i+1:10:end),[3,LagN+1]);
    imagesc((Weight_Lag))
    xlabel('Lag')
    ylabel('Subband')
    colorbar
end

for i=1:10%sum(ch)
    figure;
    set(gcf, 'Position', [300, 300, 1200, 600]);
    Weight_Lag=reshape(WeightPos(i+1:10:end),[3,LagN+1]);
    imagesc((Weight_Lag))
    xlabel('Lag')
    ylabel('Subband')
    colorbar
end


%% using xcorr function with lags for one channel: 

% between kinematics and subbands 
Fi=5;
Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
LagN=500;
Colors={'r','b','y','m','c','g','k',[0.8500 0.3250 0.0980]};

for ch=90:94%[90:94,106:110]
    figure;
    set(gcf, 'Position', [300, 300, 700, 600]);
    for band=[1,3,6]
        
        Hilbert_band=abs(hilbert(HG_Finger5SubBands(band).EnvDeltaWOutMean(:,ch)));
        [R,Lags]=xcorr(Output_DynamicsPos,Hilbert_band,LagN,'normalized');
        plot(Lags,R,'color',Colors{band},'linewidth',1.5)
        hold on
    end
    %legend('band1','band2','band3','band4','band5','band6','band7','band8')
    legend('band1','band3','band6')
    xlabel('Lag')
    ylabel('Rcorr')
    %ylim([0.5,1])
    title(['Finger: ',num2str(Fi),'; Ch: ',num2str(ch)])
    set(gca,'fontsize',14)
end

% between LFO and subbands of HG-LFO
 
Fi=5;
LagN=400;
Colors={'r','b','y','m','c','g','k',[0.8500 0.3250 0.0980]};
figure;
set(gcf, 'Position', [100, 100, 900, 900]);
k=0;
for ch=106:110%[90:94,106:110]
    k=k+1;
    subplot(3,2,k)
    for band=[4,8]
        
        Hilbert_band=abs(hilbert(HG_Finger5SubBands(band).EnvDeltaWOutMean(:,ch)));
        [R,Lags]=xcorr(Hilbert_band,LFO_signals(Fi).Hilbert(:,ch),LagN,'normalized');
        plot(Lags,R,'color',Colors{band},'linewidth',1.5)
        [value Index]=max(R);
        hold on
        vline(Index-LagN,'color',Colors{band})
        hold on
    end
    %legend('band1','band2','band3','band4','band5','band6','band7','band8')
    if k==1
        legend('band4','peak4','band8','peak8')
    end
    xlabel('Lag')
    ylabel('Rcorr')
    xticks(-400:200:400)
    %ylim([0.5,1])
    title(['Finger: ',num2str(Fi),'; Ch: ',num2str(ch)])
    set(gca,'fontsize',14)
    
end

HighQualityFigs('Xcorr_LFO_HGLFO_Finger5_1')

%% using xcorr function with lags for whole grid:
Fi=5;
Output_DynamicsVel=Kin_Vel{1,Fi}(1:end);
Output_DynamicsPos=Kin_Pos{1,Fi}(1:end);
LagN=500;

for band=1:8
    figure;
    set(gcf, 'Position', [100, 100, 1000, 800]);
    suptitle(['xcorr (Velocity-Subband',num2str(band),'); lag500; Finger: ',num2str(Fi)]);
    
    for i=1:256
        subplot(16,16,i)
        r=ceil(i/16);
        c=mod(i,16);
        if c==0
            c=8;
        end
        Hilbert_band=abs(hilbert(HG_Finger5SubBands(band).EnvDeltaWOutMean(:,Ch_num_3(r,c))));
        [R,Lags]=xcorr(Output_DynamicsVel,Hilbert_band,LagN,'normalized');
        plot(Lags,R,'linewidth',1.5)
        xticks('')
        yticks('')
        hold on
        vline(0)
    end
    
end


%% Map of predicted weights for channels across subbands of HG and lags


