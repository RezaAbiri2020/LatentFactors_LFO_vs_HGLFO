%Hall paper method
% is it possible to predict the hg-lfo from surrounding lfos
% is it possible to predict the lfo from surrounding hg-lfos

clear all
close all
clc
%% General parameters that should be specified at the begining

% choosing brain area to do analysis on that part and do referencing
% options:
%Selected_Chs;  case 1 all electrodes
%Selected_PMChs; case 2 
%Selected_SMChs; case 3
%Selected_PMSMChs; case 4
%Selected_HandChs; case 5
Brain_part=2;


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
%% loading the processed ECoG data
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\LFO_signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Direct_LFO_Signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Avg_LFO_Signals.mat')



%% 
Nch_Record=256;
% grid layout
Ch_num_1=1:Nch_Record;
Ch_num_2=reshape(Ch_num_1,[16,16]);
% Ch_num_3 is the grid layout
Ch_num_3=rot90(rot90(Ch_num_2));
[R,C]=size(Ch_num_3);

% predict the hg-lfo using surrounding lfos 
for Fi=5
    %HGLFO=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    HGLFO=(HG_Avg_LFO_Signals(Fi).DeltaofEnv);
    R2_hglfo=zeros(R,C);
    for r=2:R-1
        for c=2%:C-1
            Output_Dynamics = HGLFO(:,Ch_num_3(r,c));
            %Output_Dynamics = LFO_signals(Fi).Hilbert(:,Ch_num_3(r,c));
           
            Input_dynamics = [LFO_signals(Fi).Filtered(:,Ch_num_3(r-1,c-1)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r-1,c)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r-1,c+1)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r,c-1)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r,c+1)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r+1,c-1)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r+1,c)),...
                LFO_signals(Fi).Filtered(:,Ch_num_3(r+1,c+1))];
              Figures=1;
            R2=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
            R2_hglfo(r,c)=R2;
        end
    end
    
end

% predict the lfo using surrounding hg-lfos 
for Fi=5
    %HGLFO=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    HGLFO=abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    R2_hglfo=zeros(R,C);
    for r=2:R-1
        for c=2%:C-1
            %Output_Dynamics = LFO_signals(Fi).Hilbert(:,Ch_num_3(r,c));
            Output_Dynamics = HGLFO(:,Ch_num_3(r,c));
            Input_dynamics = [HGLFO(:,Ch_num_3(r-1,c-1)),...
                HGLFO(:,Ch_num_3(r-1,c)),...
                HGLFO(:,Ch_num_3(r-1,c+1)),...
                HGLFO(:,Ch_num_3(r,c-1)),...
                HGLFO(:,Ch_num_3(r,c+1)),...
                HGLFO(:,Ch_num_3(r+1,c-1)),...
                HGLFO(:,Ch_num_3(r+1,c)),...
                HGLFO(:,Ch_num_3(r+1,c+1))];
              Figures=1;
            R2=MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
            R2_hglfo(r,c)=R2;
        end
    end
    
end



%using glmfit
b = glmfit(Input_dynamics,Output_Dynamics);
y_new=[ones(length(Input_dynamics),1),Input_dynamics]*b;
plot(Output_Dynamics)
hold on
plot(y_new)


% predict the hg-lfo using all lfos in that related brain areas 
for ch=[86, 87, 88]
    for Fi=5
        %HGLFO=abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
        HGLFO=(HG_Avg_LFO_Signals(Fi).DeltaofEnv);
        Output_Dynamics = HGLFO(:,ch);
        
        Modified_area = Brain_Area;
        Modified_area(ch) = 0;
        Input_dynamics = LFO_signals(Fi).Filtered(:,Modified_area);
        Figures = 1;
        R2 = MultipleRegFunc(Fi,Input_dynamics,Output_Dynamics,Figures);
        
    end
end
%HighQualityFigs('LFO_decoder_PM_ch87')
