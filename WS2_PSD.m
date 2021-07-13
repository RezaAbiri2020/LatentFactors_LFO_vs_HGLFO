
% showing examples of power spectral density to confirm the LFO activity is
% important during motor control, as well

% looking at power spectrom of channels for choosing LFO band


%% initial parameters

close all;
clear all;
clc;

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

%% loading and breaking raw ECoG data into trials
load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientTwo/ECoGData/ECoG_data.mat');

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

 load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientTwo/ImagingData/TDT_elecs_all.mat')
 
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

%% FFT

% for each finger, plot
figure;
for Fi=1:5
    
    clear P1_MoveAvg
    clear f
    % modification to sync
    ECoGdata_Modified=lfp(1:Nch_Record,trial_start(Fi):trial_stop(Fi));
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
    
    % observe FFT for channels
    k=0;
    for ch=1:256
        if Brain_Area(ch)
            k=k+1;
            data=ECoGdata_2(ch,:);
            L=length(data);
            wq=fft(data);
            P2 = abs(wq/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            % filtering smoothing P1
            P1_MoveAvg(k,:)= movmean(P1,20);
            f(k,:) = Fs*(0:(L/2))/L;
%             figure;
%             plot(f(k,:),P1(:),'b')
%             hold on
%             plot(f(k,:),P1_MoveAvg(k,:),'r')
%             title('Single-Sided Amplitude Spectrum of X(t)')
%             xlabel('f (Hz)')
%             ylabel('|P1(f)|')
            
        end
    end
    
    
    % generating figure with shaded area for 95 CI
    
    % calculating the values for 95 CI
    
    N = size(P1_MoveAvg,1); % Number of Channel In Data Set
    Mean = nanmean(P1_MoveAvg); % Mean Of All channel At Each Value
    SEM = nanstd(P1_MoveAvg)/sqrt(N); % Compute Standard Error Of The Mean Of All Experiments At Each Value Of x
    CI95 = tinv([0.025 0.975], N-1);% Calculate 95% Probability Intervals Of t-Distribution
    CI95_Chs = bsxfun(@times, SEM , CI95(:)); %Calculate 95% Confidence Intervals Of All Experiments At Each Value Of x
    
    % plot for the whole band less than (Fs/2)
%     figure;
%     ciplot(Mean(1,:)+CI95_Chs(1,:),Mean(1,:)-CI95_Chs(1,:),f(1,:),'b')
%     alpha(.4)
%     hold on
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
%     vline(4,'m')
%     xlim([0, 50])
%     set(gca,'FontSize',12)
    %HighQualityFigs('WS2_PSD_1')
    
    
    % zoom in plot for the band less than (50 hz)
    subplot(1,5,Fi)
    ciplot(Mean(1,1:3131)+CI95_Chs(1,1:3131),Mean(1,1:3131)-CI95_Chs(1,1:3131),f(1,1:3131),'b')
    alpha(.4)
    hold on
    title(['PSD for finger: ', num2str(Fi)])
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    vline(1,'m')
    vline(10,'m')
    xlim([0, 50])
    set(gca,'FontSize',12)
    %HighQualityFigs('WS2_PSD_2')
    
end

%HighQualityFigs('WS2_PSD_Fingers')
