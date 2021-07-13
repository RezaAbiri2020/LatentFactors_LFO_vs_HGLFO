
% showing examples of power spectral density to confirm the LFO activity is
% important during motor control, as well

% looking at power spectrom of channels for choosing LFO band

close all;
clear all;
clc;

%% General parameters that should be specified at the begining

% choosing brain area; could have effect on median referencing(?)
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
load('/media/reza/WindowsDrive/ECoGLeapMotion/DataPatientThree/ECoGData/EC176_ecog_kinematics_25Hz.mat');

if 0
    for i=257:260
        subplot(8,8,i-256)
        plot(raw_fullFs{1,1}(1:1000,i))
        
    end
    
    plot(raw_fullFs{1,1}(1:1000,64+2))
    
end

% Final list of bad channels by observation
BadChs=[49, 64+2, 64+50, 192+49];

All_Index=ones(size(raw_fullFs{1,1},2),1);
All_Index(BadChs')=0;
Selected_Chs=logical(All_Index(1:260,1));

%% Determining brain section for analysis e.g.: median ref

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
Nch_Record = length(Selected_Chs);

% for each finger, plot
figure;
for Fi=1:5
    
    clear P1_MoveAvg
    clear f
    ECoGdata_1 = zscore(raw_fullFs{1,Fi+1})';
    
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
    for ch=1:260
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
    %alpha(1)
    hold on
    title(['PSD for finger: ', num2str(Fi)])
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    vline(1,'m')
    vline(10,'m')
    xlim([0, 50])
    set(gca,'FontSize',12)
        
end

%HighQualityFigs('WS3_PSD_Fingers(1hz-10hz)')