
% choose finger
Fi=5;
amp=abs(hilbert(HG_DeltaofAvgHG(Fi).Finger));
End=length(HG_DeltaofAvgHG(Fi).Finger);

%amp=Delta_Fingers(Fi).Hilbert;
%End=length(Delta_Fingers(Fi).Hilbert);

BrainArea=Brain_Area;

%% set a figure for all plots 
Ch_num_1=1:256;
Ch_num_2=reshape(Ch_num_1,[16,16]);
Ch_num_3=rot90(rot90(Ch_num_2));

Central_Sulcus=[0,3.5;2,3.5;2.5,4;4.5,5;4.5,6;5.5,7;7.5,8;10.5,9;10.5,10;11.5,11;13.5,12;14,12.5;15,12.5;17,12.5];
set(figure,'position',[100,100,850,850]);
%suptitle('PC1 for ERP (8PCs~90%); Band: Mu Rhythms; Analysis Area: M1+S1')
suptitle(['Hilbert; finger: ',num2str(Fi)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finger: for Delta of Avg-HG

i=1;
All_W=amp(i,:);
All_W=reshape(All_W,[16,16]);
All_W=rot90(rot90(All_W));

All_handle=imagesc(All_W);
% Adjusting the transparency for individual channel image
BrainArea=double(BrainArea);
BrainArea=reshape(BrainArea,[16,16]);
BrainArea=rot90(rot90(BrainArea));
alpha(BrainArea)

All_ax=gca;
All_q=0.08; %0.08; %0.001; % % %0.01; Delta_q=5; %0.2; % % %0.7;
All_ax.CLim=[-All_q All_q];

%colorbar 
hold on;
for col=1:16
    for row=1:16 
        plot(col,row,'ok')
        hold on
    end 
end

hold on;
plot(Central_Sulcus(:,2),Central_Sulcus(:,1),'-k','linewidth',3)
text(12,4,'M1','FontSize',40)
text(2,12,'S1','FontSize',40)
All_time=text(1,17,[num2str(0),'     Sample'],'FontSize',12);
%title ('Delta of Avg-HG')
title ('HG-LFO')
axis off 

%% plotting the figure while updating the handle

vobj=VideoWriter(['Hilberts_HGLFO_WholeGrid_F',num2str(Fi)], 'Motion JPEG AVI');
%vobj=VideoWriter('PC1_MuRythm', 'Archival');

vobj.FrameRate=60;
%vobj.Quality=75;


open(vobj);

for i=1:5:End 
    
    %%%% HG 
    All_W=amp(i,:);
    All_handle=UpdatePlot(All_handle,All_W,BrainArea);
    All_ax.CLim=[-All_q All_q];
    set(All_time,'String',[num2str(i),'   Sample'],'FontSize',12)
         
    drawnow
    %%%%%%%%%%%%%%%%%%%
    F=getframe(gcf);
    writeVideo(vobj, F);
     
end

close(vobj)

%% update plot; function
function handle=UpdatePlot(handle,Weights,BrainArea)
W2=reshape(Weights,[16,16]);
W3=rot90(rot90(W2));
handle.CData=W3;
handle.AlphaData=BrainArea;

end 


