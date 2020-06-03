function HighQualityFigs(name)
% Defaults for this blog post

%For generating the target figure size


% 20.41 inch is the size for the width of monitor 
pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
width1 = pos(3);
width=(width1/2560)*20.41;
height1 = pos(4);
height=(height1/2560)*20.41;
% alw = 10;    % AxesLineWidth
% fsz = 50;      % Fontsize
% lw = 30;      % LineWidth
% msz = 30;       % MarkerSize
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

% Save the file as PNG
print(name,'-dpng','-r300');
end