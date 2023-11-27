clear all;
clc;

mean_tau = 2;
x = 0:0.01:2*mean_tau;

i=0;
NE = [1,10,100];


hf = figure;

for shape = NE
rate = shape/mean_tau;
pdf = ((rate^shape) * (x.^(shape -1)) .* exp(-rate*x)/gamma(shape));
i=i+1;
ha(i) = subplot(length(NE),1,i)
plot(x,pdf,'LineStyle','-','LineWidth',2);hold on;yticks([]);
set(gca,'FontSize',20)
yticks([0 0.3 0.6 0.9 1.2 1.5 1.8 2.1]);
line([mean_tau, mean_tau], [0, max(pdf)], 'Color', [0.1,0.1,0.1],'LineStyle','--',LineWidth=2);




CV = 1/sqrt(shape);
str = strcat('cv = ',num2str(CV),' N_E =', num2str(shape));

%annotation('textbox',dim,'String',str,'FitBoxToText','on');

%ha(i).Position = ha(i).Position + [0 0 0 0.5]; 

end




pos = get(ha, 'position');
dim = cellfun(@(x) x.*[1 1 0.5 0.5], pos, 'uni',0);

% 2 points after decimal


for i = 1: length(NE)

CV = 1/sqrt(NE(i));
annotation(hf, 'textbox', dim{i}, 'String',  strcat('CV = ','  ',num2str(CV,2), '; ',' N_E =   ',num2str(NE(i)) )  ,'Position', ha(i).Position + [0.45 0 0 0],'EdgeColor','none', 'FontSize',20,'verticalalignment', 'top','FitBoxToText','on') ;

end


han=axes(hf,'visible','off'); 
%han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'PDF of \tau ');
xlabel(han,'\tau (hours)');
set(gca,'FontSize',20)



% saving figure

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 28]); %x_width=10cm y_width=15cm
saveas(gcf,'cv-ne.png')