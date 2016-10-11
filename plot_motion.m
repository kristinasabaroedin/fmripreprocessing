function [] = plot_motion(mov)

% [] = plot_motion(mov)
%
% This function will plot the translation and rotation parameters following
% realignment and save them to a pdf file.
%
% ------
% INPUTS
% ------
%
% mov   - an N x 6 matrix, where N = number of time points. Each column
%       represents a different motion parameter as per the following order:
%       1 - x translation; 2 - y translation; 3 - z translation; 4 - pitch;
%       5 - roll; 6 - yaw.
%
% -------
% OUTPUTS
% -------
%
% pdf file contating separate plots for both
%
%==========================================================================

figure(10)
subplot(2,1,1)
plot(mov(:,1));
hold on
plot(mov(:,2),'g');
plot(mov(:,3),'r');
title('translation','fontsize',15,'fontweight','bold')
legend({'x','y','z'})
ylabel('mm')
xlabel('time (volumes)')

subplot(2,1,2)
plot(mov(:,4));
hold on
plot(mov(:,5),'g');
plot(mov(:,6),'r');
title('rotation','fontsize',15,'fontweight','bold');
legend({'pitch','roll','yaw'})
ylabel('degrees')
xlabel('time (volumes)')

set(gcf,'color','white');

saveas(gcf,['prepro_report_motion_',date],'pdf')
