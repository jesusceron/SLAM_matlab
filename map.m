% Map shown in the "Trajectory ZUPT-KF" figure

% Map shown in the main UI

rectangle('Position',[5.5 2 12.6 6.6]) % Whole the seminar room
% %stairs
% rectangle('Position',[-0.65 0.3125 1.3 0.3125])
% rectangle('Position',[-0.65 0.6250 1.3 0.3125])
% rectangle('Position',[-0.65 0.9375 1.3 0.3125])
% rectangle('Position',[-0.65 1.2500 1.3 0.3125])
% rectangle('Position',[-0.65 1.5625 1.3 0.3125])
% rectangle('Position',[-0.65 1.8750 1.3 0.3125])
% rectangle('Position',[-0.65 2.1875 1.3 0.3125])
% 
% rectangle('Position',[-0.65 2.5000 1.3 0.3125])
% rectangle('Position',[-0.65 2.8125 1.3 0.3125])
% rectangle('Position',[-0.65 3.1250 1.3 0.3125])
% rectangle('Position',[-0.65 3.4375 1.3 0.3125])
% rectangle('Position',[-0.65 3.7500 1.3 0.3125])
% rectangle('Position',[-0.65 4.0625 1.3 0.3125])
% rectangle('Position',[-0.65 4.3750 1.3 0.3125])
% rectangle('Position',[-0.65 4.6875 1.3 0.3125])
% rectangle('Position',[-0.65 5 1.3 0.3125])

%zones without use:
rectangle('Position',[5.5 2 4.3 2.9], 'FaceColor','k')
rectangle('Position',[5.5 6.3 1.5 2.2], 'FaceColor','k')

% Room and living room
rectangle('Position',[11.8 2 0.8 1.8], 'FaceColor','k') % dining table
rectangle('Position',[13.9 2 0.6 1.4], 'FaceColor','k') % wall btw dining and living room
rectangle('Position',[15.1 2 1.4 0.6], 'FaceColor','k') % Living table
rectangle('Position',[13.9 4.5 4.2 0.4], 'FaceColor','k') % wall btw room and living room
rectangle('Position',[16.5 5.9 1.6 1.6], 'FaceColor','k') % bed

% Kitchen
rectangle('Position',[8.1 5.9 3.1 0.6], 'FaceColor','k') % Wall btw bathroom-kitchen and main corridor
rectangle('Position',[12.5 5.9 1.4 0.6], 'FaceColor','k') % Wall btw kitchen and corridor
rectangle('Position',[9.5 6.5 0.8 1.6], 'FaceColor','k') % Wall btw bathroom and kitchen
rectangle('Position',[13.1 6.5 0.8 1.6], 'FaceColor','k') % Wall btw kitchen and room
rectangle('Position',[10.3 7.9 2.8 0.6], 'FaceColor','k') % Tables of the kitchen

% % Plot fixed Beacons
% plot(13, 4.4,'^','markerfacecolor','y','MarkerSize',10)  % Room
% %viscircles([13, 4.4], 4,'Color','y','LineStyle','--');
% plot(8, 6.2,'^','markerfacecolor','k','MarkerSize',10)  % Kitchen
% %viscircles([8, 6.2], 4,'Color','k','LineStyle','--');
% plot(4.5, 6.8,'^','markerfacecolor',[0.4940, 0.1840, 0.5560],'MarkerSize',10)  % Bathroom
% %viscircles([4.5, 6.8], 4,'Color',[0.4940, 0.1840, 0.5560],'LineStyle','--');
% plot(8.4, 1.3,'^','markerfacecolor',[0.3010, 0.7450, 0.9330],'MarkerSize',10)  % Dining table
% %viscircles([8.4, 1.3], 4,'Color',[0.3010, 0.7450, 0.9330],'LineStyle','--');
% plot(12.3, 1,'^','markerfacecolor',[0.4660, 0.6740, 0.1880],'MarkerSize',10)  % Living room
% %viscircles([12.3, 1], 4,'Color',[0.4660, 0.6740, 0.1880],'LineStyle','--');
% 
% %plot(2, 3.6,'^','markerfacecolor',[0, 0.4470, 0.7410],'MarkerSize',10)  % Door
% %viscircles([2, 3.6], 4,'Color',[0, 0.4470, 0.7410],'LineStyle','--');
% 
% 
% % Plot 'activity' Beacons
% plot(2, 3.6,'s','markerfacecolor','r','MarkerSize',8)  % Door
% plot(5.2, 6.4,'s','markerfacecolor','g','MarkerSize',8)  % Toiler lid
% plot(3.7, 6.4,'s','markerfacecolor','b','MarkerSize',8)  % Broom
% plot(9.6, 5.4,'s','markerfacecolor','c','MarkerSize',8)  % Pitcher
% plot(5.5, 4.4,'s','markerfacecolor','m','MarkerSize',8)  % Hair brush
