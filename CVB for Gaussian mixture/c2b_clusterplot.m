
%/////////////////////////////////////////////// data

h0 = figure;
hold on

RGBblack    = [0,0,0];
RGBdarkBlue = [.08,.17,.55];
RGBdarkCyan = [0,.45,.74];
RGBdarkGreen= [0,.5,0];
RGBdarkViolet=[.49,.18,.56];
RGBdarkRed  = [.6,.2,0];

RGB{1} = RGBblack;       %RGB black
RGB{2} = RGBdarkBlue; %RGB dark Blue
RGB{3} = RGBdarkCyan;   %RGB dark Cyan
RGB{4} = RGBdarkGreen;      %RGB dark Green
RGB{5} = RGBdarkViolet; %RGB dark Violet
RGB{6} = RGBdarkRed;     %RGB dark Red

[Gauss.x,Gauss.y,Gauss.F] = Func_contourGauss([0 0],para.SIGMA);

v = [.003 .05]; %Gauss levels
% v = [.003 .01]; % Gauss levels
% v = [.003 .01 .1]; %Gauss levels

for k=1:setting.K
    contour(Gauss.x+para.MU(1,k),Gauss.y+para.MU(2,k),Gauss.F,v,'--','LineWidth',0.5,'Color',RGB{k});

    Lk = data.L==k;

    plot(data.X(1,Lk),data.X(2,Lk),'+','Color',RGB{k});
end

ax = gca;
plot(ax.XLim,[-1 -1]*Radius + para.offset(2),'k:','Marker','none')
plot(ax.XLim,[ 1  1]*Radius + para.offset(2),'k:','Marker','none')
plot([-1 -1]*Radius + para.offset(1), ax.YLim,'k:','Marker','none')
plot([ 1  1]*Radius + para.offset(1), ax.YLim,'k:','Marker','none')

Ypos = [ 1  1]*Radius + para.offset(2);
plot(1,1,'k.','MarkerSize',15)
annotation('arrow',[0.515 0.69],[0.515 0.7])
ht = text(Ypos(1)/4,Ypos(2)/4 + 0.5,'Radius');
set(ht,'Rotation',40)
set(ht,'FontSize',14)

movegui(h0,'south')
%ax.YAxisLocation = 'origin';

initMU(:,:,1) = [-1 0; 0 1; 1 0; 0 -1]';

%///////////////////////////////////////////////

h1(1) = plot(initMU(1,:),initMU(2,:),'or','MarkerSize',10,'DisplayName','(initial)');

h1(2) = plot(kmean_MU(1,:),kmean_MU(2,:),'sk','MarkerSize',10,'DisplayName','k-means');

h1(3) = plot(emL_MU(1,:),emL_MU(2,:),'xb','MarkerSize',10,'DisplayName','EM_1');

h1(4) = plot(emMU_MU(1,:),emMU_MU(2,:),'^b','MarkerSize',10,'DisplayName','EM_2');

h1(5) = plot(VB_MU(1,:),VB_MU(2,:),'db','MarkerSize',10,'DisplayName','VB');

h1(6) = plot(CVB1_MU(1,:),CVB1_MU(2,:),'hc','MarkerSize',10,'DisplayName','CVB_1');

h1(7) = plot(CVB2_MU(1,:),CVB2_MU(2,:),'pm','MarkerSize',10,'DisplayName','CVB_2');

h1(8) = plot(CVB3_MU(1,:),CVB3_MU(2,:),'*r','MarkerSize',10,'DisplayName','CVB_3');

legend(h1)