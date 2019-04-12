function RosePlots(f,rnge,pref,c,ifnm,varargin)
%SGR 4/5/19: separate individual metrics within clustered criteria groups
%example: average fish clustered/used in BarrOpt, but now separated to
%visualize changes in each separate species

%gauss weighing? Only if multiple scenarios sampled
if any(strcmp('gauss',varargin))
    n=gaussSample(rnge);
else
    n=1;
end

%inputs: f (the fitness scores for your PPF input metrics)
%norm to own region for viz purposes
h=max(f,[],1);
l=min(f,[],1);
if any(strcmp('normrange',varargin))
    a=(f-l)./(h-l);
    if any(f(1,:)<0)
        a(:,f(1,:)<0)=1-(l(f(1,:)<0)-f(:,f(1,:)<0))./(l(f(1,:)<0)-h(f(1,:)<0));
    end
else %norm by max (or min for cost) value
    a=f./h;
    if any(min(f,[],1)<0)
        a(:,min(f,[],1)<0)=1-(-f(:,min(f,[],1)<0)./-l(min(f,[],1)<0));
    end
end

if length(pref)==1
    scens=pref;
else
    if nargin<4
        idxRank=MultiRank(f);
    elseif nargin>=4
        idxRank=MultiRank(f,pref,varargin);%'leastSquares');
        c=[0,0,0];
    end
    scens=idxRank(range);
end

%nm = number of criteria
nc=size(a,2); 
%s = spacing between criteria in polar plot. 4 criteria = 90 degree spacing.
s=360/nc; 
%theta = criteria quantities interpolated to smooth curve of rose plot
theta=0:0.25:360; %small steps are important in order to produce the nice 'rounded bar chart' look, Johan Rockstrom-style
%t = rose plot coordinates in degrees.
t=[0 s/2:s:360 360]; %this theta has a degree for each metric, used for interpolation

r=zeros(1,nc+2);
for k=1:length(scens)
    %     if k>1 && length(scens)>2
    r=r+([a(scens(k),1) a(scens(k),:) a(scens(k),end)].*n(k));%rho of criteria, including ends
end

rho=interp1(t,r,theta,'nearest'); %interpolate the values to the finer resolution theta in order to make a nice looking plot
rho=[rho 0];
theta=[theta 0];

figure
hold on

plotGrid(s,theta)
plotGauge(theta,rho,c)

axis off
hold off

ifnm=['D:\FoD\PPF\MCDA-PPF\RosePlots\' ifnm];
print(gcf,ifnm,'-dpng','-r300')

function plotGrid(s,theta)
for i=0:s:360
    xv=cos(deg2rad(i));
%     xvc=cos(deg2rad(i+s/2));
    yv=sin(deg2rad(i));
%     yvc=sin(deg2rad(i+s/2));
    %if length(crit)<ct
    %    continue
%     else
        %text(xvc*0.75,yvc*0.75,crit{ct},'horizontalAlignment','center','verticalAlignment','middle','Color','red','FontSize',12,'FontWeight','bold')
        %ct=ct+1;
%     end
    plot([0 xv],[0 yv],'Color',[0 0 0],'LineWidth',1)
end
for i=0:0.25:1
    xu=i*cos(deg2rad(theta));
    yu=i*sin(deg2rad(theta));
    plot(xu,yu,'Color',[0 0 0],'LineWidth',1)
end

function plotGauge(theta,rho,c)
%set gauge area
%now define and plot the gauge
[x,y]=pol2cart(deg2rad(theta),rho);
gauge=fill(x,y,c,'LineStyle','none');
plot(x(1:end-1),y(1:end-1),'k')
%     hold off
daspect([1 1 1])
alpha(gauge,0.6)
%     if k>1
z=(x.^2+y.^2).^(1/2);
fmin=min(z)-.5;
fmax=max(z)+0.9;
z=(z-fmin)./(fmax-fmin).*100;
gauge.FaceVertexAlphaData=z';
gauge.FaceAlpha='interp';
gauge.AlphaDataMapping='direct';
% axis([-1.25 1.25 -1.25 1.25])
set(gca,'children',flipud(get(gca,'children')))

function n=gaussSample(rnge)
%n = gauss distribution sampling pattern
n=normpdf(1:length(rnge),1,std(1:length(rnge)));
n=n./sum(n);