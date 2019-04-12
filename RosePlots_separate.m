function RosePlots_separate(f,rnge,pref,c,ifnm,varargin)
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

%determine the scenario(s) and criteria quantities used to generate plot.
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


nc=size(a,2); %nc = number of criteria
%spacing between criteria in rose plot. this is a hack.
dm=4;%3;
all=1;
if dm==4
    s_top=[0 45 90 90 90.01 105 120 135 150 165 180];% 180 180.1 225 270 270 270.1 360]; %4D
    s_bottom=[180.1 225 270 270 270.1 360];
    if all==1
        s=[s_top s_bottom]; %4D
        t=[0 diff(s)./2+s(1:end-1) 360];
        r=sum([a(scens,1) a(scens,1:2) a(scens,2) a(scens,3) a(scens,3:8) a(scens,8) a(scens,9) a(scens,9:10) a(scens,10) a(scens,11) a(scens,end)].*n,1);%rho of criteria, including ends
        theta=0:0.25:360;%%%theta = criteria quantities interpolated to smooth curve of rose plot
    elseif all==2
        s=s_top;
        t=[0 diff(s)./2+s(1:end-1) 180];
        r=sum([a(scens,1) a(scens,1:2) a(scens,2) a(scens,3) a(scens,3:8) a(scens,8)].*n,1);%rho of criteria, including ends
        theta=0:0.25:180;%%%theta = criteria quantities interpolated to smooth curve of rose plot
    elseif all==3
        s=s_bottom;
        t=[180 diff(s)./2+s(1:end-1) 360];
        r=sum([a(scens,9) a(scens,9:10) a(scens,10) a(scens,11) a(scens,11) a(scens,end)].*n,1);
        theta=180:0.25:360;%%%theta = criteria quantities interpolated to smooth curve of rose plot
    end
elseif dm==3 %TODO
    s_top=[0 60 120 120 120.1 140 160 180 200 220 240];% 180 180.1 225 270 270 270.1 360]; %4D
    s_bottom=[240 360];
    if all==1
        s=[s_top s_bottom]; %4D
        sadd=6;
        t=[0 diff(s)./2+s(1:end-1) 360];
    elseif all==2
        s=s_top;
        sadd=4;
        t=[0 diff(s)./2+s(1:end-1) 240];
    elseif all==3
        s=s_bottom;
        sadd=2;
        t=[240 diff(s)./2+s(1:end-1) 360];
    end
%     t=[0 diff(s)./2+s(1:end-1) 360];%[0 diff(s)./2+s(1:end-1) 240 240.1];%360];[240 360];%
    r=zeros(1,nc+sadd); %zeros(1,sadd); %initialize rho values of criteria (fractional based on normalization preference)
    for k=1:length(scens) %sample across range of scenarios?
        r=r+([a(scens(k),1) a(scens(k),1:2) a(scens(k),2) a(scens(k),3) a(scens(k),3:8) a(scens(k),8) a(scens(k),9) a(scens(k),9:10) a(scens(k),10) a(scens(k),11) a(scens(k),11) a(scens(k),end)].*n(k));%rho of criteria, including ends
    end
end
% s=[0 60 120 120 120.1 140 160 180 200 220 240 240 240.1 360];%[240 360];%[0 60 120 120 120.1 140 160 180 200 220 240];% 240 240.1 360]; %
% sadd=8;%6;%2;%4;%6;
% s=360/nc; %s = spacing between criteria in polar plot. 4 criteria = 90 degree spacing.
% t=[0 s/2:s:360 360]; %t = rose plot coordinates in degrees, plus start/endpoints.
% t=[0 diff(s)./2+s(1:end-1) 360];%[0 diff(s)./2+s(1:end-1) 240 240.1];%360];[240 360];%
% t=[0 diff(s)./2+s(1:end-1) 360];%[0 diff(s)./2+s(1:end-1) 240 240.1];%360];[240 360];%
% r=zeros(1,nc+sadd); %zeros(1,sadd); %initialize rho values of criteria (fractional based on normalization preference)
% for k=1:length(scens) %sample across range of scenarios?
%     r=r+([a(scens(k),1) a(scens(k),1:2) a(scens(k),2) a(scens(k),3) a(scens(k),3:8) a(scens(k),8) a(scens(k),9) a(scens(k),9:10) a(scens(k),10) a(scens(k),11) a(scens(k),11) a(scens(k),end)].*n(k));%rho of criteria, including ends
% end
% for k=1:length(scens) %sample across range of scenarios?
% %       r=r+([a(scens(k),end) a(scens(k),end)].*n(k));%rho of criteria, including ends  
% %      r=r+([a(scens(k),1) a(scens(k),1:2) a(scens(k),2) a(scens(k),3) a(scens(k),3:8) a(scens(k),8) a(scens(k),end)].*n(k));%rho of criteria, including ends  
% %      r=r+([a(scens(k),1) a(scens(k),1:2) a(scens(k),2) a(scens(k),3) a(scens(k),3:8) a(scens(k),8) a(scens(k),9) a(scens(k),9) a(scens(k),end)].*n(k));%rho of criteria, including ends
%       r=r+([a(scens(k),1) a(scens(k),1:2) a(scens(k),2) a(scens(k),3) a(scens(k),3:8) a(scens(k),8) a(scens(k),9) a(scens(k),9:10) a(scens(k),10) a(scens(k),11) a(scens(k),11) a(scens(k),end)].*n(k));%rho of criteria, including ends
% end

rho=interp1(t,r,theta,'nearest'); %interpolate the values to the finer resolution theta in order to make a nice looking plot
rho=[rho 0];
theta=[theta 0];

%figure
hold on

plotGrid(s,theta,dm)
plotGauge(theta,rho,c)

axis off
%hold off
ifnm=['D:\FoD\PPF\MCDA-PPF\RosePlots\' ifnm];
% set(gcf,'color','none')
print(gcf,ifnm,'-dpng','-r300')
%set(gcf,'color','w')

function plotGrid(s,theta,dm)
%s=unique(round(s));
for i=1:length(s)
    xv=cos(deg2rad(s(i)));
    yv=sin(deg2rad(s(i)));
    
    if dm==4
        if (ismember(s(i),0:90:360) && dm==4) || (ismember(s(i),0:120:360) && dm==3)
            plot([0 xv],[0 yv],'Color',[0 0 0 1],'LineWidth',5)
        else
            plot([0 xv],[0 yv],'Color',[0 0 0 .33],'LineWidth',2)
        end
    end
end
for i=0:0.25:1
    xu=i*cos(deg2rad(theta));
    yu=i*sin(deg2rad(theta));
    plot(xu,yu,'Color',[0 0 0 1],'LineWidth',3)
end

function plotGauge(theta,rho,c)
%set gauge area
%now define and plot the gauge
[x,y]=pol2cart(deg2rad(theta),rho);
gauge=fill(x,y,c,'LineStyle','none');
plot(x(1:end-1),y(1:end-1),'k')
%     hold off
daspect([1 1 1])
%alpha(gauge,0.6)
% z=(x.^2+y.^2).^(1/2);
% fmin=min(z)-.5;
% fmax=max(z)+0.9;
% z=(z-fmin)./(fmax-fmin).*100;
% gauge.FaceVertexCData=z';
% gauge.FaceVertexAlphaData=z';
% gauge.FaceAlpha='interp';
% gauge.AlphaDataMapping='direct';
% axis([-1.25 1.25 -1.25 1.25])
set(gca,'children',flipud(get(gca,'children')))

function n=gaussSample(rnge)
%n = gauss distribution sampling pattern
n=normpdf(1:length(rnge),1,std(1:length(rnge)));
n=n./sum(n);