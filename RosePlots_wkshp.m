function RosePlots_wkshp(f,rnge,pref,c,ifnm,varargin)
%polarPlot metrics from the MOGA
%used for generating rose plots for workshop maps
%unfortunately this had to be built form a deprecated RosePlots base code,
%and it's pretty messy code.
%inputs:
%f: criteria quantities per scenario (row = scenario, col = criteria)
%rnge: sample multiple scenarios to get a weighted average? No? set to 1.
%pref: either a scalar, indicating the scenario index, or array with
%fractions for each criteria, their sum must equal 1 or it will do that for
%you.
%c: color of rose, in rgb array, fractional [R,G,B], default is black [0,0,0]
%ifnm: name of png save file
%varargin: if you're running MultiRank through RosePlots, here's where you
%add your MultiRank arguments (such as weightedSum, normax, etc)

%gauss weighing? for weighted averaging if range isn't scalar
gauss=1;
if length(rnge)>1
    if gauss
        n=normpdf(1:length(rnge),1,std(1:length(rnge)));
        n=n./sum(n);
    else
        n=zeros(1,rnge)+1/rnge;
    end
else
    n=1;
end

%inputs: f (the fitness scores for your PPF input metrics)
%norm to values to own watershed FOR VISUALIZATION NOT RANKING
h=max(f,[],1);
l=min(f,[],1);
%     a=f./h;
if any(strcmp('normrange',varargin))
    a=(f-l)./(h-l);
    if any(f(1,:)<0)
        a(:,f(1,:)<0)=1-(l(f(1,:)<0)-f(:,f(1,:)<0))./(l(f(1,:)<0)-h(f(1,:)<0)); %indexing only criteria with negative values normalized relative to minimum
    end
else %norm by max (or min for cost) value, default
    a=f./h;
    if any(min(f,[],1)<0)
        a(:,min(f,[],1)<0)=1-(-f(:,min(f,[],1)<0)./-l(min(f,[],1)<0)); %indexing only criteria with negative values normalized relative to minimum
    end
end
%a=(f-l)./(h-l); %normalize metrics in their own watershed
nm=size(a,2); %count the number of metrics

s=360/nm; %spacing between metrics in polar plot, in degrees
%Use standard terms for polar reference frame: rho (magnitude), and theta
%(orientation)
theta=0:0.25:360; %small steps are important in order to produce the nice 'rounded bar chart' look, Johan Rockstrom-style
t=[0 s/2:s:360 360]; %this theta has a degree for each metric, used for interpolation
%for rho, grab a row from a to plot

if length(pref)==1
    scens=pref;
    plotGauge(a,c,scens,nm,theta,t,s,n,ifnm) %subfunction tht plots the rose
else
    if nargin<3
        idxRank=MultiRank(f);
    elseif nargin>=3
        idxRank=MultiRank(f,pref,varargin);%'leastSquares');
        %         idxRank=MultiRank(f,pref);
        c=[0,0,0];
    else
        idxRank=pref;
    end
    idxRank(1)=idxRank(1);
    scens=idxRank(rnge);
    plotGauge(a,c,scens,nm,theta,t,s,n,ifnm)%subfunction tht plots the rose
end

%subfunction tht plots the rose
function plotGauge(a,c,scens,nm,theta,t,s,n,ifnm)

figure %generate figure by default
rn=zeros(1,nm+2);
for k=1:length(scens) %pick the scenario, default is initial scenario: c = length(f)-1;
    %     if k>1 && length(scens)>2
    rn=rn+([a(scens(k),1) a(scens(k),:) a(scens(k),end)].*n(k));
    if k==length(scens)
        r=rn;
    else
        continue
    end
    rho=interp1(t,r,theta,'nearest'); %interpolate the values to the finer resolution theta in order to make a nice looking plot
    rho=[rho 0]; %close the circle
    theta=[theta 0]; %close the circle
    
    % figure
    hold on
    
    %     if k>1
    %set spokes between metrics
    %labels, hardcoded for workshop:
    crit=({'Fish','RiverRec','Storage','Cost','Breach','Property','Electricity','CO2','Heritage','Historic','Town','Aesthetic'});
    ct=1;
    
    %plot spokes separating criteria, one by one, and add the criteria labels
    for i=0:s:360
        xv=cos(deg2rad(i));
        xvc=cos(deg2rad(i+s/2));
        yv=sin(deg2rad(i));
        yvc=sin(deg2rad(i+s/2));
        if length(crit)<ct
            continue
        else
            text(xvc*0.75,yvc*0.75,crit{ct},'horizontalAlignment','center','verticalAlignment','middle','Color','red','FontSize',12,'FontWeight','bold')
            ct=ct+1;
        end
        plot([0 xv],[0 yv],'Color',[0 0 0],'LineWidth',1) 
    end
    %set/plot concentric rings for 25,50,75,100% levels
    for i=0:0.25:1
        xu=i*cos(deg2rad(theta));
        yu=i*sin(deg2rad(theta));
        plot(xu,yu,'Color',[0 0 0],'LineWidth',1)
        
    end
    
    %     end
    
    %now define and plot the gauge
    [x,y]=pol2cart(deg2rad(theta),rho);
    %fill the rose (I used to call it gauge, hence the name)
    gauge=fill(x,y,c,'LineStyle','none');
    plot(x(1:end-1),y(1:end-1),'k') %plot perimeter of rose
    daspect([1 1 1]) %make concentric circles and rose equant
    
    %create transparency gradient, this is not critical, just pretty
    alpha(gauge,0.6) 
    z=(x.^2+y.^2).^(1/2);
    fmin=min(z)-.5;
    fmax=max(z)+0.9;
    z=(z-fmin)./(fmax-fmin).*100;
    gauge.FaceVertexAlphaData=z';
    gauge.FaceAlpha='interp';
    gauge.AlphaDataMapping='direct';
    
end

set(gca,'children',flipud(get(gca,'children'))) %flip everything, for some reason this original code added everything backwards...
axis off %remove axis to keep only background color
hold off %not necessary, but further figur eedits would overwrite everything
ifnm=['D:\FoD\PPF\MCDA-PPF\RosePlots\' ifnm]; %save directory and filename
print(gcf,ifnm,'-dpng','-r300') %print current (gcf=get current figure handle) rose plot to the directory
