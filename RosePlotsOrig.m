function RosePlotsOrig(f,rnge,pref,c,ifnm,varargin)
%polarPlot metrics from the PPF analysis
%used for SPC: scenario pairwise comparison
cd('D:\FoD\PPF\Figures\fingerprints')
%gauss weighing?
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
if ~iscell(f)
    %norm to own watershed
    h=max(f,[],1);
    l=min(f,[],1);
%     a=f./h;
    if any(strcmp('normrange',varargin))
        a=(f-l)./(h-l);
        if any(f(1,:)<0)
            a(:,f(1,:)<0)=1-(l(f(1,:)<0)-f(:,f(1,:)<0))./(l(f(1,:)<0)-h(f(1,:)<0));
        end
    else %norm by max (or min for cost) value
        a=f./h;
        if any(f(1,:)<0)
            a(:,f(1,:)<0)=1-(-f(:,f(1,:)<0)./-l(f(1,:)<0));
        end
    end
    %a=(f-l)./(h-l); %normalize metrics in their own watershed
    nm=size(a,2); %count the number of metrics

else
    %norm to entire region
    nm=size(f{2,2},2);
    hn=zeros(size(f,2)-1,nm); h=zeros(1,nm); %l=zeros(1,nm);
    for i=2:1:size(f,2)
        try
            hn(i-1,:)=max(f{2,i},[],1);
            %             h(hn(i)>h)=hn(i)(hn(i)>h); %Commout to normalize for
            %             watershed, rather than for entire study region.
            %             h=h+max(f{2,i},[],1);
            %             l=l+min(f{2,i},[],1);
        end
    end
end

s=360/nm; %spacing between metrics in polar plot, in degrees
%Use standard terms for polar reference frame: rho (magnitude), and theta
%(orientation)
theta=0:0.25:360; %small steps are important in order to produce the nice 'rounded bar chart' look, Johan Rockstrom-style
t=[0 s/2:s:360 360]; %this theta has a degree for each metric, used for interpolation
%for rho, grab a row from a to plot


if iscell(f)
    for i=2:1:size(f,2)
        if ~isempty(f{2,i})
            %             a=(f{2,i}-l)./(h-l);
            a=f{2,i}./hn(i-1,:);
            %a=f{2,i}./h;
            col=f{3,i};
            na=f{1,i};
            if nargin<3
                idxRank=MultiRank(f{2,i});
            else
                idxRank=MultiRank(f{2,i},pref);
            end
            %             scens=[length(a)-1; idxRank(1:rnge)];
            scens=idxRank(rnge);
            plotGauge(a,col,scens,nm,theta,t,s,n,na)
        end
    end
else
    if nargin<3
        idxRank=MultiRank(f);
    elseif nargin>=3
        idxRank=MultiRank(f,pref,varargin);%'leastSquares');
%         idxRank=MultiRank(f,pref);
        col=[0,0,0];
    else
        idxRank=pref;
        col=c;
    end
    idxRank(1)=idxRank(1);
    scens=idxRank(rnge);
    plotGauge(a,col,scens,nm,theta,t,s,n,ifnm)
end


function plotGauge(a,col,scens,nm,theta,t,s,n,ifnm,na)

figure
% scens=[length(f)-1 1];
colors=[0.2 0.2 0.2;col];%distinguishable_colors(length(scens)-1,{'w','k'});
rn=zeros(1,nm+2);
%colors=[0 0 0;colors];
for k=1:length(scens) %pick the scenario, default is initial scenario: c = length(f)-1;
    c=scens(k);
    %     if k>1 && length(scens)>2
    rn=rn+([a(c,1) a(c,:) a(c,end)].*n(k));
    if k==length(scens)
        r=rn;
    else
        continue
    end
    %     else
    %         r=[a(c,1) a(c,:) a(c,end)];
    %     end
    rho=interp1(t,r,theta,'nearest'); %interpolate the values to the finer resolution theta in order to make a nice looking plot
    rho=[rho 0];
    theta=[theta 0];
    
    % figure
    hold on
    
%     if k>1
        %set spokes between metrics
        %labels:
        crit=({'Power','Storage','property','cost','hazard','fish','river rec'});
        ct=1;
        si=s;
        for i=0:s:360
            xv=cos(deg2rad(i));
            xvc=cos(deg2rad(i+si/2));
            yv=sin(deg2rad(i));
            yvc=sin(deg2rad(i+si/2));
            if length(crit)<ct
                continue
            else
                text(xvc*0.75,yvc*0.75,crit{ct},'horizontalAlignment','center','verticalAlignment','middle','Color','red','FontSize',16,'FontWeight','bold')
                ct=ct+1;
            end
            plot([0 xv],[0 yv],'Color',[0 0 0],'LineWidth',1)
        end        
        %set gauge area
        for i=0:0.25:1
            xu=i*cos(deg2rad(theta));
            yu=i*sin(deg2rad(theta));
            plot(xu,yu,'Color',[0 0 0],'LineWidth',1)
            
        end

%     end
    
    %now define and plot the gauge
    [x,y]=pol2cart(deg2rad(theta),rho);
%     if k==1, co=1; else, 
        co=2;
    gauge=fill(x,y,colors(co,:),'LineStyle','none');
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
    
    %     end
end

% axis([-1.25 1.25 -1.25 1.25])
set(gca,'children',flipud(get(gca,'children')))
axis off
hold off
% if nargin>8
%     print(gcf,sprintf('LINKED_%s_equal',na),'-dpng','-r300')
% else
print(gcf,ifnm,'-dpng','-r300')
% end
% close