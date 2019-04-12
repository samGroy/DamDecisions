function idxRank = MultiRank(f,pref,varargin)
%MultiRank
%rank results from multiobjective PPF

%5/29/18: checked to make sure code can handle a matrix of prefs, each
%sensitive to the frequency of certain decision alternatives used. It
%seems that it'll work based on the current architecture.
%normalize quantities by dividing all values by the maximum value. If
%objective is to minimize the quantity, subtract by 1. The min objective
%assumes that quantities are reported as negative values. The only one
%should be cost. Properties are positive, and hazard is framed in the fitfn
%code as "number of removed hazardous dams".
numM=size(f,2);
ranks=zeros(length(f),1);

%inputs from CritRose may be cell in cell. Flatten these.
if iscell(varargin{1,1})
    varargin=vertcat(varargin{:});
end

%check if user gave preference values.
if ~exist('pref','var')
    rk=1/numM;
    pref=zeros(1,numM)+rk;
    fprintf('You have not defined any preferences. Default uniform preference: %s\n',num2str(rk))
end    
%check if user gave all preferences necessary for assessment.
if size(pref,2)~=numM
    error('Uh oh. You need to make sure every metric has a rank. One rank per metric.')
end

%if user requested specific budget target: (assume cost is at the end of evertything)
if any(strcmpi('leastSquaresCost',varargin))
    ranks=sqrt(sum(power(f(:,end)-pref(end),2),2));
    [~,costRank]=sort(ranks,'ascend');  
    %absolute acceptable cost envelope ($2.5M over- or $5M under-budget)
    %provides fewer options, but stays closer to budget limits
    costRank=costRank(((f(costRank,end)-pref(end))<=2.5e6) & ((f(costRank,end)-pref(end))>=-5e6),:);
    %relative acceptable cost envelope (5% over- or 10% under-budget)
    %provides more options, but generally overruns budget more, esp
    %apparent on rose plots at >100M budgets
%     costRank=costRank(((f(costRank,end)-pref(end))./pref(end)<=0.05) & ((f(costRank,end)-pref(end))./pref(end)>=-0.1),:);
    pref=pref(1:end-1);
    f=f(:,1:end-1);
end

%Check for the ranking method and normalization method selected by the user.
%Options are:
%weightedProduct: weighted product model
%leastSquares: least squares to select a "closest fit" scenario.
%weightedSum: weighted sum model is default if no method provided by user.
h=max(f,[],1);
l=min(f,[],1);
if any(strcmp('normrange',varargin))
    a=(f-l)./(h-l);
    if any(f(1,:)<0)
        a(:,f(1,:)<0)=1-(l(f(1,:)<0)-f(:,f(1,:)<0))./(l(f(1,:)<0)-h(f(1,:)<0));
    end
elseif any(strcmp('stansum',varargin))%standardize by sum
    a=f./sum(f,1);
elseif any(strcmp('stanz',varargin))%standardize by standard deviations
    a=zscore(f);
elseif any(strcmp('none',varargin))%raw values (use only in weighted product!!!!)
    a=f;
else% normax is the default option.
    a=f./h;
    if any(sum(f,1)<0)
        a(:,sum(f,1)<0)=1-(-f(:,sum(f,1)<0)./-l(sum(f,1)<0));
    end
end

if any(strcmpi('leastSquares',varargin)) %least squares useful for seeking scenario from absolute values.
    ranks=sqrt(sum(power(a(:,pref~=0)-pref(pref~=0),2),2));
    [~,idxRank]=sort(ranks,'ascend');
elseif any(strcmpi('weightedProduct',varargin))%nonlinear approach to selection
    pref = prefUnityCheck(pref);
    if any(a(1,:)<0)
        a(:,a(1,:)<0)=-1./(a(:,a(1,:)<0)-1);
    end
    %far as I know, sum of preferences for weighted product don't need to
    %sum to unity? Is that correct?? Probably not...
    ranks=prod(a.^pref,2);
    [~,idxRank]=sort(ranks,'descend'); %must descend, highest ranked value at top
else %weighted sum is default
    pref = prefUnityCheck(pref);
    %     for i=1:numM
    %         ranks=ranks+a(:,i).*pref(i);
    %     end
    ranks=sum(a.*pref,2);
    [~,idxRank]=sort(ranks,'descend');
end

if any(strcmpi('leastSquaresCost',varargin))
    idxRank=idxRank(ismember(idxRank,costRank));
end

function pref = prefUnityCheck(pref)%subfunction, make sure all prefs sum to unity
if any(abs(sum(pref,2)-1)>0.001)
    fprintf('Preferences do not sum to unity. Normalizing.\n')
    pref=pref./sum(pref);
    fprintf('New preferences are: %s\n',num2str(pref))
end
