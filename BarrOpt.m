function [x, f, xu, fu, exitflag, output, population, score] = ...
    BarrOpt(BarrIndex,populationSize,varargin)
%1-18-19 update: Organized for general barriers dams and culverts
    %BarrIndex: index of barriers and their downstream neighbor(s)
    %populationSize: number of scenarios you want to solve for
    %simultaneously. Larger number = more accurate PPF, but longer solve
    %time.
    %varargin may contain:  ipop, a population of initial scenarios
    %                       v, contains all criteria names and data
    
    %this version always assumes that like criteria are aggregated (e.g., fish)

tic
if any(strcmp('randomswarm',varargin))
    rs=1;
else
    rs=0;
end

%Set decision variables / find the dams being used to find optimal
%trade-offs % Number of decision variables
numberOfWatersheds=sum(BarrIndex(:,1)>=1e9);
if size(BarrIndex,2)>2 %if island dams exist in wshd, remove the extra data at the end when counting for x
    numberOfVariables = size(BarrIndex,1)-numberOfWatersheds-1;
else %remove the wshd ID at the end from x
    numberOfVariables = size(BarrIndex,1)-numberOfWatersheds;
end

%Always place ipop before v
if length(varargin)>1
    ipop=varargin{1};
    v=varargin{end};
else
    ipop=0;
    v=varargin;
end

%Objective Variables setup: get ranges of the different variables
%skip dependent variables passage, ales, ag, if they exist

%OutNames=VarNames(~ismember(VarNames,[{'passage'},{'ag'}]));
% end

%Report OVs to user
fprintf('Number of Criteria = %i:\n',length(v))
for i=1:length(v)
    fprintf('%i) %s\n',i,v{1,i})
end

%Norm metrics to their max (or min....)
%Necessary because cross-unit values are summed for the cluster of 5
%criteria: example; hydropower and reservir storage
%Do this only for cases where cross-metrics are summed for their criterion
%Diadromous Fish Habitat. Problem is I don't know the max value because it
%requires the network calc
%first copy v to un-normalize later
v_cpy=v;
if any(strcmp('Diadromous Fish Habitat',v(1,:)))
    try
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.BkTHabA=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.BkTHabA./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.BkTHabA(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
    try
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.AlRm_5d_sm=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.AlRm_5d_sm./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.RhHabA(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.RhHabA=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.RhHabA./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.RhHabA(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
    try
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.SalmoHabA=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.SalmoHabA./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.SalmoHabA(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
    try
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.ShadHabA=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.ShadHabA./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.ShadHabA(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
    try
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.SLamHabA=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.SLamHabA./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.SLamHabA(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
    try
        v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.EelHabA=...
            v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.EelHabA./...
            sum(v{2,strcmp('Diadromous Fish Habitat',v(1,:))}.EelHabA(numberOfVariables+1:numberOfVariables+numberOfWatersheds));
    catch
    end
end
%Dam Utilities
if any(strcmp('Dam Utilities',v(1,:)))
    try
        v{2,strcmp('Dam Utilities',v(1,:))}.Cpcty_kW=...
            v{2,strcmp('Dam Utilities',v(1,:))}.Cpcty_kW./...
            sum(v{2,strcmp('Dam Utilities',v(1,:))}.Cpcty_kW(1:numberOfVariables));
    catch
    end
    try
        v{2,strcmp('Dam Utilities',v(1,:))}.RsvrVolEst=...
            v{2,strcmp('Dam Utilities',v(1,:))}.RsvrVolEst./...
            sum(v{2,strcmp('Dam Utilities',v(1,:))}.RsvrVolEst(1:numberOfVariables));
    catch
    end
end
%Recreational Boating:problem is I don't know what the max value is because
%it relies on the network calc.
if any(strcmp('Recreational Boating',v(1,:)))
    try
        v{2,strcmp('Recreational Boating',v(1,:))}.FRRU=...
            v{2,strcmp('Recreational Boating',v(1,:))}.FRRU./...
            sum(v{2,strcmp('Recreational Boating',v(1,:))}.FRRU(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
    try
        v{2,strcmp('Recreational Boating',v(1,:))}.FLRU=...
            v{2,strcmp('Recreational Boating',v(1,:))}.FLRU./...
            sum(v{2,strcmp('Recreational Boating',v(1,:))}.FLRU(numberOfVariables+1:numberOfVariables+1+numberOfWatersheds-1));
    catch
    end
end

%if random swarm
if rs
    f=nan(populationSize,length(v));
    for i=1:populationSize
        x=round(rand(1,numberOfVariables));
        f(i,:)=BarrOpt_fitfn(x,BarrIndex,v);
    end
    exitflag=[];
    output=[];
    population=[];
    score=[];
else
    fitFn = @(x)BarrOpt_fitfn(x,BarrIndex,v);
    
    %double OV vector necessary when decision is not binary
    %     if v{end}(8)>2
    %         options = gaoptimset('PopulationType','doubleVector',...%'bitstring',...
    %             'PopulationSize',populationSize,...
    %             'Generations',500,...
    %             'CrossoverFraction',0.8,...
    %             'UseParallel',true);%*****************
    %         if length(ipop)>1
    %             options.InitialPopulation = ipop;
    %         end
    %         lb=zeros(1,numberOfVariables);
    %         ub=ones(1,numberOfVariables).*(v{end}(8)-0.001);
    %         [x, f, exitflag, output, population, score] = gamultiobj(fitFn,...
    %             numberOfVariables, [], [], [], [], lb, ub,options);
    %     else
    options = gaoptimset('PopulationType','bitstring',...
        'PopulationSize',populationSize,...
        'Generations',500,...
        'CrossoverFraction',0.8,...
        'UseParallel',true);%*****************
    if length(ipop)>1
        options.InitialPopulation = ipop;
    end
    %         options.PlotFcns = @gaplotpareto;
    [x, f, exitflag, output, population, score] = gamultiobj(fitFn,...
        numberOfVariables, [], [], [], [], [], [],options);
    ub=1;
    %     end
end

%extrema may not be pareto efficient. Add manually.
extra_x = zeros(2,numberOfVariables);
extra_x(1,:) = ub;
x(end+1:end+2,:)=extra_x;
% if ismember('ag',VarNames)
f(end+1,:) = BarrOpt_fitfn(extra_x(1,:),BarrIndex,v);
f(end+1,:) = BarrOpt_fitfn(extra_x(2,:),BarrIndex,v);

%now time to keep only unique PPF points
[~,id]=unique(f,'rows');
xu=x(id,:);
fu=zeros(size(xu,1),size(f,2));
for i=1:size(xu,1)
    fu(i,:) = -BarrOpt_fitfn(xu(i,:),BarrIndex,v_cpy);
end

f=-f; %think positively.

%convert f matrix to f table
% fields=cell(1,size(v,2));
% for i=1:size(v,2)
%     fields(i)={strrep(v{1,i},' ','_')};
% end
% f=array2table(f);
% f.Properties.VariableNames=fields;
toc
