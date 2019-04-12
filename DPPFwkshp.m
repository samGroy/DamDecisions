function [x, f, OutNames, exitflag, output, population, score] = ...
    DPPFwkshp(DamIndex,populationSize,blocklist,varargin)
%3-4-19 version of DPPF tailored purely for Emma's workshop. Requires
%doubleVector argument variables to cover multiple decision alternatives.
%These doubleVector AV's are floored in DPPFwkshp_fitfn to fit the
%integer-based model for selecting decision alternatives. The MOGA seems to
%accept this approach, but it is slower than binary. There is no built in
%option for integers. It's much faster/mroe efficient than the summed
%concatenated binary approach.

tic %time it

%Set decision variables / find the dams being used to find optimal
%trade-offs % Number of decision variables
numberOfWatersheds=sum(DamIndex(:,1)>=1e9); %All wshd nodes have DID codes >=1e9
if size(DamIndex,2)>2 %if island dams exist in wshd, remove the extra data at the end when counting for x
    numberOfVariables = size(DamIndex,1)-numberOfWatersheds-1; %island dam DID (and importantly, habitat data) added at end
else %remove the wshd ID at the end from x
    numberOfVariables = size(DamIndex,1)-numberOfWatersheds;
end


ipop=0;
for i=1:length(varargin)
    if iscell(varargin{i})
        VarNames=varargin{i}(1,1:end);
        v=varargin{i}(2,1:end);
    elseif isa(varargin{i},'double')
        ipop=varargin{i};
    end
end
varargin=[];
%Objective Variables setup: get ranges of the different variables
%skip dependent variables passage, ales, ag, if they exist

% If aggregating, check which parts are included
if ismember('ag',VarNames)
    %     fi=zeros(length(varargin),1);
    %fish aggregation
    if ~ismember('alewifeBlueback',VarNames)
        v{strcmp('ag',VarNames)>0}(1)=0;
        %     else
        %         fi=fi+ismember(VarNames,'alewifeBlueback');
    end
    if ~ismember('shad',VarNames)
        v{strcmp('ag',VarNames)>0}(2)=0;
        %     else
        %         fi=fi+ismember(VarNames,'shad');
    end
    if ~ismember('salmon',VarNames)
        v{strcmp('ag',VarNames)>0}(3)=0;
        %     else
        %         fi=fi+ismember(VarNames,'salmon');
    end
    if ~ismember('lakeRecreation',VarNames)
        v{strcmp('ag',VarNames)>0}(4)=0;
        %     else
        %         r=r+ismember(VarNames,'lakeRecreation');
    end
    if ~ismember('riverRaft',VarNames) || ~ismember('riverCanoe',VarNames)
        v{strcmp('ag',VarNames)>0}(5)=0;
        v{strcmp('ag',VarNames)>0}(7)=1;
        %     else
        %         r=r+ismember(VarNames,'riverRecreation');
    end
end

%User may want to merge services that are similar, to simplify
%results/analysis
if v{strcmp('ag',VarNames)>0}(7)>0 %merge fish species as fish biomass
    OutNames=VarNames(~ismember(VarNames,[{'shad'},{'salmon'},{'alewifeBlueback'},{'passage'},{'ag'}]));
    if sum(v{strcmp('ag',VarNames)}(1:3)>0)
        OutNames{end+1}='fishBiomass';
    end
    if v{strcmp('ag',VarNames)>0}(7)>1 %merge river recreation by small and large craft
        OutNames=OutNames(~ismember(OutNames,[{'riverCanoe'},{'riverRaft'}]));
        if v{strcmp('ag',VarNames)>0}(7)>2
            OutNames=OutNames(~ismember(OutNames,{'lakeRecreation'}));
            OutNames{end+1}='recreation';
        else
            OutNames{end+1}='riverRecreation';
        end
        %         if sum(varargin{strcmp('ag',VarNames)}(4:5)>0)
        %             OutNames{end+1}='recreation';
        %         end
    end
else
    OutNames=VarNames(~ismember(VarNames,[{'passage'},{'ag'}]));
end

OVcount=length(OutNames);

%List the objective variable (OV) names, and manage the index values where you
%expect to write them out.
%Report OVs to user
fprintf('OV count = %i:\n',OVcount)
for i=1:length(OutNames)
    fprintf('%i) %s\n',i,OutNames{i})
end

%Here define the fitness function, used to calculate criteria quantities,
%that are then used within MOGA ranking algorithm to determine if they are
%efficient relative to scenarios in the pool.
fitFn = @(x)DPPFwkshp_fitfn(x,DamIndex,v,VarNames,OVcount);

%Some trial and error shows populationSize is the only option variable that
%improves/completes representation of the pareto frontier. The second
%important thing is to manually re-run the MOGA algo several times with the
%prior population (x). For some reason, sometimes the MOGA just needs a fresh start...
options = gaoptimset('PopulationType','doubleVector',...%'bitstring',... %Very important difference between two (binary) and more than two (integer) decision alternatives
    'PopulationSize',populationSize,...
    'Generations',500,...%MOGA often quits between 100-200, def worth overestimating though
    'CrossoverFraction',0.8,...%the high ranked/elite population fraction used in the crossover function, to make children scenarios
    'UseParallel',true);%*****************%Up to 4x faster than serial code. Do it if you have multiple cores!

%If user provides an initial population of scenarios, use it. Time to
%solution and level of efficiency is often much greater than starting from
%scratch with pseudorandom values.
if length(ipop)>1
    options.InitialPopulation = ipop;
end

%upper and lower bounds of argument variables. Only important for
%doubleVector/integer simulations with >2 decision alts
lb=zeros(1,numberOfVariables);
lb(blocklist==1)=1;
ub=ones(1,numberOfVariables).*(v{end}(8)-0.001);
ub(blocklist==1)=1;

%This is the call for the MOGA, genetic algorithm with multiobjective option
[x, f, exitflag, output, population, score] = gamultiobj(fitFn,...
    numberOfVariables, [], [], [], [], lb, ub,options);

%At the end, be sure that the MOGA found the extrema. If it didn't, just
%run them manually and append to the criteria quantities (f) and the
%population (x)
extra_x = zeros(2,numberOfVariables);
extra_x(1,:) = ub;
x(end+1:end+2,:)=extra_x;

f(end+1,:) = DPPFwkshp_fitfn(extra_x(1,:),DamIndex,v,VarNames,OVcount);
f(end+1,:) = DPPFwkshp_fitfn(extra_x(2,:),DamIndex,v,VarNames,OVcount);

%MOGA works to minimize variables, so OV signs are flipped the whole time.
%Flip back for post-analysis (e.g., positive values for fish, hydropower, etc., negative
%values for cost.
f=-f;
toc %report how long it took
