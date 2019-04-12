function [DamIndex, v] = ...
    DPPFwkshp_prep(watershed, basedir, varargin) %%%%%,costRemove, lakeRecreationRm, alewifeBluebackRm
%Use this script to select the watershed/dams you want to PPF
%WARNING: this function is messy as hell and not up to date with my current
%design. It's built from the deprecated DPPF_prep.m code used prior to 2018 PNAS
%paper, modified specifically for the 2019 winter workshop

%user inputs
%watershed: user provides name for the watershed (e.g., Penobscot,
%Presumpscot, Union, etc)
%basedir: base directory where the data are held. I use 'D:\FoD\PPF'. The
%code then look for basedir + '/DamIndex/' and basedir + '/DamData/'
%varargin: list of filters for dam selection, rules
    %v: list of individual criteria quantities for each dam
    %x, or xx, or ipop: initial population provided by user
    %blocklist: for 2019 workshop: must include dams that are included in
    %calculation of criteria quants, but are excluded for decision
    %alternatives (i.e., always kept as-is (1)

%outputs:
%DamIndex: list of dams and immediate downstream neighbor. watersheds with
%islands will have one to a few dams with two downstream neighbors
%v: individual criteria values for each dam in the simulation.
DamDataF = sprintf('%s/DamData/',basedir);
DamIndexF = sprintf('%s/DamIndex/',basedir);

%initialize v, normally I automate sizing, but it's hard-coded for the
%workshop with 16 (!!!) criteria
v=cell(2,16);
D=readtable(sprintf('%sCaseStudyDamsONLY_02182019_SK_SGR.xlsx',DamDataF)); %Hardcoded name of input data

%Dam Indexing: keep only the ones in the data list
DamIndex=dlmread(sprintf('%s%sDownDamList.txt',DamIndexF,watershed),',');
%filter to keep only dams in the watershed of interest
[~,ind1]=ismember(D.DID,DamIndex(:,1));
x=zeros(length(DamIndex),1);
x(ind1)=1; %Dams with x==1 are kept
ActiveDamIndex=DPPF_idx(x,DamIndex,0); %DPPF_idx is a function used to crop dam list
DamIndex=ActiveDamIndex(:,1:2);
[~,ind2]=ismember(DamIndex(:,1),D.DID);
D=D(ind2,:); %crop data list to the dams in the specified region/watershed

%criteria
%fish. Populate the passage values with the multiple alternative values.
v{1,15}='passage';
v{2,15}=[D.PassUp.*D.PassDn ...
    D.PassUpSal.*D.PassDnSal ...
    (1-D.PassUp.*D.PassDn).*0.5+D.PassUp.*D.PassDn ...
    (1-D.PassUpSal.*D.PassDnSal).*0.5+D.PassUpSal.*D.PassDnSal];

%Functional habitat units for the fish, converted to sq km.
v{1,1}='salmon';
v{2,1}=D.SUM_AtlSal./1e6;
v{1,2}='shad';
v{2,2}=D.SUM_ShadMF./1e6;
v{1,3}='alewifeBlueback';
v{2,3}=[D.SUM_AlBkMF./1e6 + D.Al5Dam ...
    D.AlRm5d_skm];

%river recreation
v{1,4}='riverCanoe'; %Ignore river rafting for this workshop. Rafting requires at least 2000cfs, available but unused for this analysis for simplicity.
v{2,4}=[D.FRRU./1e6 ...%12 %functional river recreation units, in sq km
    D.fdr300cfs ...%25 %number of recreation days available given natural flows, with 300cfs threshold (generally good for canoes)
    D.fdy300cfs]; %33 %number of additional recreation days made possible by upstream dam releases, with 300cfs threshold

%reservoir storage in cubic km
v{1,5}='storage';
v{2,5}=D.RsvrVolEst;

%costs: everything costs money in the workshop version. for PNAS, there was
%only one cost array, for removal. This version has 5 for the 5 decision
%alternatives.
v{1,6}='cost';
v{2,6}=D{:,208:212};

%breach damage potential
v{1,7}='hazard';
v{2,7}=D.Hazard;
v{2,7}(v{2,7}>0)=v{2,7}(v{2,7}>0)-1; %alter the hazard scores, 0-2 rather than 1-3. Makes for cleaner post-analysis (safe dams: score of zero)

%number of properties
v{1,8}='property';
v{2,8}=D.StrcRD200; %using only count of properties within 200 m of dam and/or reservoir

%power
v{1,9}='power';
v{2,9}=D{:,215:216}; %there are two power options, not including remove dam (0 kW), with normal and potential improvement
%DamDat=[DamDat D{:,215:216}];

%CO2
v{1,10}='CO2';
v{2,10}=D{:,219}; %CO2 offset thanks to hydropower dams

%cultural 'values'
v{1,11}='indigenous';
v{2,11}=D{:,232:236};
v{1,12}='historic';
v{2,12}=D{:,237:241};
v{1,13}='town';
v{2,13}=D{:,227:231};
v{1,14}='aesthetic';
v{2,14}=D{:,222:226};

D=[]; %clean out the dataset, everything important is moved to v

%ag: constants that are used to standardize and aggregate similar ES
%ag1-3: carrying capacity of river herring, shad, salmon
%ag4: drainage area of the watershed in consideration
%ag5-6: New England lake and river rec densities, respectively
ag=[21.3 24.216 3.002 0 0 0 1 2]; %ag(7) triggers aggregation (1) or not (0)

%ag(8): code indicates the number of decision alternatives to use, based on
%user varargin inputs:
if any(strcmp('ExpPow',varargin))
    ag(8)=ag(8)+1;
end
if any(strcmp('ExpFish',varargin))
    ag(8)=ag(8)+1;
end
if any(strcmp('ExpPow',varargin)) && any(strcmp('ExpFish',varargin))
    ag(8)=ag(8)+1;
end
v{1,16}='ag';
v{2,16}=ag;
