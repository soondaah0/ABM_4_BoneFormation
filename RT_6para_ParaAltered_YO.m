close all;
clear all;
%
% THIS PROGRAM COMPUTES PERIBFR AND ENDOBFR FOR A PARA LIST IN YOUNG AND
% FOR AGED BONE CELL NETWORKS EXPOSED TO DISTINCT PROTOCOLS.  This program
% computes PeriBFR under the null hypothesis that a 12 para ABM (6 paras
% for young and 6 for aged) is sufficient to describe periosteal BFR
% induced by 10 protocols in young animals and by 7 protocols in aged
% animals.  



%% in this program, the fact that sigmoid equations are non zero at below
%% the min threshold has been accounted for by setting the Ca signaling
%% value representive of zero as 0.67 (which is the value of the expression
%% (100*1/(1+exp(5)).  This is true regardless of min and max thresholds
% %% 
% In this program, the osteocyte is the mechanosensor and all cells %
% communicate with all others.  While the model retains the previous
% hypothesis re the real-time signaling module (i.e., 6 paras Temax, Temin,
% Tcmax, Tcmin, Rt, M0), but Tcmax and Tcmin have been assigned arbitrary
% fixed values.  the remaining module has been completely modified.
%
% Now, once real-time signaling is obtained for each protocol in
% osteoblasts, a geometric moving average with an 'infinite' backwards
% window is computed to obtain downstream 'gene' expression that is
% activated due to real-time Ca2+ signaling where 'gene(t) =
% gene(t-1)*Decay+Ca(t)' and where 0<Decay<1.  The gene expression, when above a threshold
% ('GeneThreshold') when integrated over the protocol duration (+extratime)
% provides the osteogenic index ('OsteoIndex') for the protocol.  The
% 'OsteoIndex' is then related to MAR via a two parameter sigmoid that
% defines MARmax and OsteoIndexMax. When a cell has MAR=0, the
% corresponding MS for that cell equals 0, else MS = 100%.  In this
% fashion, MS and MAR are computed for each cell and used to calculate
% tissue level MS, MAR and BFR.
%
%

%%  Initial Parameter and Protocal Input


%  In this cell, the protocals and parameters are input, as are the
%  location of cells.  In addition, the neighborhood is set up for a radius
%  of 5 (need to change once code is finished.  The max time for each
%  protocal is also calculated (might want to do this on an individual
%  basis inside the prot loop instead???

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input the variable sizes here
%   This program will only accept a single parameter list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reading the list of parameter and their values for the ABM
% Assumes the parameters are listed as follows:
% Temax,Temin,Rt,M0,Decay,GeneThreshold,OsteoIndexMax,MARmax

%  Due to prior investigation, TCmin and TCmax parameters were robust (or
%  in Thomas's wording irrelevant) so they are arbitrarily specified as a
%  non-zero value.
Tcmin = 2.5;% TCmax is really irrelevant....

paralist=dlmread('para.txt');
npara=size(paralist,1)
nABM=size(paralist,2) % Number of columns representing number of parameters (no use as of yet)

% input the parameter which is the same for the 8 para ABM between young
% and old ABMs
Samepara=input('Para number in 6 para ABM which is the same between young and aged ABMs?');

if Samepara > 0
    paralist(:,Samepara)=(paralist(:,Samepara)+paralist(:,Samepara+nABM/2))/2;
    paralist(:,Samepara+nABM/2)=paralist(:,Samepara);
end

% reading the protocols
% PROTOCOLS ARE SPECIFIED AS FOLLOWS
% Strain, Total Cycles, Rest Interval, Cycle Grouping (following which, rest is inserted)!
% e.g., 1250,50, 10,5 would implie a protocol inducing 1250 me strain for 50 cycles with a 10-s rest
% inserted between every 5 load cycles.

%  Because this program requires information from both young and aged data
%  sets, we are reading these data sets individual and using them when
%  necessary.

%  Reading old Protocol list
prot_old=dlmread('prot_old7.txt');  %i.e. 'loading waveforms'
nprot_old=size(prot_old,1);  %number of protocals to examine

%  Reading young protocol list
prot_young=dlmread('prot_young.txt');
nprot_young=size(prot_young,1);

%  Combining young and aged protocol lists
prot = [prot_young;prot_old];
nprot = nprot_young + nprot_old;
%nprot=1
% extratime=input('how many extra seconds should real-time signaling be computed for?');
extratime=200;
prottime=((prot(:,4)+prot(:,3)).*prot(:,2)./prot(:,4));
maxtime=((prot(:,4)+prot(:,3)).*prot(:,2)./prot(:,4))+extratime;

%
% program to read cell co-ordinates
% reading osteoblasts; remember that endocortical then periosteal
% osteoblast


%  Reading in cell locations for young animals.
d_young =dlmread('oblast young.txt');
nob_young=size(d_young,1)/2;
obxy_young=reshape(d_young,nob_young,2); % Takes array 'd' and creates a two column matrix with x-cord (column 1) and y-cord (column 2)

e_young=dlmread('ocytes young.txt');
noc_young=size(e_young,1)/2;
ocxy_young=reshape(e_young,noc_young,2); % Same as procedure for obxy above
ncell_young=nob_young+noc_young;

% cellxy order- endo oblasts, peri oblast, ocytes
figure(2);
cellxy_young=[obxy_young; ocxy_young];  % 'cellxy' is a combination of obxy and ocxy with obxy appearing first
plot(cellxy_young(1:nob_young/2,1),cellxy_young(1:nob_young/2,2),'-cd',cellxy_young(nob_young/2+1:nob_young,1),...
    cellxy_young(nob_young/2+1:nob_young,2),'-rd',cellxy_young(nob_young+1:ncell_young,1),cellxy_young(nob_young+1:ncell_young,2),'b.')
hold on  % allows subsequent plotting commands to occur in existing window

% Setting Plot Dimensions
axis([min(cellxy_young(:,1))*1.25,max(cellxy_young(:,1))*1.25,min(cellxy_young(:,2))*1.25,max(cellxy_young(:,2))*1.25])
axis equal %resizes the figure window appropriately


%  Reading in cell loctions for aged animals.
d_old =dlmread('oblast aged.txt');
nob_old=size(d_old,1)/2;
obxy_old=reshape(d_old,nob_old,2); % Takes array 'd' and creates a two column matrix with x-cord (column 1) and y-cord (column 2)

%reading osteocytes for aged
e_old=dlmread('ocytes aged.txt');
noc_old=size(e_old,1)/2;
ocxy_old=reshape(e_old,noc_old,2); % Same as procedure for obxy above
ncell_old=nob_old+noc_old;

% cellxy order- endo oblasts, peri oblast, ocytes
figure(1);
cellxy_old=[obxy_old; ocxy_old];  % 'cellxy' is a combination of obxy and ocxy with obxy appearing first
plot(cellxy_old(1:nob_old/2,1),cellxy_old(1:nob_old/2,2),'-cd',cellxy_old(nob_old/2+1:nob_old,1),...
    cellxy_old(nob_old/2+1:nob_old,2),'-rd',cellxy_old(nob_old+1:ncell_old,1),cellxy_old(nob_old+1:ncell_old,2),'b.')
hold on  % allows subsequent plotting commands to occur in existing window

% Setting Plot Dimensions
axis([min(cellxy_old(:,1))*1.25,max(cellxy_old(:,1))*1.25,min(cellxy_old(:,2))*1.25,max(cellxy_old(:,2))*1.25])
axis equal %resizes the figure window appropriately




%rad=input('what is the canalicular length? ');
rad=5;


%  This section of code calculates the network for young animals.

% computing the network, i.e., what cells are connected to what
% first determining the length between a given cell and all other cells in
% the network using pdist and then saving cell ids if within a specified
% radius
%
distance_young=squareform(pdist(cellxy_young));  %  Measure distance from cell to cell and outputs it in a square matrix
distance_young(1:nob_young,1:nob_young)=0.;
agentset_young = distance_young > 0 & distance_young < rad; % 1 if connected 0 if not
agentset_young = double(agentset_young);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specifying osteoblast to osteoblast connectivity at endo and peri
% surfaces
for i=2:nob_young/2-1
    agentset_young(i,i-1)=1;
    agentset_young(i,i+1)=1;
    agentset_young(i+nob_young/2,i+nob_young/2-1)=1;
    agentset_young(i+nob_young/2,i+nob_young/2+1)=1;
end
% osteoblast connectivities for ob number 1 and nob/2 on endocortical surface
agentset_young(1,2)=1;
agentset_young(1,nob_young/2)=1;
agentset_young(nob_young/2,1)=1;
agentset_young(nob_young/2,nob_young/2-1)=1;
% osteoblast connectivities for ob number nob/2+1 and nob on periosteal surface
agentset_young(1+nob_young/2,nob_young/2+2)=1;
agentset_young(1+nob_young/2,nob_young)=1;
agentset_young(nob_young,1+nob_young/2)=1;
agentset_young(nob_young,nob_young-1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

agentcount_young = sum(agentset_young,2); % Number of connections for each cell
dummyagentcount_young = agentcount_young;
dummyagentcount_young(dummyagentcount_young == 0) = 1;  % This steps allows for calculation of cells w/out neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  This section of code calculates the network for aged animals.

% computing the network, i.e., what cells are connected to what
% first determining the length between a given cell and all other cells in
% the network using pdist and then saving cell ids if within a specified
% radius
%
distance_old=squareform(pdist(cellxy_old));  %  Measure distance from cell to cell and outputs it in a square matrix
distance_old(1:nob_old,1:nob_old)=0.;
agentset_old = distance_old > 0 & distance_old < rad; % 1 if connected 0 if not
agentset_old = double(agentset_old);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specifying osteoblast to osteoblast connectivity at endo and peri
% surfaces
for i=2:nob_old/2-1
    agentset_old(i,i-1)=1;
    agentset_old(i,i+1)=1;
    agentset_old(i+nob_old/2,i+nob_old/2-1)=1;
    agentset_old(i+nob_old/2,i+nob_old/2+1)=1;
end
% osteoblast connectivities for ob number 1 and nob/2 on endocortical surface
agentset_old(1,2)=1;
agentset_old(1,nob_old/2)=1;
agentset_old(nob_old/2,1)=1;
agentset_old(nob_old/2,nob_old/2-1)=1;
% osteoblast connectivities for ob number nob/2+1 and nob on periosteal surface
agentset_old(1+nob_old/2,nob_old/2+2)=1;
agentset_old(1+nob_old/2,nob_old)=1;
agentset_old(nob_old,1+nob_old/2)=1;
agentset_old(nob_old,nob_old-1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

agentcount_old = sum(agentset_old,2); % Number of connections for each cell
dummyagentcount_old = agentcount_old;
dummyagentcount_old(dummyagentcount_old == 0) = 1;  % This steps allows for calculation of cells w/out neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%performing the real-time Ca to cell BFR by calling the following.
%%  Strain Calculation Portion (for all protocals!!!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the Moment of Inertia and Area variables for Strain
% Calculations...

% setting up the Moment of Inertia and Area variables for Strain
% Calculations in young bones
%
Ix_young = 4024620.0d+00;
Iy_young = 3278040.0d+00;
Ixy_young = -202161.0d+00;
Ixyx_young = Ix_young*Iy_young-Ixy_young^2;
Area_young = 4239.0d+00;

%%THESE VALUES ARE FOR AGED BONE- AVE CROSS_SECTION!!!
%
Ix_old = 4.02047d+006;
Iy_old = 3.53954d+006;
Ixy_old = 97445.7d+00;
Ixyx_old = Ix_old*Iy_old-Ixy_old^2;
Area_old = 3733.0d+00;


%determining the strains at the location of each cell based upon peak
%strains induced by a given protocol "prot"
%first determining the force, moment boundary conditions
% THESE VALUES HAVE ALSO BEEN CORRECTED AS APPROPRIATE FOR AGED BONES
%
%  Fz,Mx,My for all protocals (young first/ aged second)
;
Fz(1:nprot_young) = 0.736251/1600.0*prot(1:nprot_young,1);
Fz(nprot_young+1:nprot) = 0.736251/1620.7522*prot(nprot_young+1:nprot,1);
Mx(1:nprot_young) = 23.765/1600.0*prot(1:nprot_young,1);
Mx(nprot_young+1:nprot) = 23.765/1620.7522*prot(nprot_young+1:nprot,1);
My(1:nprot_young) = -112.079/1600.0*prot(1:nprot_young,1);
My(nprot_young+1:nprot) = -112.079/1620.7522*prot(nprot_young+1:nprot,1);

% determining strains and strain induced Ca responses in OSTEOCYTIC
% cells.  Assumes that osteocytes are mechanosensors and osteoblasts
% cannot sense mechanical strain.


%
%%  Input parameters into variables

rpBFR_Para=zeros(npara,nprot);
rpMAR_Para=zeros(npara,nprot);
rpMS_Para=zeros(npara,nprot);

for ipara=1:npara

    para(1:nABM)=paralist(ipara,1:nABM);

    temax=para(1);
    temin=0.0d+00;
    if (temin >= temax) % to prevent singularities
        temin=.9999*temax;
    end



    Rt=para(2);
    CaStore0=para(3);

    alphae=(temax-temin)/10.0d+00;
    alphaStore=CaStore0/10.0d+00;

    %Camin,Camax,MARmax,

    Decay=para(4);% decay factor for gene expression
    GeneThreshold=0.0d+00; %Gene amplitude threshold beyond which gene expression has influence upon MAR

    OsteoIndexMax=para(5);
    MARmax=para(6);

    MARCell=zeros(nprot,max(nob_young,nob_old));


    %%  Protocal Calcium Response (for all protocals!!!)

    %  Only ocytes respond to strain
    %  For each protocal CaStrain=CaStrain(iprot,:)
    CaStrain=zeros(nprot,max(ncell_young,ncell_old));
    Realtime=zeros(nprot,max(maxtime),max(nob_young,nob_old));
    RealtimeGene=zeros(nprot,max(maxtime),max(nob_young,nob_old));
    
    %realenergy=zeros(nprot,max(maxtime),max(ncell_young,ncell_old));
    %%  Protocal Loop
    nob=nob_young;
    ncell=ncell_young;
    noc=ncell-nob;
    cellxy=cellxy_young;
    agentset=agentset_young;
    dummyagentcount=dummyagentcount_young;

    Ix=Ix_young;
    Iy=Iy_young;
    Ixy=Ixy_young;
    Ixyx=Ixyx_young;
    Area=Area_young;
    
    strain = zeros(nprot,max(ncell_young,ncell_old));

    
    for iprot=1:nprot

        if iprot == nprot_young+1
            % setting all bone anatomy and network connections to aged bone
            % values
            nob=nob_old;
            ncell=ncell_old;
            noc=ncell-nob;
            cellxy=cellxy_old;
            agentset=agentset_old;
            dummyagentcount=dummyagentcount_old;

            Ix=Ix_old;
            Iy=Iy_old;
            Ixy=Ixy_old;
            Ixyx=Ixyx_old;
            Area=Area_old;
            
            % Reseting all exlicit para variables to the old ABM values.
            temax=para(1+6);
            temin=0.0d+00;
            Rt=para(2+6);
            CaStore0=para(3+6);
            Decay=para(4+6);% decay factor for gene expression
            GeneThreshold=0.0d+00; %Gene amplitude threshold beyond which gene expression has influence upon MAR
            OsteoIndexMax=para(5+6);
            MARmax=para(6+6);
            
            if (temin >= temax) % to prevent singularities
                temin=.9999*temax;
            end
            alphae=(temax-temin)/10.0d+00;
            alphaStore=CaStore0/10.0d+00;
        end
        
        strain(iprot,1:ncell)=abs(Fz(iprot)./Area+cellxy(1:ncell,2).*(Mx(iprot)*Iy-My(iprot)*Ixy)/Ixyx-...
        cellxy(1:ncell,1).*(My(iprot)*Ix-Mx(iprot)*Ixy)/Ixyx)*1.0d+06;
  
        CaStrain(iprot,nob+1:ncell)=100./(1.+exp(-(strain(iprot,nob+1:ncell)-...
            (temax+temin)/2.0)/alphae));%Ocyte is the mechanosensor

        %  Initiating for storage purposes
        CaCell=zeros(maxtime(iprot),nob);

        Gene=zeros(maxtime(iprot)+extratime/2,nob);

        Energy=ones(1,ncell).*CaStore0;
        CaComm=zeros(1,ncell);  %  Initialize communication equal to zero for t=1
        CaActual=zeros(1,ncell);%  This will be fed into CaCell and used for communication

        %%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%
        % The first index it CaCell is t=1
        % That means there is no t=0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        keystim = 0;     %   If Keystim equals 1 stimulus is OFF
        %   If Keystim equals 0 stimulus is ON

        keycycle = 0;    %   Keycycle used to track cycle groupin


        %   First Time Loop is for Strain portion
        %   Second Time Loop is for Communication during 'Extra Time'
        for t=1:prottime(iprot)

            % replenishing the stores if the cell is quiet
            Energy(CaActual<=1.0d-06)=Energy(CaActual<=1.0d-06)+Rt*1./(1+exp((Energy(CaActual<=1.0d-06)-CaStore0/2)/alphaStore));
            Energy(Energy>CaStore0) = CaStore0;
            %realenergy(iprot,t,1:ncell)=Energy(1:ncell);

            %%%%% keystim equal to 0 means stimulus on
            if keystim == 0

                %new communication try
                CaComm=transpose((agentset*(transpose(CaActual)))./dummyagentcount);
                CaComm(CaComm <= Tcmin) = 0;


                % This sets CaActual to the higher value between strain and
                % communication induced calcium
                CaActual=max([CaStrain(iprot,1:ncell);CaComm],[],1);
                % This sets CaActual to the minimum between the desired
                % response and what is allowable due to energy levels
                CaActual=min([CaActual;Energy],[],1);
            else
                CaComm=transpose((agentset*(transpose(CaActual)))./dummyagentcount);
                CaComm(CaComm <= Tcmin) = 0;
                % This sets CaActual to the minimum between the desired
                % response and what is allowable due to energy levels
                CaActual=min([CaComm;Energy],[],1);
            end

            %%%%%% Setting small signals equal to zero%%%%%
            CaActual(CaActual < 0.67) = 0.0d+00; % correction for the sigmoid value when under min threshold
            % value of 100*1/(1+exp(5))

            % Depleting energy due to response
            Energy = Energy - CaActual;
            Energy(Energy< 0) = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Used to track when stimulus is applied %%%%

            keycycle = keycycle + 1;
            if (keycycle == prot(iprot,4))
                keystim = 1;
            end
            if (keycycle == prot(iprot,3) + prot(iprot,4))
                keystim = 0;
                keycycle = 0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            CaCell(t,:)=CaActual(1:nob);
            if t>1
                Gene(t,:)=Decay*Gene(t-1,:)+CaCell(t,:);%Computing the gene transcription
                %occuring in the real-time based upon geometrically weighted moving
                %averages with an infinite backward window
            end

            Realtime(iprot,t,:)=CaActual(1:nob);
            RealtimeGene(iprot,t,:)=Gene(t,:);

        end


        %%%%% Extra Time %%%%%%
        %Rt=0.01*Rt
        for t=prottime(iprot)+1:maxtime(iprot)
            % replenishing the stores if the cell is quiet

            Energy(CaActual<=1.0d-06)=Energy(CaActual<=1.0d-06)+Rt*1./(1+exp((Energy(CaActual<=1.0d-06)-CaStore0/2)/alphaStore));
            Energy(Energy>CaStore0)=CaStore0;
            %realenergy(iprot,t,1:ncell)=Energy(1:ncell);

            CaComm=transpose((agentset*(transpose(CaActual)))./dummyagentcount);
            CaComm(CaComm <= Tcmin) = 0;

            % This sets CaActual to the minimum between the desired
            % response and what is allowable due to energy levels
            CaActual=min([CaComm;Energy],[],1);

            %%%%%% Setting small signals equal to zero%%%%%
            CaActual(CaActual < 0.67) = 0.0d+00;% correction for the sigmoid value when under min threshold

            % Depleting energy due to response
            Energy = Energy - CaActual;
            Energy(Energy< 0) = 0;

            CaCell(t,:)=CaActual(1:nob);
            Gene(t,:)=Decay*Gene(t-1,:)+CaCell(t,:);%Computing the gene transcription
            %occuring in the real-time based upon geometrically weighted moving
            %averages with an infinite backward window

            Realtime(iprot,t,:)=CaActual(1:nob);
            RealtimeGene(iprot,t,:)=Gene(t,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Computing Gene transcription for a bit after the protocol ends
        %for all but the largest protocol (typically rest-inserted,
        %which really does not require any more computation than that
        %allowed by the 'extratime' number...
        endtime=max(maxtime);
        if maxtime(iprot)<max(maxtime)
            endtime=maxtime(iprot)+extratime/2;
            CaCell=[CaCell;zeros(extratime/2,nob)];
        end

        for t=maxtime(iprot)+1:endtime
            Gene(t,:)=Decay*Gene(t-1,:)+CaCell(t,:);%Computing the gene transcription
            %occuring in the real-time based upon geometrically weighted moving averages with an infinite backward window
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Gene=Gene-GeneThreshold;
        Gene(Gene<0)=0.;
        OsteoIndex=zeros(1,nob);
        OsteoIndex=sum(Gene,1);
        MARCell(iprot,OsteoIndex>0)=MARmax*1./(1+exp(-(OsteoIndex(OsteoIndex>0)-OsteoIndexMax/2)/0.1/OsteoIndexMax));
        %
        %
    end
    % %


    BFRPeri=zeros(nprot,1);%periosteal tissue level BFR
    BFREndo=zeros(nprot,1);%endocortical tissue level BFR
    MARPeri=zeros(nprot,1);%periosteal tissue level MAR
    MAREndo=zeros(nprot,1);%endocortical tissue level MAR
    MSPeri=zeros(nprot,1);%periosteal tissue level MS
    MSEndo=zeros(nprot,1);%endocortical tissue level MS

    for iprot=1:nprot
        if nnz(MARCell(iprot,(nob/2 + 1):nob)) ~= 0
            MSPeri(iprot)=nnz(MARCell(iprot,(nob/2 + 1):nob))/nob*2*100;% MS is simply the no of cells that have nonzero MAR
            % Previously (before 04/08/2010), the MARPeri was computed
            % in-correctly where the MSPeri...was not enclosed within
            % brackets; this BUG has now been corrected
            MARPeri(iprot)=sum(MARCell(iprot,1+nob/2:nob))/(MSPeri(iprot)*nob/2/100.); %MAR is the non-zero ave of MAR in each cell
            BFRPeri(iprot)=MSPeri(iprot)*MARPeri(iprot)/100.;
        end
        %
        if nnz(MARCell(iprot,1:nob/2)) ~= 0
            MSEndo(iprot)=nnz(MARCell(iprot,1:nob/2))/nob*2*100;% MS is simply the no of cells that have nonzero MAR
            MAREndo(iprot)=sum(MARCell(iprot,1:nob/2))/(MSEndo(iprot)*nob/2/100.); %MAR is the non-zero ave of MAR in each cell
            BFREndo(iprot)=MSEndo(iprot)*MAREndo(iprot)/100.;
        end
    end
    %
    rpBFR_Para(ipara,:)=BFRPeri(:);
    rpMAR_Para(ipara,:)=MARPeri(:);
    rpMS_Para(ipara,:)=MSPeri(:);
end

toc

hold off
