function [hmat,pout,clusts]=clust_mass_1d_2f_multi(RR,p,pclust,Niter)

% This function performs the within-units cluster mass test to test whether
% multiple conditions across two factors (NgroupsA and Ngroups B) have
% distinct firing rates during adjacent points in time. It was used in this
% paper: https://www.biorxiv.org/content/10.1101/2020.11.11.378224v1. A
% companion mansucript will describe statistical properties, motivation,
% and generalization of this test. This function should be considered to
% assumes there is no interaction of the two factors. Dealing with the
% interaction effect requires a different permutation strategy. This will
% be addressed in an updated function in the companion manuscript.
% 
% Inputs:
% 
% RR: Ntimes x Ntrials x NgroupsA x NgroupsB array
%   Ntrials can refer trials if all conditions occur within the same trial
%   or neurons if all conditions are tested in each neuron. In this case each
%   Ntimes x 1 x NgroupsA x NgroupsB sub-array refers to the PSTH of each
%   neuron for each combination of conditions.
%
% p: scalar
%   p-value that defines the threshold to define clusters
%
% pcust: scalar
%   p-value that defines significance for each cluster
%
% Niter: scalar
%   Number of permutations
%
%
% Outputs:
% h: Ntimes x 1 logical vector
%   True for time points in significant clusters, false otherwise
%
% pout: Ntimes x 3 vector
%   p-value for each point in clusters that are tested. These are the
%   p-values of each cluster, not separate p-values for each time point.
%   Time points not in a cluster are set to NaN. The first row contains the
%   p-values for the effect corresponding to NgroupsA and the second row
%   for the effect corresponding to NgroupsB. The third row provides
%   interaction effect p-values, but these are overly conservative. In
%   general permutation across both factors should not be done if there is
%   an interaction and this requires an approximate test. This will be
%   addressed in a companion paper with an updated function.
%
% clusts: Nclusts x 3 struct
%   Structure containing p-values and time-points of each cluster for each
%   main effect and the interaction effect.

Ntimes=size(RR,1);
Ntrials=size(RR,2);
NgroupsA=size(RR,3);
NgroupsB=size(RR,4);
Ntrialsvec=sum(sum(reshape(~isnan(RR),Ntimes,Ntrials,NgroupsA*NgroupsB),3)==NgroupsA*NgroupsB,2);
%number of units/trials at each time point with complete data for both factors

grandmean=mean(reshape(RR,Ntimes,Ntrials*NgroupsA*NgroupsB),2,'omitnan');
groupAmean=permute(mean(reshape(permute(RR,[1,3,2,4]),Ntimes,NgroupsA,NgroupsB*Ntrials),3,'omitnan'),[1,3,2]);
groupBmean=permute(mean(reshape(permute(RR,[1,4,2,3]),Ntimes,NgroupsB,NgroupsA*Ntrials),3,'omitnan'),[1,3,4,2]);
trialsmean=mean(reshape(RR,Ntimes,Ntrials,NgroupsA*NgroupsB),3,'omitnan');

groupABmean=mean(RR,2,'omitnan');
groupAtrialsmean=mean(RR,4,'omitnan');
groupBtrialsmean=mean(RR,3,'omitnan');

ssA=Ntrialsvec.*NgroupsB.*sum((bsxfun(@minus,groupAmean,grandmean)).^2,3,'omitnan');
ssB=Ntrialsvec.*NgroupsA.*sum((bsxfun(@minus,groupBmean,grandmean)).^2,4,'omitnan');
sstrials=NgroupsA*NgroupsB.*sum(bsxfun(@minus,trialsmean,grandmean).^2,2,'omitnan');

ssAB=Ntrialsvec.*sum(reshape((bsxfun(@minus,bsxfun(@minus,groupABmean,groupAmean),...
    groupBmean)+grandmean).^2,Ntimes,NgroupsA*NgroupsB),2,'omitnan');

ssAtrials=NgroupsB*sum(reshape((bsxfun(@minus,bsxfun(@minus,groupAtrialsmean,groupAmean),...
    trialsmean)+grandmean).^2,Ntimes,NgroupsA*Ntrials),2,'omitnan');

ssBtrials=NgroupsA*sum(reshape((bsxfun(@minus,bsxfun(@minus,groupBtrialsmean,groupBmean),...
    trialsmean)+grandmean).^2,Ntimes,NgroupsB*Ntrials),2,'omitnan');

sstot=sum(reshape((RR-grandmean).^2,Ntimes,Ntrials*NgroupsA*NgroupsB),2,'omitnan');

ssE=sstot-ssA-ssB-sstrials-ssAB-ssAtrials-ssBtrials;


dfA=NgroupsA-1;
dfB=NgroupsB-1;
dfAB=(NgroupsA-1)*(NgroupsB-1);
dfAtrials=(NgroupsA-1)*(Ntrialsvec-1);
dfBtrials=(NgroupsB-1)*(Ntrialsvec-1);
dfE=(NgroupsA-1)*(NgroupsB-1)*(Ntrialsvec-1);

FA=(ssA./dfA)./(ssAtrials./dfAtrials);
FB=(ssB./dfB)./(ssBtrials./dfBtrials);
FAB=(ssAB./dfAB)./(ssE./dfE);
%Fstat=[FA';FB';FAB'];

thresholdA=finv(1-p,dfA,dfAtrials);
thresholdB=finv(1-p,dfB,dfBtrials);
thresholdAB=finv(1-p,dfAB,dfE);


sig_indA=FA>thresholdA;
sig_indB=FB>thresholdB;
sig_indAB=FAB>thresholdAB;
sig_ind=[sig_indA';sig_indB';sig_indAB'];

Fstat=[ssAtrials';ssBtrials';ssE']; %we use the error sums of squares for our cluster statistics
for j=1:3
    
    clusts{j}=regionprops(sig_ind(j,:),'PixelList');
    Fclust{j}=zeros(1,length(clusts{j}));
    clustscell=cell(length(clusts{j}),1);
    
    if size(clusts{j},1)>0
        for ii=1:length(clusts{j})
            Fclust{j}(ii)=sum(arrayfun(@(x) Fstat(j,x),clusts{j}(ii).PixelList(:,1)));
            clustscell{ii}=clusts{j}(ii).PixelList(:,1);
        end
    end
end

maxIter=10;
iterblocks=0:maxIter:Niter;
iterblocks(end)=Niter;

Fclustrandmax=NaN*ones(3,Niter);
pout=NaN.*ones(3,Ntimes);
hmat=NaN*ones(3,Ntimes);

for iter=1:(size(iterblocks,2)-1)
    niter=iterblocks(iter+1)-iterblocks(iter);
    permsA=zeros(niter,Ntrials,NgroupsA);
    for i=1:niter
        for j=1:Ntrials
            permsA(i,j,:)=randperm(NgroupsA);
        end
    end
    
    RRrandA=zeros(Ntimes,Ntrials,NgroupsA,NgroupsB,niter);
    
    for m=1:Ntrials
        for k=1:NgroupsA
            for i=1:niter
                RRrandA(:,m,k,:,i)=RR(:,m,permsA(i,m,k),:); %permute labels of condition A without permuting those on condition B
            end
        end
    end
    
    permsB=zeros(niter,Ntrials,NgroupsB);
    for i=1:niter
        for j=1:Ntrials
            permsB(i,j,:)=randperm(NgroupsB);
        end
    end
    
    RRrandB=zeros(Ntimes,Ntrials,NgroupsA,NgroupsB,niter);
    
    for m=1:Ntrials
        for k=1:NgroupsB
            for i=1:niter
                RRrandB(:,m,:,k,i)=RR(:,m,:,permsB(i,m,k));
            end
        end
    end
    
    permsAB=zeros(niter,Ntrials,NgroupsA*NgroupsB);
    for i=1:niter
        for j=1:Ntrials
            permsAB(i,j,:)=randperm(NgroupsA*NgroupsB); %permute labels of conditions A and B together, this results in an overly conservative test
        end
    end
    
    RRrandAB=zeros(Ntimes,Ntrials,NgroupsA,NgroupsB,niter);
    
    for m=1:Ntrials
        for k=1:(NgroupsA*NgroupsB)
            for i=1:niter
                RRrandAB(:,m,mod(k+NgroupsA-1,NgroupsA)+1,ceil(k/NgroupsA),i)=RR(:,m,mod(permsAB(i,j,k)+NgroupsA-1,NgroupsA)+1,ceil(permsAB(i,m,k)/NgroupsA));
            end
        end
    end
    
    %grandmeanrand is equal to grandmean
    
    groupAmeanrandA=permute(mean(reshape(permute(RRrandA,[1,5,3,2,4]),Ntimes,niter,NgroupsA,NgroupsB*Ntrials),4,'omitnan'),[1,4,3,5,2]);
    trialsmeanrandA=permute(mean(reshape(permute(RRrandA,[1,5,2,3,4]),Ntimes,niter,Ntrials,NgroupsA*NgroupsB),4,'omitnan'),[1,3,4,5,2]);
    groupAtrialsmeanrandA=mean(RRrandA,4,'omitnan');
    
    ssArandA=Ntrialsvec.*NgroupsB.*sum((bsxfun(@minus,groupAmeanrandA,grandmean)).^2,3);
    ssAtrialsrandA=NgroupsB*sum(sum((bsxfun(@minus,bsxfun(@minus,groupAtrialsmeanrandA,groupAmeanrandA),...
        trialsmeanrandA)+grandmean).^2,2,'omitnan'),3,'omitnan');
    
    FArandA=squeeze((ssArandA./dfA)./(ssAtrialsrandA./dfAtrials));
    
    groupBmeanrandB=permute(mean(reshape(permute(RRrandB,[1,5,4,2,3]),Ntimes,niter,NgroupsB,NgroupsA*Ntrials),4,'omitnan'),[1,4,5,3,2]);
    trialsmeanrandB=permute(mean(reshape(permute(RRrandB,[1,5,2,3,4]),Ntimes,niter,Ntrials,NgroupsA*NgroupsB),4,'omitnan'),[1,3,4,5,2]);
    groupBtrialsmeanrandB=mean(RRrandB,3,'omitnan');
    
    ssBrandB=Ntrialsvec.*NgroupsA.*sum((bsxfun(@minus,groupBmeanrandB,grandmean)).^2,4);
    ssBtrialsrandB=NgroupsA*sum(sum((bsxfun(@minus,bsxfun(@minus,groupBtrialsmeanrandB,groupBmeanrandB),...
        trialsmeanrandB)+grandmean).^2,2,'omitnan'),4,'omitnan');
    
    FBrandB=squeeze((ssBrandB./dfB)./(ssBtrialsrandB./dfBtrials));
    
    groupAmeanrandAB=permute(mean(reshape(permute(RRrandAB,[1,5,3,2,4]),Ntimes,niter,NgroupsA,NgroupsB*Ntrials),4,'omitnan'),[1,4,3,5,2]);
    groupBmeanrandAB=permute(mean(reshape(permute(RRrandAB,[1,5,4,2,3]),Ntimes,niter,NgroupsB,NgroupsA*Ntrials),4,'omitnan'),[1,4,5,3,2]);
    trialsmeanrandAB=permute(mean(reshape(permute(RRrandAB,[1,5,2,3,4]),Ntimes,niter,Ntrials,NgroupsA*NgroupsB),4,'omitnan'),[1,3,4,5,2]);
    groupABmeanrandAB=mean(RRrandAB,2,'omitnan');
    groupAtrialsmeanrandAB=mean(RRrandAB,4,'omitnan');
    groupBtrialsmeanrandAB=mean(RRrandAB,3,'omitnan');
    
    ssArandAB=Ntrialsvec.*NgroupsB.*sum((bsxfun(@minus,groupAmeanrandAB,grandmean)).^2,3);
    ssBrandAB=Ntrialsvec.*NgroupsA.*sum((bsxfun(@minus,groupBmeanrandAB,grandmean)).^2,4);
    sstrialsrandAB=NgroupsA*NgroupsB*sum((bsxfun(@minus,trialsmeanrandAB,grandmean)).^2,2,'omitnan');
    
    ssAtrialsrandAB=NgroupsB*sum(sum((bsxfun(@minus,bsxfun(@minus,groupAtrialsmeanrandAB,groupAmeanrandAB),...
        trialsmeanrandAB)+grandmean).^2,2,'omitnan'),3,'omitnan');
    ssBtrialsrandAB=NgroupsA*sum(sum((bsxfun(@minus,bsxfun(@minus,groupBtrialsmeanrandAB,groupBmeanrandAB),...
        trialsmeanrandAB)+grandmean).^2,2,'omitnan'),4,'omitnan');
    ssABrandAB=Ntrialsvec.*sum(sum((bsxfun(@minus,bsxfun(@minus,groupABmeanrandAB,groupAmeanrandAB),...
        groupBmeanrandAB)+grandmean).^2,3,'omitnan'),4,'omitnan');
    
    sstotrandAB=permute(sum(reshape(permute(bsxfun(@minus,RRrandAB,grandmean).^2,[1,5,2,3,4]),Ntimes,niter,NgroupsA*NgroupsB*Ntrials),3,'omitnan'),[1,3,4,5,2]);
    
    ssErandAB=sstotrandAB-ssArandAB-ssBrandAB-sstrialsrandAB-ssAtrialsrandAB-ssBtrialsrandAB-ssABrandAB;
    
    FABrandAB=squeeze((ssABrandAB./dfAB)./(ssErandAB./dfE));
    
    sig_ind_randA=FArandA>thresholdA;
    sig_ind_randB=FBrandB>thresholdB;
    sig_ind_randAB=FABrandAB>thresholdAB;
   
    sig_ind_rand=permute(cat(3,sig_ind_randA,sig_ind_randB,sig_ind_randAB),[3,1,2]);
    
    Fstatrand=permute(cat(2,ssAtrialsrandA,ssBtrialsrandB,ssErandAB),[2,1,5,3,4]);
    
    %Fstatrand=permute(cat(3,FArandA,FBrandB,FABrandAB),[3,1,2]);
    
    
    
    for j=1:3
        
        if size(clusts{j},1)>0
            
            for i=1:niter
                
                clusts_rand=regionprops(sig_ind_rand(j,:,i),'PixelList');
                
                if ~isempty(clusts_rand)
                    Fclust_rand=zeros(1,length(clusts_rand));
                    for ii=1:length(clusts_rand)
                        Fclust_rand(ii)=sum(arrayfun(@(x) Fstatrand(j,x,i),clusts_rand(ii).PixelList(:,1)));
                    end
                    Fclustrandmax(j,iterblocks(iter)+i)=max(Fclust_rand);
                else
                    Fclustrandmax(j,iterblocks(iter)+i)=0;
                end
                
            end
        end
    end
end

for j=1:3
    if size(clusts{j},1)>0
        pclusts=1-sum(bsxfun(@gt,Fclust{j}',Fclustrandmax(j,:)),2)/Niter;
    else
        pclusts=[];
    end

coords=[];
for ii=1:length(pclusts)
    clusts{j}(ii).p=pclusts(ii);
    pout(j,clusts{j}(ii).PixelList(:,1))=pclusts(ii);
    
    if clusts{j}(ii).p<pclust
        coords=[coords;clusts{j}(ii).PixelList];
    end
end

if isempty(coords)
    h=false(1,Ntimes);
else
    h=sparse(coords(:,1),coords(:,2),1,Ntimes,1);
    
    h=logical(full(h));
end
hmat(j,:)=h;
end
end