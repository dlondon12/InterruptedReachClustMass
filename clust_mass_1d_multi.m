function [h,pout,clusts,Fstat,RRrandout,Fclustrandmax,hposthoc,poutposthoc,clustspair,pairs]=clust_mass_1d_multi(RR,p,pclust,Niter,varargin)

% This function performs the within-units cluster mass test to test whether
% multiple conditions (Ngroups) have distinct firing rates during adjacent
% points in time. It was used in this paper:
% https://www.biorxiv.org/content/10.1101/2020.11.11.378224v1. A companion
% mansucript will describe statistical properties, motivation, and
% generalization of this test.
% 
% Inputs:
% 
% RR: Ntimes x Ntrials x Ngroups array
%   Ntrials can refer trials if all conditions occur within the same trial
%   or neurons if all conditions are tested in each neuron. In this case each
%   Ntimes x 1 x Ngroups sub-matrix refers to the PSTH of each neuron for
%   each condition.
%
% p: scalar
%   p-value that defines the threshold to define clusters
%
% pclust: scalar
%   p-value that defines significance for each cluster
%
% Niter: scalar
%   Number of permutations
%
% posthoc: scalar: 1, optional parameter
%   Flag that determines whether to run follow up pairwise tests on clusters
%   that are statistically significant. Set to 1 to run post-hoc tests.
%   This feature was not used in the paper above.
%
%
% Outputs:
% h: Ntimes x 1 logical vector
%   True for time points in significant clusters, false otherwise
%
% pout: Ntimes x 1 vector
%   p-value for each point in clusters that are tested. These are the
%   p-values of each cluster, not separate p-values for each time point.
%   Time points not in a cluster are set to NaN.
%
% clusts: Nclusts x 1 struct
%   Structure containing p-values and time-points of each cluster
%
% Fstat: Ntimes x 1 vector
%   F-statistics used to generate the clusters
%
% RRrandout: Ntimes x Ntrials x Ngroups x maxiter array
%   Permuted data for the last iteration of permutations (see code below).
%   To output all permutations set maxiter equal to Niter.
%
% Fclustrandmax: Niter x 1 vector
%   Null distribution of cluster statistics
% 
% hposthoc: [Ngroups choose 2] x Ntimes logical matrix
%   For each pair of conditions, this is true if that timepoint is in a
%   significant cluster identified by post-hoc testing. Otherwise false.
%
% poutposthoc: [Ngroups choose 2] x Ntimes matrix
%   Same as pout but for post-hoc tests
%
% clustspair: [Ngroups choose 2] x 1 struct
%   Contains strutures with information about each post-hoc cluster
%   analogous to the clusts structure
%
% pairs: [Ngroups choose 2] x 2 matrix
%   Rows indicate which pairs of conditions were compared for each post-hoc
%   test. Each row corresponds to each row of the other post-hoc outputs.


if nargin==5
    posthoc=varargin{1};
else
    posthoc=0;
end

Ntimes=size(RR,1);
Ntrials=size(RR,2);
ntrialstime=mean(sum(~isnan(RR),2),3); % number of trials/units at each point in time, allows missing data (NaN) at certain points in time
ngroupstime=mean(sum(~isnan(RR),3),2); % number of groups/conditions at each point in time
Ngroups=size(RR,3);

grandmean=mean(reshape(RR,Ntimes,Ntrials*Ngroups),2,'omitnan');
groupmean=mean(RR,2,'omitnan');
unitmean=mean(RR,3,'omitnan');

% ssbetween=Ntrials.*sum(bsxfun(@minus,squeeze(groupmean),grandmean).^2,2);
% sswithin=sum(sum(bsxfun(@minus,RR,groupmean).^2,3,'omitnan'),2,'omitnan');
% ssunit=Ngroups.*sum(bsxfun(@minus,unitmean,grandmean).^2,2);

% Support for NaN values in this code which adjusts sums of squares for
% number of non-missing data at each time point.
ssbetween=ntrialstime.*sum(bsxfun(@minus,squeeze(groupmean),grandmean).^2,2);
sswithin=sum(sum(bsxfun(@minus,RR,groupmean).^2,3,'omitnan'),2,'omitnan');
ssunit=ngroupstime.*sum(bsxfun(@minus,unitmean,grandmean).^2,2);

sserror=sswithin-ssunit;
Fstat=(ssbetween./(Ngroups-1))./(sserror./((Ntrials-1)*(Ngroups-1)));
% Fstat=(ssbetween./(ngroupstime-1))./(sserror./((ntrialstime-1).*(ngroupstime-1)));

threshold=finv(1-p,Ngroups-1,(Ntrials-1)*(Ngroups-1)); % use provided p-value to set the threshold for determining clusters of points to test

sig_ind=Fstat>threshold; %points which define clusters

% Next we permute the data and repeat the above. Doing this for a large
% number of permutations at a time in vectorized fashion improves
% computational speed until memory limitations play a role.

maxIter=100; % maximum number of permutations to perform at a time
iterblocks=0:maxIter:Niter; % number of permuations in each iteration 
iterblocks(end)=Niter;

for iter=1:(size(iterblocks,2)-1)
    niter=iterblocks(iter+1)-iterblocks(iter); % number of permutations in this iteration
    RRrand=zeros(Ntimes,Ntrials,Ngroups,niter);
    perms=zeros(niter,Ntrials,Ngroups);
    
    for i=1:niter
        for j=1:Ntrials
            perms(i,j,:)=randperm(Ngroups); % permute the group labels separately within each trial/unit
        end
    end
    
    for m=1:Ntrials
        for k=1:Ngroups
            for i=1:niter
                RRrand(:,m,k,i)=RR(:,m,perms(i,m,k)); % form the permuted data
            end
        end
    end
    RRrandout=RRrand;
    %grandmeanrand is equal to grandmean
    groupmeanrand=mean(RRrand,2,'omitnan');
    unitmeanrand=mean(RRrand,3,'omitnan');
    
    ssbetweenrand=squeeze(Ntrials.*sum(bsxfun(@minus,squeeze(groupmeanrand),grandmean).^2,2));
    sswithinrand=squeeze(sum(sum(bsxfun(@minus,RRrand,groupmeanrand).^2,3,'omitnan'),2,'omitnan'));
    ssunitrand=squeeze(Ngroups.*sum(bsxfun(@minus,unitmeanrand,grandmean).^2,2));
    
    sserrorrand=sswithinrand-ssunitrand;
    Fstatrand=(ssbetweenrand./(Ngroups-1))./(sserrorrand./((Ntrials-1)*(Ngroups-1)));
    
    sig_ind_rand=Fstatrand>threshold;
    
    
    clusts=regionprops(sig_ind,'PixelList'); % find adjacent point and form clusters
    Fclust=zeros(1,length(clusts));
    
    Fstat=sserror; % we will use the sum of the error sum of squares as the test statistic for each cluster
    Fstatrand=sserrorrand;
    
    for ii=1:length(clusts)
        Fclust(ii)=sum(arrayfun(@(x) Fstat(x),clusts(ii).PixelList(:,2))); % this is the test statistic for each cluster
    end
    
    
    % Next we calculate the maximum test statistic on each permutation. If 
    % there are no clusters on a permutation, we set the maximum test
    % statistic to zero.
    
    for i=1:niter
        
        clusts_rand=regionprops(sig_ind_rand(:,i),'PixelList');
        
        if ~isempty(clusts_rand)
            Fclust_rand=zeros(1,length(clusts_rand));
            for ii=1:length(clusts_rand)
                Fclust_rand(ii)=sum(arrayfun(@(x) Fstatrand(x,i),clusts_rand(ii).PixelList(:,2)));
            end
            Fclustrandmax(iterblocks(iter)+i)=max(Fclust_rand);
        else
            Fclustrandmax(iterblocks(iter)+i)=0;
        end
        
    end
end

pclusts=1-sum(bsxfun(@gt,Fclust',Fclustrandmax),2)/Niter; % p-value of each cluster
pout=NaN.*ones(1,Ntimes);

coords=[];
for ii=1:length(pclusts)
    clusts(ii).p=pclusts(ii);
    pout(clusts(ii).PixelList(:,2))=pclusts(ii);
    
    if clusts(ii).p<pclust
        coords=[coords;clusts(ii).PixelList];
    end
end

if isempty(coords)
    h=false(1,Ntimes);
else
    h=sparse(coords(:,1),coords(:,2),1,1,Ntimes);
    
    h=logical(full(h));
end


% post-hoc tests (analogous to above)

if posthoc==1
    pairs=nchoosek(1:Ngroups,2);
    Fclustrandmat=zeros(length(pairs),Niter);
    clustspair=cell(1,length(pairs));
    Fclustpair=clustspair;
    Ngroups=2;
    Ntimes=nnz(h);
    for r=1:length(pairs)
        RRpair=RR(:,:,pairs(r,:));
        RRpair=RRpair(h,:,:);
        
        grandmean=mean(reshape(RRpair,Ntimes,Ntrials*Ngroups),2,'omitnan');
        groupmean=mean(RRpair,2,'omitnan');
        unitmean=mean(RRpair,3,'omitnan');
        
        ssbetween=Ntrials.*sum(bsxfun(@minus,squeeze(groupmean),grandmean).^2,2);
        sswithin=sum(sum(bsxfun(@minus,RRpair,groupmean).^2,3,'omitnan'),2,'omitnan');
        ssunit=Ngroups.*sum(bsxfun(@minus,unitmean,grandmean).^2,2);
        
        sserror=sswithin-ssunit;
        Fstat=(ssbetween./(Ngroups-1))./(sserror./((Ntrials-1)*(Ngroups-1)));
        
        threshold=finv(1-p,Ngroups-1,(Ntrials-1)*(Ngroups-1));
        
        sig_ind=Fstat>threshold;
        
        perms=zeros(Niter,Ntrials,Ngroups);
        for i=1:Niter
            for j=1:Ntrials
                perms(i,j,:)=randperm(Ngroups);
            end
        end
        
        RRrand=zeros(Ntimes,Ntrials,Ngroups,Niter);
        
        for m=1:Ntrials
            for k=1:Ngroups
                for i=1:Niter
                    RRrand(:,m,k,i)=RRpair(:,m,perms(i,m,k));
                end
            end
        end
        
        %grandmeanrand is equal to grandmean
        
        groupmeanrand=mean(RRrand,2,'omitnan');
        unitmeanrand=mean(RRrand,3,'omitnan');
        
%         ssbetweenrand=squeeze(Ntrials.*sum(bsxfun(@minus,squeeze(groupmeanrand),grandmean).^2,2));
%         sswithinrand=squeeze(sum(sum(bsxfun(@minus,RRrand,groupmeanrand).^2,3,'omitnan'),2,'omitnan'));
%         ssunitrand=squeeze(Ngroups.*sum(bsxfun(@minus,unitmeanrand,grandmean).^2,2));
        
        ssbetweenrand=squeeze(ntrialstime.*sum(bsxfun(@minus,squeeze(groupmeanrand),grandmean).^2,2));
        sswithinrand=squeeze(sum(sum(bsxfun(@minus,RRrand,groupmeanrand).^2,3,'omitnan'),2,'omitnan'));
        ssunitrand=squeeze(ngroupstime.*sum(bsxfun(@minus,unitmeanrand,grandmean).^2,2));
        
        sserrorrand=sswithinrand-ssunitrand;
%         Fstatrand=(ssbetweenrand./(Ngroups-1))./(sserrorrand./((Ntrials-1)*(Ngroups-1)));
        
        Fstatrand=(ssbetweenrand./(ngroupstime-1))./(sserrorrand./((ntrialstime-1)*(ngroupstime-1)));
        
        sig_ind_rand=Fstatrand>threshold;
        
        
        clustspair{r}=regionprops(sig_ind,'PixelList');
        Fclustpair{r}=zeros(1,length(clustspair{r}));
        
        for ii=1:length(clustspair{r})
            Fclustpair{r}(ii)=sum(arrayfun(@(x) sserror(x),clustspair{r}(ii).PixelList(:,2)));
        end
        
        for i=1:Niter
            
            clusts_rand=regionprops(sig_ind_rand(:,i),'PixelList');
           
            
            if ~isempty(clusts_rand)
                Fclust_rand=zeros(1,length(clusts_rand));
                for ii=1:length(clusts_rand)
                    Fclust_rand(ii)=sum(arrayfun(@(x) sserrorrand(x,i),clusts_rand(ii).PixelList(:,2)));
                end
                Fclustrandmat(r,i)=max(Fclust_rand);
            else
                Fclustrandmat(r,i)=0;
            end
            
        end
        
    end
    
    Fclustrandmat=max(Fclustrandmat);
    poutpair=NaN.*ones(length(pairs),Ntimes);
    hpair=false(length(pairs),Ntimes);
    for r=1:length(pairs)
        pclusts=1-sum(bsxfun(@gt,Fclustpair{r}',Fclustrandmat),2)/Niter;
        
        
        coords=[];
        for ii=1:length(pclusts)
            clustspair{r}(ii).p=pclusts(ii);
            poutpair(r,clustspair{r}(ii).PixelList(:,2))=pclusts(ii);
            
            if clustspair{r}(ii).p<pclust
                coords=[coords;clustspair{r}(ii).PixelList];
            end
        end
        
        if ~isempty(coords)
            htemp=sparse(coords(:,1),coords(:,2),1,1,Ntimes);
            
            hpair(r,:)=logical(full(htemp));
        end
    end
    hposthoc=false(length(pairs),size(RR,1));
    poutposthoc=NaN.*ones(length(pairs),size(RR,1));
    hposthoc(repmat(h,3,1))=hpair;
    poutposthoc(repmat(h,3,1))=poutpair;
end
end