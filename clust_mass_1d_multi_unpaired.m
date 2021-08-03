function [h,pout,clusts]=clust_mass_1d_multi_unpaired(RR,p,pclust,Niter)
% This function performs the across trials cluster mass test to test whether
% multiple conditions (Ngroups) have distinct firing rates during adjacent
% points in time. It was used in this paper:
% https://www.biorxiv.org/content/10.1101/2020.11.11.378224v1. A companion
% mansucript will describe statistical properties, motivation, and
% generalization of this test.
% 
% Inputs:
% 
% RR: Ntimes x Ntrials x Ngroups array
%   Firing rate on each trial on each condition. If the number of trials
%   are unequal in the different conditoins the Ntimes x 1 x 1 sub-vector
%   should contain all NaNs. These trials are analyzed in an unpaired
%   manner (i.e., each trial contains only one condition).
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

Ntimes=size(RR,1);
Ntrials=size(RR,2);
Ngroups=size(RR,3);

ntrials=squeeze(sum(~isnan(RR),2));
ntrialstype=max(ntrials);
tottrials=sum(ntrialstype);

grandmean=mean(reshape(RR,Ntimes,Ntrials*Ngroups),2,'omitnan');
groupmean=mean(RR,2,'omitnan');

ssbetween=sum(ntrials.*bsxfun(@minus,squeeze(groupmean),grandmean).^2,2,'omitnan');
sstot=sum(sum(bsxfun(@minus,RR,grandmean).^2,2,'omitnan'),3,'omitnan');
sserror=sstot-ssbetween;
Fstat=(ssbetween./(Ngroups-1))./(sserror./(sum(ntrials,2)-Ngroups));

threshold=finv(1-p,(Ngroups-1)*ones(Ntimes,1),sum(ntrials,2)-Ngroups);

sig_ind=Fstat>threshold;

count=0;
subind=zeros(tottrials,3);
for k=1:Ngroups
    subind((count+1):(count+ntrialstype(k)),3)=k;
    subind((count+1):(count+ntrialstype(k)),2)=1:ntrialstype(k);
    count=count+ntrialstype(k);
end


linind=repmat(subind,Ntimes,1);
linind(:,1)=reshape(repmat(1:Ntimes,tottrials,1),tottrials*Ntimes,1);
linind=sub2ind([Ntimes,Ntrials,Ngroups],linind(:,1),linind(:,2),linind(:,3));

maxIter=1000;
iterblocks=0:maxIter:Niter;
iterblocks(end)=Niter;

for iter=1:(size(iterblocks,2)-1)
    niter=iterblocks(iter+1)-iterblocks(iter);
    RRrand=NaN*ones(Ntimes,Ntrials,Ngroups,iter);
    %perms=zeros(niter,Ntrials,Ngroups);
    perms=zeros(niter,tottrials);
    
    for i=1:niter
        %for j=1:Ntrials
        %    perms(i,j,:)=randperm(Ngroups);
        %end
        perms(i,:)=randperm(tottrials);
    end    
    
    for i=1:niter
        subindrand=subind(perms(i,:),:);
        linindrand=repmat(subindrand,Ntimes,1);
        linindrand(:,1)=reshape(repmat(1:Ntimes,tottrials,1),tottrials*Ntimes,1);
        linindrand=sub2ind([Ntimes,Ntrials,Ngroups],linindrand(:,1),linindrand(:,2),linindrand(:,3));
        RRranditer=NaN*ones(Ntimes,Ntrials,Ngroups,iter);
        RRranditer(linind)=RR(linindrand);
        RRrand(:,:,:,i)=RRranditer;
    end
    
%     for m=1:Ntrials
%         for k=1:Ngroups
%             for i=1:niter
%                 RRrand(:,m,k,i)=RR(:,m,perms(i,m,k));
%             end
%         end
%     end
    RRrandout=RRrand;
    %grandmeanrand is equal to grandmean
%    RRrand=RRrandin(:,:,:,(iterblocks(iter)+1):iterblocks(iter+1));
    groupmeanrand=mean(RRrand,2,'omitnan');

    ssbetweenrand=squeeze(sum(bsxfun(@times,bsxfun(@minus,squeeze(groupmeanrand),grandmean).^2,ntrials),2));
    sstotrand=squeeze(sum(sum(bsxfun(@minus,RRrand,grandmean).^2,2,'omitnan'),3,'omitnan'));
    
    sserrorrand=sstotrand-ssbetweenrand;
    Fstatrand=(ssbetweenrand./(Ngroups-1))./(sserrorrand./(sum(ntrials,2)-Ngroups));
    
    sig_ind_rand=Fstatrand>threshold;
    
    
    clusts=regionprops(sig_ind,'PixelList');
    Fclust=zeros(1,length(clusts));
    
    Fstat=sserror;
    Fstatrand=sserrorrand;
    
    for ii=1:length(clusts)
        Fclust(ii)=sum(arrayfun(@(x) Fstat(x),clusts(ii).PixelList(:,2)));
    end
    
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

pclusts=1-sum(bsxfun(@gt,Fclust',Fclustrandmax),2)/Niter;
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

end