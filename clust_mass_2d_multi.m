function [h,pout,clusts]=clust_mass_2d_multi(RR,p,pclust,Niter)

% This function performs the two-dimensional extension to the within-units 
% cluster mass test to test whether multiple conditions (Ngroups) have
% distinct firing rates during adjacent points in two dimensions. For
% example these dimensions can be time (Ntimes) and frequency (Nf). It was used in this
% paper: https://www.biorxiv.org/content/10.1101/2020.11.11.378224v1. A
% companion mansucript will describe statistical properties, motivation,
% and generalization of this test.
% 
% Inputs:
% 
% RR: Ntimes x Nf x Ntrials x Ngroups array
%   Ntrials can refer trials if all conditions occur within the same trial
%   or neurons if all conditions are tested in each neuron. In this case each
%   Ntimes x 1 x Ngroups sub-matrix refers to the PSTH of each neuron for
%   each condition.
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
% pout: Ntimes x Nf matrix
%   p-value for each point in clusters that are tested. These are the
%   p-values of each cluster, not separate p-values for each point. Points
%   not in a cluster are set to NaN.
%
% clusts: Nclusts x 1 struct
%   Structure containing p-values and time-frequency points of each cluster


Ntimes=size(RR,1);
Nf=size(RR,2);
Ntrials=size(RR,3);
Ngroups=size(RR,4);

grandmean=mean(reshape(RR,Ntimes,Nf,Ntrials*Ngroups),3);
groupmean=mean(RR,3);
unitmean=mean(RR,4);

ssbetween=Ntrials.*sum(bsxfun(@minus,squeeze(groupmean),grandmean).^2,3);
sswithin=sum(sum(bsxfun(@minus,RR,groupmean).^2,4),3);
ssunit=Ngroups.*sum(bsxfun(@minus,unitmean,grandmean).^2,3);

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

maxiter=10;

iterblocks=0:maxiter:Niter;
iterblocks(end)=Niter;

sserrorrand=NaN*ones(Ntimes,Nf,Niter);
Fstatrand=sserrorrand;
for iter=1:(size(iterblocks,2)-1)
    niter=iterblocks(iter+1)-iterblocks(iter);
RRrand=single(zeros(Ntimes,Nf,Ntrials,Ngroups,niter));

for m=1:Ntrials
    for k=1:Ngroups
        for i=1:niter
            RRrand(:,:,m,k,i)=RR(:,:,m,perms(iterblocks(iter)+i,m,k));
        end
    end
end

%grandmeanrand is equal to grandmean

groupmeanrand=mean(RRrand,3);
unitmeanrand=mean(RRrand,4);

ssbetweenrand=squeeze(Ntrials.*sum(bsxfun(@minus,squeeze(groupmeanrand),grandmean).^2,3));
sswithinrand=squeeze(sum(sum(bsxfun(@minus,RRrand,groupmeanrand).^2,4),3));
ssunitrand=squeeze(Ngroups.*sum(bsxfun(@minus,unitmeanrand,grandmean).^2,3));

sserrorrand(:,:,(iterblocks(iter)+1):iterblocks(iter+1))=sswithinrand-ssunitrand;
Fstatrand(:,:,(iterblocks(iter)+1):iterblocks(iter+1))=(ssbetweenrand./(Ngroups-1))./(sserrorrand(:,:,(iterblocks(iter)+1):iterblocks(iter+1))./((Ntrials-1)*(Ngroups-1)));
end

sig_ind_rand=Fstatrand>threshold;


clusts=regionprops(sig_ind,'PixelList');
Fclust=zeros(1,length(clusts));
clustscell=cell(length(clusts),1);

for ii=1:length(clusts)
    Fclust(ii)=sum(arrayfun(@(x,y) sserror(x,y),clusts(ii).PixelList(:,2),clusts(ii).PixelList(:,1)));
    clustscell{ii}=[clusts(ii).PixelList(:,2),clusts(ii).PixelList(:,1)];
end

for i=1:Niter
    
    clusts_rand=regionprops(sig_ind_rand(:,:,i),'PixelList');
    
    if ~isempty(clusts_rand)
        Fclust_rand=zeros(1,length(clusts_rand));
        for ii=1:length(clusts_rand)
            Fclust_rand(ii)=sum(arrayfun(@(x,y) sserrorrand(x,y,i),clusts_rand(ii).PixelList(:,2),clusts_rand(ii).PixelList(:,1)));
        end
        Fclustrandmax(i)=max(Fclust_rand);
    else
        Fclustrandmax(i)=0;
    end
    
end

pclusts=1-sum(bsxfun(@gt,Fclust',Fclustrandmax),2)/Niter;
pout=NaN.*ones(length(clusts),Nf,Ntimes);

inds=[];
vals=[];
coords=[];
for ii=1:length(pclusts)
    clusts(ii).p=pclusts(ii);
    inds=[inds;clusts(ii).PixelList(:,1),clusts(ii).PixelList(:,2)];
    vals=[vals;pclusts(ii)*ones(size(clusts(ii).PixelList,1),1)];
    
    if clusts(ii).p<pclust
        coords=[coords;clusts(ii).PixelList];
    end
end
pout=accumarray(inds,vals,[Nf,Ntimes],[],NaN);

if isempty(coords)
    h=false(Nf,Ntimes);
else
    h=sparse(coords(:,1),coords(:,2),1,Nf,Ntimes);
    
    h=logical(full(h));
end


end