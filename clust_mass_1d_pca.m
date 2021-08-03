function [hmat,pout,clusts,Fclust,Fclustrandmax]=clust_mass_1d_pca(RR,RRrand,p,pclust)
% This function performs the cluster mass test to test whether multiple
% conditions (Ngroups) have distinct population firing rates
% along particular dimensions (Ndimensions) in state space during adjacent
% points in time. It was used in this paper:
% https://www.biorxiv.org/content/10.1101/2020.11.11.378224v1. A companion
% mansucript will describe statistical properties, motivation, and
% generalization of this test.
% 
% Inputs:
% 
% RR: Ntimes x Ndimensions x Ngroups array
%   Population responses along each state space dimension analyzed
%   marginalized by condition
%
% RRrandA: Ntimes x Ndimensions x Ngroups x Niter array:
%   Permutation state space responses
%
% p: scalar
%   p-value that defines the threshold to define clusters
%
% pclust: scalar
%   p-value that defines significance for each cluster
%
% Outputs:
% h: Ntimes x 3 logical matrix
%   True for time points in significant clusters, false otherwise
%
% pout: Ntimes x 3 matrix
%   p-value for each point in clusters that are tested. These are the
%   p-values of each cluster, not separate p-values for each time point.
%   Time points not in a cluster are set to NaN.
%
% clusts: Nclusts x 3 struct
%   Structure containing p-values and time-points of each cluster for each
%   main effect and the interaction effect.
% 
% Fclust: Nclusts x 1 vector
%   Cluster statistic for each cluster
%
% Fclustrandmax: Niter x 1 vector
%   Null distribution of cluster statistics

Ntimes=size(RR,1);
Ndims=size(RR,2);
Ngroups=size(RR,3);
Niter=size(RRrand,4);

grandmean=mean(RR,3,'omitnan');

ss=sum(sum((bsxfun(@minus,RR,grandmean)).^2,3),2);

ssrand=squeeze(sum(sum((bsxfun(@minus,RRrand,grandmean)).^2,3),2));

ssrandsort=sort(ssrand,2);

threshold=ssrandsort(:,round((1-p)*Niter));

sig_ind=ss>threshold;

sig_ind_rand=bsxfun(@gt,ssrand,threshold);

Fstat=ss;
Fstatrand=ssrand;

Fclustrandmax=NaN*ones(1,Niter);
pout=NaN.*ones(1,Ntimes);
hmat=NaN*ones(1,Ntimes);
    
clusts=regionprops(sig_ind,'PixelList');
Fclust=zeros(1,length(clusts));
clustscell=cell(length(clusts),1);

if size(clusts,1)>0
for ii=1:length(clusts)
    Fclust(ii)=sum(arrayfun(@(x) Fstat(x),clusts(ii).PixelList(:,2)));
    clustscell{ii}=clusts(ii).PixelList(:,1);
end

for i=1:Niter
    
    clusts_rand=regionprops(sig_ind_rand(:,i),'PixelList');
    
    if ~isempty(clusts_rand)
        Fclust_rand=zeros(1,length(clusts_rand));
        for ii=1:length(clusts_rand)
            Fclust_rand(ii)=sum(arrayfun(@(x) Fstatrand(x,i),clusts_rand(ii).PixelList(:,2)));
        end
        Fclustrandmax(i)=max(Fclust_rand);
    else
        Fclustrandmax(i)=0;
    end
    
end

pclusts=1-sum(bsxfun(@gt,Fclust',Fclustrandmax),2)/Niter;
else
    pclusts=[];
end

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
hmat=h;

end
