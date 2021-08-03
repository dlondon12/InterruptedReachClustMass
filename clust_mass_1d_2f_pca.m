function [hmat,pout,clusts,hposthoc,poutposthoc,clustspair]=clust_mass_1d_2f_pca(RR,RRrandA,RRrandB,RRrandAB,RRrandpairA,p,pclust)
% This function performs the cluster mass test to test whether multiple
% conditions (NgroupsA, Ngroups B) have distinct population firing rates
% along particular dimensions (Ndimensions) in state space during adjacent
% points in time. It was used in this paper:
% https://www.biorxiv.org/content/10.1101/2020.11.11.378224v1. A companion
% mansucript will describe statistical properties, motivation, and
% generalization of this test.
% 
% Inputs:
% 
% RR: Ntimes x Ndimensions x NgroupsA x Ngroups B array
%   Population responses along each state space dimension analyzed
%   marginalized by condition
%
% RRrandA, RRrandB, RRrandAB, RRrandpairA: Ntimes x Ndimensions x NgroupsA
% x Ngroups B x Niter arrays:
%   Permutation state space responses for each main effect and interaction
%   effect. RRrandpairA is used for post-hoc testing which was not used in
%   the paper. Interaction effect was not used in the paper.
%
% p: scalar
%   p-value that defines the threshold to define clusters
%
% pclust: scalar
%   p-value that defines significance for each cluster
%
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


Ntimes=size(RR,1);
Ndims=size(RR,2);
NgroupsA=size(RR,3);
NgroupsB=size(RR,4);
Niter=size(RRrandA,5);

grandmean=mean(reshape(RR,Ntimes,Ndims,NgroupsA*NgroupsB),3,'omitnan');
groupAmean=mean(RR,4,'omitnan');
groupBmean=mean(RR,3,'omitnan');

ssA=NgroupsB*sum(sum((bsxfun(@minus,groupAmean,grandmean)).^2,3),2);
ssB=NgroupsA*sum(sum((bsxfun(@minus,groupBmean,grandmean)).^2,4),2);
ssAB=sum(sum(reshape((bsxfun(@plus,bsxfun(@minus,bsxfun(@minus,RR,groupAmean),groupBmean),grandmean)).^2,Ntimes,Ndims,NgroupsA*NgroupsB),3,'omitnan'),2,'omitnan');

groupAmeanrandA=mean(RRrandA,4,'omitnan');
groupBmeanrandB=mean(RRrandB,3,'omitnan');
groupAmeanrandAB=mean(RRrandAB,4,'omitnan');
groupBmeanrandAB=mean(RRrandAB,3,'omitnan');

ssArandA=squeeze(NgroupsB*sum(sum((bsxfun(@minus,groupAmeanrandA,grandmean)).^2,3),2));
ssBrandB=squeeze(NgroupsA*sum(sum((bsxfun(@minus,groupBmeanrandB,grandmean)).^2,4),2));
ssABrandAB=squeeze(sum(sum(sum((bsxfun(@plus,bsxfun(@minus,bsxfun(@minus,RRrandAB,groupAmeanrandAB),groupBmeanrandAB),grandmean)).^2,4,'omitnan'),3,'omitnan'),2,'omitnan'));

ssArandAsort=sort(ssArandA,2);
ssBrandBsort=sort(ssBrandB,2);
ssABrandABsort=sort(ssABrandAB,2);

thresholdA=ssArandAsort(:,round((1-p)*Niter));
thresholdB=ssBrandBsort(:,round((1-p)*Niter));
thresholdAB=ssABrandABsort(:,round((1-p)*Niter));

sig_indA=ssA>thresholdA;
sig_indB=ssB>thresholdB;
sig_indAB=ssAB>thresholdAB;

sig_ind_randA=bsxfun(@gt,ssArandA,thresholdA);
sig_ind_randB=bsxfun(@gt,ssBrandB,thresholdB);
sig_ind_randAB=bsxfun(@gt,ssABrandAB,thresholdAB);

Fstat=[ssA';ssB';ssAB'];
Fstatrand=permute(cat(3,ssArandA,ssBrandB,ssABrandAB),[3,1,2]);

sig_ind=[sig_indA';sig_indB';sig_indAB'];
sig_ind_rand=permute(cat(3,sig_ind_randA,sig_ind_randB,sig_ind_randAB),[3,1,2]);

Fclustrandmax=NaN*ones(3,Niter);
pout=NaN.*ones(3,Ntimes);
hmat=NaN*ones(3,Ntimes);
for j=1:3
    
    clusts{j}=regionprops(sig_ind(j,:),'PixelList');
    Fclust=zeros(1,length(clusts{j}));
    clustscell=cell(length(clusts{j}),1);
    
    if size(clusts{j},1)>0
        for ii=1:length(clusts{j})
            Fclust(ii)=sum(arrayfun(@(x) Fstat(j,x),clusts{j}(ii).PixelList(:,1)));
            clustscell{ii}=clusts{j}(ii).PixelList(:,1);
        end
        
        for i=1:Niter
            
            clusts_rand=regionprops(sig_ind_rand(j,:,i),'PixelList');
            
            if ~isempty(clusts_rand)
                Fclust_rand=zeros(1,length(clusts_rand));
                for ii=1:length(clusts_rand)
                    Fclust_rand(ii)=sum(arrayfun(@(x) Fstatrand(j,x,i),clusts_rand(ii).PixelList(:,1)));
                end
                Fclustrandmax(j,i)=max(Fclust_rand);
            else
                Fclustrandmax(j,i)=0;
            end
            
        end
        
        pclusts=1-sum(bsxfun(@gt,Fclust',Fclustrandmax(j,:)),2)/Niter;
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

pairs=nchoosek(1:NgroupsA,2);

for r=1:3
    RRpair=RR(:,:,pairs(r,:),:);
    RRpair=RRpair(logical(hmat(1,:)),:,:,:);
    RRpair=squeeze(mean(RRpair,4,'omitnan'));
    
    RRrandpairA_temp=RRrandpairA(:,:,:,:,:,r);
    RRrandpairA_temp=RRrandpairA_temp(logical(hmat(1,:)),:,:,:,:);
    RRrandpairA_temp=squeeze(mean(RRrandpairA_temp,4,'omitnan'));
    
    [~,~,clustspair{r,1},Fclustpair{r},Fclustrandmax(r,:)]=clust_mass_1d_pca(RRpair,RRrandpairA_temp,0.1,0.05);
end

Fclustrandmax=max(Fclustrandmax);
poutpair=NaN.*ones(length(pairs),nnz(hmat(1,:)));
hpair=false(length(pairs),nnz(hmat(1,:)));

for r=1:length(pairs)
    pclusts=1-sum(bsxfun(@gt,Fclustpair{r}',Fclustrandmax),2)/Niter;
    
    coords=[];
    for ii=1:length(pclusts)
        clustspair{r,1}(ii).p=pclusts(ii);
        poutpair(r,clustspair{r,1}(ii).PixelList(:,2))=pclusts(ii);
        
        if clustspair{r,1}(ii).p<pclust
            coords=[coords;clustspair{r,1}(ii).PixelList];
        end
    end
    
    if ~isempty(coords)
        htemp=sparse(coords(:,1),coords(:,2),1,1,nnz(hmat(1,:)));
        
        hpair(r,:)=logical(full(htemp));
    end
end

hposthoc{1}=false(length(pairs),size(RR,1));
poutposthoc{1}=NaN.*ones(length(pairs),size(RR,1));
hposthoc{1}(repmat(logical(hmat(1,:)),3,1))=hpair;
poutposthoc{1}(repmat(logical(hmat(1,:)),3,1))=poutpair;

for r=1:NgroupsA
    RRpair=squeeze(RR(:,:,r,:));
    RRpair=RRpair(logical(hmat(3,:)),:,:);
    
    RRrandpairB=squeeze(RRrandB(:,:,r,:,:));
    RRrandpairB=RRrandpairB(logical(hmat(3,:)),:,:,:);
    
    [~,~,clustspair{r,3},Fclustpair{r},Fclustrandmax(r,:)]=clust_mass_1d_pca(RRpair,RRrandpairB,0.1,0.05);
end

for r=1:NgroupsB
    RRpair=squeeze(RR(:,:,:,r));
    RRpair=RRpair(logical(hmat(3,:)),:,:);
    
    RRrandA_temp=squeeze(RRrandA(:,:,:,r,:));
    RRrandA_temp=RRrandA_temp(logical(hmat(3,:)),:,:,:);
    
    [~,~,clustspair{r+NgroupsA,3},Fclustpair{r+NgroupsA},Fclustrandmax(r+NgroupsA,:)]=clust_mass_1d_pca(RRpair,RRrandA_temp,0.1,0.05);
end

Ncomps=[NgroupsB.*ones(1,NgroupsA),NgroupsA.*ones(1,NgroupsB)];

Fclustrandmax=bsxfun(@times,Fclustrandmax,1./Ncomps');
Fclustpair=arrayfun(@(x,y) x{:}./y,Fclustpair,Ncomps,'Uniformoutput',0);

Fclustrandmax=max(Fclustrandmax);
poutpair=NaN.*ones(NgroupsA+NgroupsB,nnz(hmat(3,:)));
hpair=false(NgroupsA+NgroupsB,nnz(hmat(3,:)));

for r=1:(NgroupsA+NgroupsB)
    pclusts=1-sum(bsxfun(@gt,Fclustpair{r}',Fclustrandmax),2)/Niter;
    
    coords=[];
    for ii=1:length(pclusts)
        clustspair{r,3}(ii).p=pclusts(ii);
        poutpair(r,clustspair{r,3}(ii).PixelList(:,2))=pclusts(ii);
        
        if clustspair{r,3}(ii).p<pclust
            coords=[coords;clustspair{r,3}(ii).PixelList];
        end
    end
    
    if ~isempty(coords)
        htemp=sparse(coords(:,1),coords(:,2),1,1,nnz(hmat(3,:)));
        
        hpair(r,:)=logical(full(htemp));
    end
end

hposthoc{3}=false(NgroupsA+NgroupsB,size(RR,1));
poutposthoc{3}=NaN.*ones(NgroupsA+NgroupsB,size(RR,1));
hposthoc{3}(repmat(logical(hmat(3,:)),NgroupsA+NgroupsB,1))=hpair;
poutposthoc{3}(repmat(logical(hmat(3,:)),NgroupsA+NgroupsB,1))=poutpair;
end