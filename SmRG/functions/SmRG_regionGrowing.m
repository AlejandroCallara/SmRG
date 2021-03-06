function [P, J] = SmRG_regionGrowing(cIM, initPos, tfFillHoles, tfSimplify,background,p_th,BoolDistSeed,use_mex)
% SmRG_regionGrowing:
%           Region growing algorithm routine based on Daniel Kellner's
%           code available at:
%
% https://it.mathworks.com/matlabcentral/fileexchange/32532-region-growing--2d-3d-grayscale-?focused=5195969&tab=function
%
% Syntax:
%           [P, J]=SmRG_regionGrowing(cIM, initPos, tfFillHoles,...
%                   tfSimplify,background,p_th)
%
% Inputs:
%           cIM:            3D grayscale image
%           initPos:        3 element array with seeds coordinates
%           tfFillHoles:    Fills enclosed holes in the binary mask  {true}
%           tfSimplify:     Reduces the number of vertices {true, if dpsimplify exists}
%           background:     max value of background
%           p_th:           probability threshold for mixture model
%           BoolDistSeed    if true computes new seeds from regional maxima
%                           of distance trasform, otherwise all segmented
%                           pixels are used as new seeds
%           use_mex         if true attemps at using mex version of
%                           SmRG_mixtureModel. At the moment the code is
%                           given only with mexw64.
%
% Outputs:
%           P:  VxN array (with V number of vertices, N number of dimensions)
%               P is the enclosing polygon for all associated pixel/voxel
%           J:  Binary mask (with the same size as the input image) indicating
%               1 (true) for associated pixel/voxel and 0 (false) for outside
%
%% check inputs

if nargin<1
    help SmRG_regionGrowing
    return
end

% too many inputs?
if nargin >8
    error('Wrong number of input arguments!')
end

% check for dataset
if ~exist('cIM', 'var')
    cIM = SmRG_importData;
end

% check for seed
if ~exist('initPos', 'var') || isempty(initPos)
    seed_tmp = SmRG_getSeed(cIM);
    initPos  = permute(seed_tmp,[2,1,3]);
end

% check for fill holes
if ~exist('tfFillHoles', 'var')
    tfFillHoles = false;
end

% check image dims
if isequal(ndims(cIM), 2)
    initPos(3) = 1;
elseif isequal(ndims(cIM),1) || ndims(cIM) > 3
    error('There are only 2D images and 3D image sets allowed!')
end

% get image dims
[nRow,nCol,nSli] = size(cIM);

% check for simplify option
if ~isempty(which('dpsimplify.m'))
    if ~exist('tfSimplify', 'var')
        tfSimplify = true;
    end
    simplifyTolerance = 1;
else
    tfSimplify = false;
end

% check for background
if ~exist('background', 'var')
    zseed = initPos(1,3);
    background = SmRG_getBackground(cIM(:,:,zseed));
end

% check for p_th
if ~exist('p_th','var')
    p_th = 95;
end

% check for BoolDistSeed
if ~exist('BoolDistSeed', 'var')
    BoolDistSeed = true;
end

% check for use_mex
if ~exist('use_mex', 'var')
    str = which('SmRG_mixtureModelFitting_mex');
    if ~isempty(str) && ispc
        use_mex = true;
    else
        use_mex = false;
    end
end

% text output with initial parameters
disp(['RegionGrowing Opening: Initial position (' num2str(initPos(1))...
    '|' num2str(initPos(2)) '|' num2str(initPos(3)) ')'])

% preallocate array
J = false(nRow, nCol, nSli);
cIMtmp=zeros(nRow,nCol,nSli);

% REGION GROWING ALGORITHM
%% CROP
% crop size is N/8xN/8x3 by default. It can be changed by changing the
% crop_delta variable. the crop_delta variable acts on the plane crop size
% by adding and subtracting its value from the seed.
%
% i.e:
%       xseed-crop_delta:xseed+crop_delta-1
%       yseed-crop_delta:yseed+crop_delta-1
%       zseed-1:zseed+1

% get dimension
crop_delta_row=64;%round(nRow/16,0);
crop_delta_col=64;%round(nCol/16,0);

% set square crop
crop = [crop_delta_row,crop_delta_col];
[crop_min,i_crop_min]=min(crop);

crop_delta=crop(i_crop_min);

% if to few points force crop to minimum allowed value crop_delta = 32
if crop_delta <32
    crop_delta = 32;
end


xv=initPos(1);yv=initPos(2);zv=initPos(3);
queue = [xv, yv, zv];
bool_x=false; bool_y=false; bool_z=false;
xinit = xv-crop_delta; if xinit<=1; xinit=1; bool_x=true;end
yinit = yv-crop_delta; if yinit<=1; yinit=1; bool_y=true;end
zinit = zv-1;  if zinit<=1; zinit=1; bool_z=true;end

% if out of bounds resize crop
if bool_x
    xend = xinit+2*crop_delta-1;
else
    xend = xv+crop_delta-1;
end
if xend>=nRow
    xend = nRow;
    xinit=xend-(2*crop_delta-1);
end
if bool_y
    yend = yinit+2*crop_delta-1;
else
    yend = yv+crop_delta-1;
end
if yend>=nCol
    yend = nCol;
    yinit=yend-(2*crop_delta-1);
end
if bool_z
    zend = zinit+2;
else
    zend = zv+1;
end
if zend>=nSli
    zend = nSli;
    zinit=zend-2;
end

% sort of seeds depending on initPos
if zinit<round(nSli/2)
    ordine='ascend';
else
    ordine='descend';
end
v_init = [xinit yinit zinit];
%% Hartigan's dip statistic
% create null for Hartigans' dip test: boot_dip will
% be used intensively during the region growing procedure
nboot =10000;
boot_dip = SmRG_nullOfUnifDist(nboot);

Vx3 = cIM(xinit:xend,yinit:yend,zinit:zend);

% dip stuff
Vdip= Vx3(:,:,2);
[row,col] = ind2sub(size(Vdip),find(Vdip==4095));

% set to nan saturated pixels, otherwise statistics may fail due to a
% "spike" on the right tail of the histogram
for ii=1:length(row)
    Vdip(row(ii),col(ii))=0/0;
end
Ne = sort(Vdip(:));
Ne = Ne(~ismissing(Ne));

% dip statistic of current crop
dip  = HartigansDipTest(Ne);
p    = round(sum(dip<boot_dip)/nboot,4);

%% check on bimodality

mask = false(size(Vx3));

% if bimodal OTSU
if p<=0.01
    if max(Vx3(:))<background
        disp('Only background -> skipping segmentation'), disp('');
        
        % set to zero values in crop (skip segmentation)
        cIMtmp(xinit:xend,yinit:yend,zinit:zend)=(zeros(size(Vx3)));
    else
        disp('the distribution was found to be bimodal -> segmenting with Otsu''s method'), disp('');
        
        % segment with Otsu's method
        im_tmp = (Vx3-min(Vx3(:)))/(max(Vx3(:))-min(Vx3(:)));
        for dc = 1:3
            sss= graythresh(im_tmp(:,:,dc));
            mask(:,:,dc)=im2bw(im_tmp(:,:,dc),sss);
        end
        mask = SmRG_pruneMask(mask, queue,v_init);
        mask = or(mask,cIMtmp(xinit:xend,yinit:yend,zinit:zend)>0);
        
        cIMtmp(xinit:xend,yinit:yend,zinit:zend)=mask.*Vx3;
    end
    % else EM
else
    disp('the distribution was found to be unimodal -> segmenting with mixture model fitting'), disp('');
    
    % get central slice of crop
    V_central_slice= Vx3(:,:,2);
    vec_cs = V_central_slice(:);
    
    %set initial values for model fitting: gaussian
    K0_double =2*(min(vec_cs))+5; % mean
    vB_double = 2*K0_double;    % variance
    
    %set initial values for model fitting: negbin
    mu_sk = 2*mean(vec_cs,'omitnan'); % mean
    v_sk = 2*mu_sk; % variance
    rk = mu_sk^2/(v_sk-mu_sk);
    
    
    
    % fit model
    if use_mex
        [p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost_mex(V_central_slice,...
            K0_double,vB_double,mu_sk,rk);
    else
        [p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost(V_central_slice,...
            K0_double,vB_double,mu_sk,rk);
    end
    
    % get threshold
    d = 0;
    p_th = 100;
    count = 0;
    p_tot(p_tot>1)=1;
    n_soglie = 2;
    trovate_2 = true;
    count_soglie = 0;
    while trovate_2
        p_th_dist = prctile(p_tot,p_th);
        d = 1- p_th_dist;
        count = count+1;
        if d~=0
            count_soglie = count_soglie+1;
            soglie(count_soglie)=p_th;
            p_dist_soglie(count_soglie)=p_th_dist;
            n_soglie = n_soglie-1;
        end
        if n_soglie == 0
            trovate_2 = false;
        end
        p_th = p_th-1;
    end
    p_th_dist = p_dist_soglie(1);
    % p_th_dist = 1;
    
    % old way to get threshold: still good
    % p_th_dist = prctile(p_tot,p_th);
    
    % get indices that satisfies condition
    i_signal = (1-p_tot)<=(1-p_th_dist);
    
    % find indices in image
    tmp = vec_cs(i_signal);
    tmp(tmp==0)=max(tmp);
    MIN = min(tmp);
    if isempty(MIN)
        mask = false(size(Vx3));
    else
        mask= Vx3>MIN;
        % do not remove what already segmented
        mask = SmRG_pruneMask(mask, queue,v_init);
        mask=or(mask,cIMtmp(xinit:xend,yinit:yend,zinit:zend)>0);
        
    end
    V_signal=double(mask).*Vx3;
    
    % set values in segmenting crop
    cIMtmp(xinit:xend,yinit:yend,zinit:zend)= V_signal;
    
end

% initialize queues
queue = [xv,yv,zv];
seed_check = [xv,yv,zv];
seed  = seed_check;
seed_fatti=[-1 -1 -1];

while size(seed,1)
    while size(queue, 1)
        
        % the first queue position determines the new values
        xv = queue(1,1); yv = queue(1,2); zv = queue(1,3);
        
        % .. and delete the first queue position
        queue(1,:) = [];
        
        % check the neighbors for the current position
        for ix = -1:1
            for jy = -1:1
                for kz = -1:1
                    if xv+ix > 0  &&  xv+ix <= nRow &&...     % within the x-bounds?
                            yv+jy > 0  &&  yv+jy <= nCol &&...% within the y-bounds?
                            zv+kz > 0  &&  zv+kz <= nSli &&...% within the z-bounds?
                            any([ix, jy, kz])       &&...     % i/j/k of (0/0/0) is redundant!
                            ~J(xv+ix, yv+jy, zv+kz) &&...     % pixelposition already set?
                            cIMtmp(xv+ix, yv+jy, zv+kz)~=0
                        
                        % current pixel is true, if all properties are fullfilled
                        J(xv+ix, yv+jy, zv+kz) = true;
                        
                        % add the current pixel to the computation queue (recursive)
                        queue(end+1,:) = [xv+ix, yv+jy, zv+kz];
                        
                    end
                end
            end
        end
    end
    %%%% start: OLD STUFF, maybe one day still useful
    %%%  cIMtmp = cIMtmp.*J;
    %%%% end: OLD stuff
    clear tmp2 tmp3 seed_to_rem
    
    seed_fatti(end+1,:)   = seed(1,:);
    seed(1,:)    = [];
    
    %     ccmask = bwconncomp(mask);
    %     for icc=1:length(ccmask.PixelIdxList)
    %         [cc1, cc2, cc3]=(ind2sub(size(mask),ccmask.PixelIdxList{1,icc}));
    %         SEMIcc= [cc1+xinit-1 cc2+yinit-1 cc3+zinit-1];
    %         if sum(ismember(SEMIcc,seed_fatti(end,:),'rows'))
    %             mask = false(size(mask));
    %             mask(ccmask.PixelIdxList{1,icc})=true;
    %         end
    %     end
    % find mask-slices with at least one 1
    s = sum(sum(mask));
    s = s(:);
    
    sfind = find(s);
    MINi = min(sfind);
    MAXi = max(sfind);
    
    % search for new seeds
    if ~BoolDistSeed
        [uot1, uot2, uot3]=ind2sub(size(mask),find(mask));
        SEMI= [uot1 uot2 uot3];
    else
        SEMI = [];
        iz  = [];
        for d3 = MINi:MAXi
            
            % compute distance transform
            Jdist = round(bwdist(~mask(:,:,d3)),0);
            
            % get regional maxima
            Jtmp  = imregionalmax(Jdist,8);
            
            % get new seeds as centroids of regional maxima
            sz      = regionprops(Jtmp);
            for i_sz = 1:length(sz)
                tmp_iz = [flip(round(sz(i_sz).Centroid)) d3];
                iz = [iz;tmp_iz];
            end
            
            % add boundaries
            % Jbound = bwperim(J(:,:,d3));
            % [boundx,boundy,boundz]=ind2sub(size(Jbound),find(Jbound));
            % SEMIbound= [boundx,boundy,ones(length(boundx),1)*d3];
            % SEMI  = cat(1,SEMI,iz,SEMIbound);
            SEMI  = cat(1,SEMI,iz);
            
            iz = [];
            
            
        end
    end
    if ~isempty(SEMI)
        SEMI = SEMI+[xinit-1 yinit-1 zinit-1];
    end
    seed = [seed;SEMI];
    matched_seeds = ismember(seed,seed_fatti,'rows');
    
    seed(matched_seeds,:)=[];
    seed=unique(seed,'rows');
    
    % sort seeds...maybe not so useful
    % [~,sort_ind]=sort(seed(:,3),ordine);
    % seed=seed(sort_ind,:);
    
    clear solo_x solo_y solo_z
    size(seed)
    if isempty(seed);break;end
    xseed = seed(1,1);
    yseed = seed(1,2);
    zseed = seed(1,3);
    queue =[xseed,yseed,zseed];
    
    %% CROP
    bool_x=false; bool_y=false; bool_z=false;
    xinit = xseed-crop_delta; if xinit<=1; xinit=1; bool_x=true;end
    yinit = yseed-crop_delta; if yinit<=1; yinit=1; bool_y=true;end
    zinit = zseed-1;  if zinit<=1; zinit=1; bool_z=true;end
    
    if bool_x
        xend = xinit+2*crop_delta-1;
    else
        xend = xseed+crop_delta-1;
    end
    if xend>=nRow
        xend = nRow;
        xinit=xend-(2*crop_delta-1);
    end
    if bool_y
        yend = yinit+2*crop_delta-1;
    else
        yend = yseed+crop_delta-1;
    end
    if yend>=nCol
        yend = nCol;
        yinit=yend-(2*crop_delta-1);
    end
    if bool_z
        zend = zinit+2;
    else
        zend = zseed+1;
    end
    if zend>=nSli
        zend = nSli;
        zinit=zend-2;
    end
    v_init = [xinit yinit zinit];
    %% Hartigan's dip statistic
    % null already created
    Vx3 = cIM(xinit:xend,yinit:yend,zinit:zend);
    
    % dip stuff
    Vdip= Vx3(:,:,2);
    [row,col] = ind2sub(size(Vdip),find(Vdip==4095));
    
    % set to nan saturated pixels, otherwise statistics may fail due to a
    % "spike" on the right tail of the histogram
    for ii=1:length(row)
        Vdip(row(ii),col(ii))=0/0;
    end
    Ne = sort(Vdip(:));
    Ne = Ne(~ismissing(Ne));
    
    % dip statistic of current crop
    dip  = HartigansDipTest(Ne);
    p    = round(sum(dip<boot_dip)/nboot,4);
    
    %% check on bimodality
    if p<=0.01
        if max(Vx3(:))<background
            disp('Only background -> skipping segmentation'), disp('');
            
            % set to zero values in crop (skip segmentation)
            mask = or(false(size(Vx3)),cIMtmp(xinit:xend,yinit:yend,zinit:zend)>0);
            cIMtmp(xinit:xend,yinit:yend,zinit:zend)=mask.*Vx3;
        else
            disp('the distribution was found to be bimodal -> segmenting with Otsu''s method'), disp('');
            
            % segment with Otsu's method
            im_tmp = (Vx3-min(Vx3(:)))/(max(Vx3(:))-min(Vx3(:)));
            for dc = 1:3
                sss= graythresh(im_tmp(:,:,dc));
                mask(:,:,dc)=im2bw(im_tmp(:,:,dc),sss);
            end
            %cIMtmp(xinit:xend,yinit:yend,zinit:zend)=binV.*Vx3;
            
            mask = SmRG_pruneMask(mask,  queue,v_init);
            mask = or(mask,cIMtmp(xinit:xend,yinit:yend,zinit:zend)>0);
            
            cIMtmp(xinit:xend,yinit:yend,zinit:zend)=mask.*Vx3;
        end
        % else EM
    else
        disp('the distribution was found to be unimodal -> segmenting with mixture model fitting'), disp('');
        
        % get central slice of crop
        V_central_slice= Vx3(:,:,2);
        vec_cs = V_central_slice(:);
        
        % reset mask
        % mask = zeros(size(V_central_slice));
        %
        % %%% test hist distances added by Nicola Vanello.
        % threshold = 10^-2;
        % if exist('stored_init_cond')
        %   proposed = zeros(1,4);
        % else
        %   stored_init_cond = [];
        %   proposed = [];%zeros(1,4);
        %   sample_values = [];
        % end
        % actual_values = vec_cs;
        % [sample_values,found,stored_init_cond,proposed]=...
        %                 SmRG_findCloser(actual_values,sample_values,stored_init_cond,...
        %                 proposed,threshold);
        %
        %
        %
        % %%%fit model added by Nicola Vanello: evaluates the distance between
        % histograms to use optimal priors. Actually not used.
        % if found
        %    K0_double   = proposed(1);
        %    vB          = proposed(2);
        %    mu_sk       = proposed(3);
        %    rk          = proposed(4);
        %    [p_tot,a,K0_double,vB,mu_sk,rk] = SmRG_mixtureModelFitting(V_central_slice,K0_double,vB,mu_sk,rk);
        % else
        %    [p_tot,a,K0_double,vB,mu_sk,rk] = SmRG_mixtureModelFitting(V_central_slice);
        %    sample_values = [sample_values;vec_cs'];
        %    stored_init_cond = [stored_init_cond;[K0_double,vB,mu_sk,rk]];
        % end
        
        % fit model
        if exist('K0_double')
            if use_mex
                [p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost_mex(V_central_slice,K0_double,vB_double,mu_sk,rk);
            else
                [p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost(V_central_slice,K0_double,vB_double,mu_sk,rk);
            end
        else
            %set initial values for model fitting: gaussian
            K0_double =2*(min(vec_cs))+5; % mean
            vB_double = 2*K0_double;    % variance
            
            %set initial values for model fitting: negbin
            mu_sk = 2*mean(vec_cs,'omitnan'); % mean
            v_sk = 2*mu_sk; % variance
            rk = mu_sk^2/(v_sk-mu_sk);
            if use_mex
                [p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost_mex(V_central_slice,K0_double,vB_double,mu_sk,rk);
            else
                [p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_newPost(V_central_slice,K0_double,vB_double,mu_sk,rk);
            end
        end
        
        % %useful in debug mode: plot the model fitting on data hist
        %figure(1), histogram(V_central_slice)
        %hold on
        %plot(50000*normpdf(1:max(V_central_slice(:)),K0_double,sqrt(vB_double)))
        %plot(100000*nbinpdf_mu(1:max(V_central_slice(:)),mu_sk,1/rk))
        %drawnow
        %hold off
        
        % set threshold
        % get threshold
        d = 0;
        p_th = 100;
        count = 0;
        p_tot(p_tot>1)=1;
        n_soglie = 2;
        trovate_2 = true;
        count_soglie = 0;
        while trovate_2
            p_th_dist = prctile(p_tot,p_th);
            d = 1- p_th_dist;
            count = count+1;
            if d~=0
                count_soglie = count_soglie+1;
                soglie(count_soglie)=p_th;
                p_dist_soglie(count_soglie)=p_th_dist;
                n_soglie = n_soglie-1;
            end
            if n_soglie == 0
                trovate_2 = false;
            end
            p_th = p_th-0.5;
        end
        p_th_dist = p_dist_soglie(1);
%         p_th_dist = 1;
        
        % old way to get threshold: still good
        % p_th_dist = prctile(p_tot,p_th);
        
        % get indices that satisfies condition
        i_signal = (1-p_tot)<=(1-p_th_dist);
        
        % find indices in image
        tmp = vec_cs(i_signal);
        tmp(tmp==0)=max(tmp);
        MIN = min(tmp);
        if isempty(MIN)
            mask = false(size(Vx3));
        else
            mask= Vx3>MIN;
            % do not remove what already segmented
            
            mask = SmRG_pruneMask(mask,  queue,v_init);
            mask = or(mask,cIMtmp(xinit:xend,yinit:yend,zinit:zend)>0);
        end
        V_signal=double(mask).*Vx3;
        
        % set values in segmenting crop
        cIMtmp(xinit:xend,yinit:yend,zinit:zend)= V_signal;
    end
end
%%% END OF REGION GROWING ALGORITHM

% loop through each slice, fill holes and extract the polygon vertices
P = [];
for cSli = 1:nSli
    if ~any(J(:,:,cSli))
        continue
    end
    
    % use bwboundaries() to extract the enclosing polygon
    if tfFillHoles
        % fill the holes inside the mask
        J(:,:,cSli) = imfill(J(:,:,cSli), 'holes');
        B = bwboundaries(J(:,:,cSli), 8, 'noholes');
    else
        B = bwboundaries(J(:,:,cSli));
    end
    
    newVertices = [B{1}(:,2), B{1}(:,1)];
    
    % simplify the polygon via Line Simplification
    if tfSimplify
        newVertices = dpsimplify(newVertices, simplifyTolerance);
    end
    
    % number of new vertices to be added
    nNew = size(newVertices, 1);
    
    % append the new vertices to the existing polygon matrix
    if isequal(nSli, 1) % 2D
        P(end+1:end+nNew, :) = newVertices;
    else                % 3D
        P(end+1:end+nNew, :) = [newVertices, repmat(cSli, nNew, 1)];
    end
end

% text output with final number of vertices
disp(['RegionGrowing Ending: Found ' num2str(length(find(J)))...
    ' pixels within the threshold range (' num2str(size(P, 1))...
    ' polygon vertices)!'])
