    % test example of SmRG that calls almost all functions contained in the
% SmRG package available at: 
% 
% When running this script a window will pop up asking for a .tiff stack file. 
%
% once the dataset is imported in matlab a new window will pop up asking
% for the user to chose the seed selection mode among manual, 
% semi-automatic or automatic.
%
% case 'manual': a window will pop up asking for the user to give manually
% the coordinates of the seed by clicking on the image. The coordinates of 
% the point given by the user will be used as starting seed for SmRG.
%
% case 'semi-automatic': a window will pop up asking for the user to select
% a slice for searching for seeds. Seeds are obtained by 2D Hough's
% trasform applied on selected slice. 
%
% case 'automatic': 3D Houhgh's transform will be applied on data searching
% for spherical objects within the dataset (potentially cell soma). For
% each sphere a seed is chosen as the center of the sphere. 

% clear
% close all
% clc
% % DATASET IMPORT
% tic
% % import data and keep memory of file name and file path
% [V, name_im, path_im] = SmRG_importData;
% V = medfilt3(V);

% get size of data
[nx,ny,nz]=size(V);

%% seed selection mode
str_mode = SmRG_seedSelectionMode;

%% seed selection
switch str_mode
    case 'none'
        
        % do nothing
        disp('mode not selected')
        
        % Jcell stores the output of Region Growing algorithm
        Jcell = {};
        
    case 'manual'
        
        disp('select seed and press Enter')
        
        % Manual seed input
        seed = SmRG_getSeed(V);
        
        % Jcell stores the output of Region Growing algorithm
        Jcell = cell(1,2);
        
    case 'semi-automatic'
        
        disp('select slice on which search for seeds')
        
        % Semi-Automatic
        % scrolls dataset until slice is chosen
        zseed = SmRG_scrollDataset(V);
        Iseed = V(:,:,zseed);
        
        % find circles on selected slice by 2D Hough's transform
        [c,r] = imfindcircles(Iseed, [15,45]);
        while isempty(c)==1
            zseed = SmRG_scrollDataset(V);
            Iseed = V(:,:,zseed);
            [c,r] = imfindcircles(Iseed, [15 ,45]);
        end
        U     = size(c,1);
        
        % Jcell stores the output of Region Growing algorithm
        Jcell = cell(U,2);
        
        figure;
        imagesc(Iseed);
        viscircles(c, r,'EdgeColor','r');
        disp('seed detection finished');
        cent = c;
    case 'automatic'

        fprintf...
            ('implementing Spherical Hough transform on data \nthis may take a while')
        
        % Automatic
        % find spheres in V by 3D Hough's transform
        
        % inputs for SphericalHough: values are chosen to fit properly on 
        % test data.
        % For customize call to SphericalHough please refer to
        % SpehircalHough help. 
        radrange = [20 25]; % radius range
        grdthres = 0.3;     % gradient threshold
        fltrLM_R = 8;       % filter radius
        multirad = 1;       % multiple radius detection threshold
        obj_cint = 0.5;     % object center intensity threshold
        
        % some filtering to data may help SphericalHough call
        Vtmp = medfilt3(V,[11, 11, 11]);
        
        % search for spheres
        tic
        [center_img,sphere_img,cent,radi]=SphericalHough(Vtmp,radrange,grdthres,fltrLM_R,multirad,obj_cint);
        toc
        U     = size(cent,1);
        
        % visualize potential seeds. 
        SmRG_visSeed(V,sphere_img);
           
        disp('seed detection finished');
        
        
        
        % Jcell stores the output of Region Growing algorithm
        Jcell = cell(U,2);
end

%% THRESHOLD CALCULATION AND REGION GROWING
prob_thresh= 100;             % for testing with PCs
% prob_thresh= 50;           % for testing with OP fibers
BoolDistSeed = true;
l = size(Jcell,1);
if l==1
    xseed = round(seed (1,2));
    yseed = round(seed (1,1));
    zseed = round(seed (1,3));
    
    % get background
    background = SmRG_getBackground(V(:,:,zseed),'otsu');
    
    % region growing
    tic
    [P, J] = SmRG_regionGrowing(V, [xseed,yseed,zseed],1,1,...
        background,prob_thresh,BoolDistSeed);
    elapsedTime = toc
    Jcell{1,1} = J;
    Jcell{1,2} = P;
    % text output with final number of vertices
    disp(['RegionGrowing Ending: Found ' num2str(length(find(J)))...
    ' pixels within the threshold range (' num2str(size(P, 1))...
    ' polygon vertices)!'])
elseif l>1
    disp('Multiple seeds were found')
    str_txt = strcat('Do you want to proceed with n=',num2str(l),' seeds? [y/n]')
    str_in=input(str_txt,'s');

    switch str_in
        case 'y'
            for abc = 1:U;
                
                xseed = round(cent (abc,2));
                yseed = round(cent (abc,1));
                
                % if automatic 
                if strcmp(str_mode,'automatic')
                  zseed = round(cent (abc,3));
                end
                % get background
                background = SmRG_getBackground(V(:,:,zseed));
                
                % region growing
                [P, J] = SmRG_regionGrowing(V, [xseed,yseed,zseed],1,1,...
                    background,prob_thresh);
                Jcell{abc,1} = J;
                Jcell{abc,2} = P;
            end
        case 'n'
            disp('thanks')
    end
    
elseif l<1
    disp('thanks')
end
toc
