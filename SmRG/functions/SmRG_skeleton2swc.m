function [swc] = SmRG_skeleton2swc(Jskel, initPos)
% SmRG_skeleton2swc: 
%           creates .swc file from 3D skeleton with first node defined by
%           initPos
%
% Syntax:
%           swc = SmRG_getBackground(Jskel,initPos)
%
% Inputs:
%           Vin: logical 3D array of skeleton
%           initPos = [x y z] coordinates.
%
% Outputs:
%           swc:  swc file 
%
%           IMPORTANT NOTE: the swc file is created without assignin radius
%           info (i.e. column 6 of swc file). We assume that any relevant
%           morphological information is available from the output of SmRG
%           itself, i.e. the binary J dataset obtained from
%           SmRG_regionGrowing
%
%% chek inputs
if nargin<1
    help SmRG_skeleton2swc
    return
end

if nargin > 2
    error('Wrong number of input arguments!')
end
if ~exist('Jskel', 'var')
    return
end
if ~exist('initPos', 'var') || isempty(initPos)
    return
end

[nRow, nCol, nSli] = size(Jskel);
if initPos(1) < 1 || initPos(2) < 1 ||...
        initPos(1) > nRow || initPos(2) > nCol
    error('Initial position out of bounds, please try again!')
end
if Jskel(initPos(2),initPos(1),initPos(3))~=0
    error('Initial position must be a point belonging to the skeleton')
end

%% initialize  
% preallocate J
J=zeros(size(Jskel));

% add the initial pixel to the queue
queue = [initPos(1), initPos(2), initPos(3)];
J(initPos(1), initPos(2), initPos(3))= 1;

% initialize swc file
swc = zeros(1,7);
swc(1,:) = [1 2 queue 1 -1];

%%% START %%%
supercount = 1;
parentcount= 0;


while size(queue, 1)
    % the first queue position determines the new values
    xv = queue(1,1);
    yv = queue(1,2);
    zv = queue(1,3);
    
    % .. and delete the first queue position
    queue(1,:) = [];
    count = 0;
    
    % check the neighbors for the current position
    for i = -1:1
        for j = -1:1
            for k = -1:1
                
                if xv+i > 0  &&  xv+i <= nRow &&...          % within the x-bounds?
                        yv+j > 0  &&  yv+j <= nCol &&...          % within the y-bounds?
                        zv+k > 0  &&  zv+k <= nSli &&...          % within the z-bounds?
                        any([i, j, k])       &&...      % i/j/k of (0/0/0) is redundant!
                        ~J(xv+i, yv+j, zv+k) &&...         % pixelposition already set?
                        Jskel(xv+i, yv+j, zv+k)==1
                    count = count+1;
                    % current pixel is true, if all properties are fullfilled
                    J(xv+i, yv+j, zv+k) = true;
                    
                    % add the current pixel to the computation queue (recursive)
                    queue(end+1,:) = [xv+i, yv+j, zv+k];
                    parent_(count,1)  = find(ismember(swc(:,3:5),[xv yv zv],'rows'));
                    
                end
            end
        end
    end
    
    % if one child
    if count==1
        supercount = supercount+1;
        parentcount= parentcount+1;
        swc(end+1,:) = [supercount 2 queue(end,:) 1 parent_];
    elseif count>1 % else
        supercount = supercount+count;
        parentcount=parentcount+1;
        swc(end+1:end+count,:) = [[supercount-count+1:supercount]' 2*ones(count,1) queue(end-count+1:end,:) ones(count,1) parent_];
    end
    parent_ = [];
end

% re-order x and y coordinates. One day I will understand the difference
swc = swc(:,[1 2 4 3 5 6 7]);
