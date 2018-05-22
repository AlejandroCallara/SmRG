function [V, name_im, path_im] = SmRG_importData
% SmRG_importData: 
%           import 3D dataset. The dataset must be a .tiff stack file.
%           Future updates will consider other formats. 
% 
% Syntax:
%           V = SmRG_importData
%           [V, name_im] = SmRG_importData
%           [V, name_im, path_im] = SmRG_importData
%
% Output:
%           V: 3D dataset
%           name_im: file name
%           path_im: file path

% dataset info
[name_im, path_im]= uigetfile('*.tif');
if name_im==0
    warning('no dataset selected')
    return
end
full_name         = strcat(path_im,name_im);
info_im = imfinfo(full_name);

% image size 
nz = numel(info_im);
ny = info_im.Width;
nx = info_im.Height;

% import data
V  = zeros(nx,ny,nz);
for k=1:nz
    A=imread(full_name,k);
    V(:,:,k)=A(:,:);
end
clear A
