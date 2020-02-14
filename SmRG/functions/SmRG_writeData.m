function SmRG_writeData(V,tipo,filename,pathname)
% SmRG_writeData: 
%           write 3D dataset to .tiff stack file.
%           Future updates will consider other formats. 
%      
% Syntax:
%           SmRG_writeData(V)
%           SmRG_writeData(V,tipo)
%
% Input:
%           V: 3D dataset double, uint8,uint16,uint32
%           tipo = one of 'uint8', 'uint16', 'uint32'


% dataset info
if nargin <1
    help SmRG_writeData
    return
end

if nargin <3
    if ~strfind(tipo,'uint')
        disp('image type is not one of uint8, uint16, uint32')
        disp('using uin8 as default')
        tipo = 'uint8';
    end
end
% image size
[nx,ny,nz]=size(V);

% output type
if ~exist('tipo')
    disp('image type not specified, using uint8 as default')
    tipo = 'uint8';
end

% export data
if nargin <4
    [filename,pathname]=uiputfile({'*.tif'});
    str_out=strcat(pathname,filename);
    fprintf('saving .tiff stack to \n')
    disp(str_out)
else 
    str_out=strcat(pathname,filename);
end



classe = class(V);

if ~strcmp(classe,'double') 
    V = double(V);
end
vtmp = (V-min(V(:)))/(max(V(:))-min(V(:)));
vtmp = single(vtmp);
c = strsplit(tipo,'uint');
switch tipo
    case 'uint8'
        N = str2num(c{1,2});
        tmp_out = uint8((2^N-1)*(vtmp));
    case 'uint16'
        N = str2num(c{1,2});
        tmp_out = uint16((2^N-1)*(vtmp));
    case 'uint32'
        N = str2num(c{1,2});
        tmp_out = uint32((2^N-1)*(vtmp));
end

for K=1:nz
   imwrite(tmp_out(:, :, K), str_out, 'WriteMode', 'append');
end
disp('done')

