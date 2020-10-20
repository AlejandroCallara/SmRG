%% example of SmRG_mixtureModelFitting_multimix.m functioning. 
% get dimension
V  = SmRG_importData;
seed = SmRG_getSeed(V);
[nRow,nCol,nSli]=size(V);
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
xseed = round(seed (1,2));
yseed = round(seed (1,1));
zseed = round(seed (1,3));
initPos =[xseed,yseed,zseed];
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

Vx3 = V(xinit:xend,yinit:yend,zinit:zend);
%%
V_central_slice = Vx3(:,:,2);
vec_cs = V_central_slice(:);
[p_tot,a,K0_double,vB_double,mu_sk,rk] = SmRG_mixtureModelFitting_multmix(vec_cs');


figure, histogram(V_central_slice)
hold on
plot(50000*normpdf(1:max(V_central_slice(:)),K0_double,sqrt(vB_double)))
for irk = 1:length(rk)
    plot(50000*nbinpdf_mu(1:max(V_central_slice(:)),mu_sk(irk),1/rk(irk)))
end
drawnow
hold off









