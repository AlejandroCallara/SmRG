function mask_out = SmRG_pruneMask(mask,currentSeed,v_init)
xi = v_init(1);
yi = v_init(2);
zi = v_init(3);
% ccmask = bwconncomp(mask,);
% for icc=1:length(ccmask.PixelIdxList)
%     [cc1, cc2, cc3]=(ind2sub(size(mask),ccmask.PixelIdxList{1,icc}));
%     SEMIcc= [cc1+xi-1 cc2+yi-1 cc3+zi-1];
%     if sum(ismember(SEMIcc,currentSeed,'rows'))
%         mask_out = false(size(mask));
%         mask_out(ccmask.PixelIdxList{1,icc})=true;
%     end
% %     mask_tmp(ccmask.PixelIdxList{1,icc})=true;
% %     SmRG_visSeed(mask_out,mask_tmp)
% end
% if ~exist('mask_out')
%     keyboard
% end
mask_inv = ~mask;
Nconn=26;
bw2 = imfill(mask_inv,currentSeed-v_init+1,Nconn);

mask_out = ~xor(mask,bw2);
