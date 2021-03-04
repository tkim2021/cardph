
function find_cardph_hilbert_sl(subjectID, fMRIfile, root_dir, orig_dir)

order = 2; 
M = order;
MB = 8; 
TR = 0.8;

data_dir=strcat(root_dir,filesep,subjectID,filesep,fMRIfile, filesep);

fsl_path = '/usr/local/fsl';
setenv('FSL_DIR', fsl_path)
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); 

imagefile = strcat(orig_dir,filesep,subjectID,filesep,fMRIfile,filesep,fMRIfile,'_orig.nii.gz');
xfmsfile = strcat(root_dir,filesep,subjectID,filesep,'MNINonLinear/xfms/standard2',fMRIfile,'.nii.gz');


outputfile = strcat(data_dir,filesep,'card2',fMRIfile,'_norm.nii.gz');
if ~exist(outputfile,'file')
    system(sprintf('/usr/local/fsl/bin/applywarp -i Avg_cardcouplingICA2_norm.nii.gz -r %s -w %s -o %s',imagefile,xfmsfile,outputfile));
end
 
nii = load_untouch_nii(outputfile); 
cardimgall = nii.img;
checkmask = zeros(size(cardimgall));
checkmask(cardimgall > 3.29)=1; % 1.960(95), 2.326(98), 2.576(99), 3.29(0.001), 3.89(0.0001), 4.41 (.00001)

fileID = fopen(strcat(data_dir,filesep,'Movement_Regressors.txt'));
Am = fscanf(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f',[12 Inf]);
fclose(fileID);
motion6 = Am(1:6,:);
% if partial_ts == 1
%     motion6 = motion6(:,1:pts);
% end

figure
subplot(2,1,1), plot(motion6(1:3,:)') % xyz translation
title('trans: blue:x,red:y,ornage:z'); grid on
subplot(2,1,2), plot(motion6(4:6,:)') % xyz rotation, deg (not rad)
title('rot: blue:x,red:y,ornage:z') 
grid on
% mcflirt returns rad, but HCP converts to deg (x 180/pi)
% saveas(gcf,strcat(data_dir,'motion'),'jpeg')

% read img
nii = load_untouch_nii(imagefile);
origimg = double(nii.img);

dims.x = size(origimg,1); dims.y = size(origimg,2); dims.z = size(origimg,3); dims.t = size(origimg,4); 
origimg1 = reshape(origimg,[],dims.t);

% brain mask
origimgmask = zeros(dims.x, dims.y, dims.z);
origimgmask(mean(origimg,4) > mean(origimg(:))) = 1;
% fill holes in the mask. 2D is necessay for imfill
BW = imbinarize(origimgmask);
for ii=1:dims.z
    origimgmask(:,:,ii) = imfill(BW(:,:,ii),'holes');
end


load('slicetiming.txt');
B = sort(unique(slicetiming));


%% MB

final_angph_array = zeros(length(B), dims.t);
s_origallam_array = zeros(length(B), dims.t);


slicemask = checkmask;
slicemask(:,:,1:9) = 0;

for jj=1:length(B) 
    % select slices for separate PCA
    sl = find(slicetiming==B(jj));
    maskimg = zeros(dims.x, dims.y, dims.z);
    maskimg(:,:,sl) = slicemask(:,:,sl); %slicemask = checkmask
    target_pixels = origimg1(logical(maskimg),:);
    
    tmp_tc_origall = double(target_pixels');
    tmp_tc_origall = detrend(tmp_tc_origall);
    tmp_tc_origallam = zeros(size(tmp_tc_origall));
    for ii=1:size(tmp_tc_origall,2)
        pp = regress(tmp_tc_origall(:,ii), motion6');
        tmp_tc_origallam(:,ii) = tmp_tc_origall(:,ii) - (pp'*motion6)';
    end
    
    [final_angph_origallam, s_origallam,~,~] = pca_finalangph2(tmp_tc_origallam, dims);
   
    final_angph_array(jj,:) = final_angph_origallam;

%     s_origallam_array(jj,:) = s_origallam;
end

for jj=1:length(B)
    sl = find(slicetiming==B(jj));
    part_img = double(origimg(:,:,sl,:));
    part_mask = origimgmask(:,:,sl);
    
    part_slicemask = slicemask(:,:,sl); %
    part_respph = ones(MB,dims.t);
    part_dims = dims; part_dims.z = MB;

    part_ph = repmat(final_angph_array(jj,:),MB,1);
    sl_origallam(jj) = retroicor_motion(order, part_dims, part_img, part_mask, part_ph, part_respph, motion6, 'c', 0);
    
end

%% slices to whole brain
MB_sl_origallam.tim_card = zeros(dims.x,dims.y,dims.z);
MB_sl_origallam.im_card = zeros(dims.x,dims.y,dims.z);
for jj=1:length(B)
    sl = find(slicetiming==B(jj));
    MB_sl_origallam.tim_card(:,:,sl) = sl_origallam(jj).tim_card;
    %MB_sl_origallam.im_card(:,:,sl) = sl_origallam(jj).im_card;
end

figure, imagesc(tile3d(MB_sl_origallam.tim_card)

% ========================================
    

%% save data
nii.hdr.dime.datatype = 16;
nii.hdr.dime.bitpix = 32;
nii.hdr.dime.dim(1) = 3;
nii.hdr.dime.dim(5) = 1;

nii.img = MB_sl_origallam.tim_card;
save_untouch_nii(nii,strcat(data_dir,'card_hilbert_sl_origallam_', fMRIfile,'.nii.gz'));

%% orig to MNI
% orig_fall = strcat(data_dir,'card_hilbert_origallam_', fMRIfile,'.nii.gz');
% orig_fslall = strcat(data_dir,'card_hilbert_sl_origallam_', fMRIfile,'.nii.gz');
% %physio_f = strcat(data_dir,'card_physiocard_', fMRIfile,'.nii.gz');
% MNI_fall = strcat(data_dir,'card_hilbert_origallam_', fMRIfile,'_MNI.nii.gz');
% MNI_fslall = strcat(data_dir,'card_hilbert_sl_origallam_', fMRIfile,'_MNI.nii.gz');
% %MNI_f = strcat(data_dir,'card_physiocard_', fMRIfile,'_MNI.nii.gz');
% rfile = strcat(root_dir,subjectID,filesep,'MNINonLinear/T1w_restore.2.nii.gz');
% xfmsfile2 = strcat(root_dir,subjectID,filesep,'MNINonLinear/xfms/',fMRIfile,'2standard.nii.gz');
% system(sprintf('/usr/local/fsl/bin/applywarp -i %s -r %s -w %s -o %s', orig_fall, rfile, xfmsfile2, MNI_fall));
% system(sprintf('/usr/local/fsl/bin/applywarp -i %s -r %s -w %s -o %s', orig_fslall, rfile, xfmsfile2, MNI_fslall));
