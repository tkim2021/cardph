
function find_cardph_hilbert(subjectID, fMRIfile, root_dir, orig_dir)

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

%% 3D->2D
origmotion6 = motion6;
target_pixels = origimg1(logical(checkmask),:);
tmp_tc_origall = double(target_pixels');

tmp_tc_origall = detrend(tmp_tc_origall);

tmp_tc_origallam = zeros(size(tmp_tc_origall));
for ii=1:size(tmp_tc_origall,2)
    pp = regress(tmp_tc_origall(:,ii), origmotion6');
    tmp_tc_origallam(:,ii) = tmp_tc_origall(:,ii) - (pp'*origmotion6)';
end

% find ph for allsl
[final_angph_origallam_all, s_origallam_all, cardph_origallam_all, respph_oriallgam_all] = pca_finalangph2(tmp_tc_origallam, dims);
H= figure; plot(s_origallam_all,'o-'); title('orig all vessel motion regree'); saveas(H,strcat(data_dir,'g_orig_allam'),'jpeg')

load('slicetiming.txt');
B = sort(unique(slicetiming));

cardph_origallam = zeros(dims.z, dims.t);
for jj=1:length(B) 
    sl = find(slicetiming==B(jj));
    cardph_origallam(sl,:) = repmat(final_angph_origallam_all+ 2*pi*B(jj)/1000,1,MB)';
end
simple_origallam_ALLsli = retroicor_motion(order, dims, double(origimg), origimgmask, cardph_origallam, respph_oriallgam_all, motion6, 'c', 0);
figure, imagesc(tile3d(simple_origallam_ALLsli.tim_card),[0 10]), title('simple origallam - old way')



%% save data
nii.hdr.dime.datatype = 16;
nii.hdr.dime.bitpix = 32;
nii.hdr.dime.dim(1) = 3;
nii.hdr.dime.dim(5) = 1;

nii.img = simple_origallam_ALLsli;
save_untouch_nii(nii,strcat(data_dir,'card_hilbert_origallam_', fMRIfile,'.nii.gz'));

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
