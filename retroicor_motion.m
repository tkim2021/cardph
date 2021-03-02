function results = retroicor_motion(order, dims, img, mask, cardph, respph, motion, physio_s, display)
M=order;

% modified from PESTICA v2.0 matlab_retroicor.m
% https://www.nitrc.org/frs/?group_id=361

xdim = dims.x; ydim = dims.y; zdim = dims.z; tdim = dims.t;

mmx = motion(1,:)';
mmy = motion(2,:)';
mmz = motion(3,:)';
degx = motion(4,:)'; % deg (rad*180/pi) 
degy = motion(5,:)';
degz = motion(6,:)';
% mcflirt returns rad, but HCP convert to deg (x 180/pi)

retima=zeros(xdim,ydim,zdim,tdim);
im_ca=zeros(xdim,ydim,zdim,M);
im_cb=zeros(xdim,ydim,zdim,M);
im_ra=zeros(xdim,ydim,zdim,M);
im_rb=zeros(xdim,ydim,zdim,M);
tim_ca=zeros(xdim,ydim,zdim,M);
tim_cb=zeros(xdim,ydim,zdim,M);
tim_ra=zeros(xdim,ydim,zdim,M);
tim_rb=zeros(xdim,ydim,zdim,M);

%dof=double(tdim-(M*2*2+1));
dof=double(tdim-(M*2*2+1+6)); % The df(Residual) is the sample size minus the number of parameters being estimated,
disp('Running RETROICOR...');

if strcmp(physio_s, 'c') || strcmp(physio_s, 'r') || strcmp(physio_s, 'm')
    pim = zeros(xdim,ydim,zdim,M*2+1);
    pmotionim = zeros(xdim,ydim,zdim,6);
end
if strcmp(physio_s, 'cr') 
    pim = zeros(xdim,ydim,zdim,M*2*2+1);
    pmotionim = zeros(xdim,ydim,zdim,6);
end
    
for z=1:zdim
  if strcmp(physio_s, 'c')
    A=[sin((1:M)'*cardph(z,:))' cos((1:M)'*cardph(z,:))' mmx mmy mmz degx degy degz ones(tdim,1)];
  end
  if strcmp(physio_s, 'r')
    A=[sin((1:M)'*respph(z,:))' cos((1:M)'*respph(z,:))' mmx mmy mmz degx degy degz ones(tdim,1)];
  end
  if strcmp(physio_s, 'cr')
    A=[sin((1:M)'*cardph(z,:))' cos((1:M)'*cardph(z,:))' sin((1:M)'*respph(z,:))' cos((1:M)'*respph(z,:))' mmx mmy mmz degx degy degz ones(tdim,1)];
  end
  if strcmp(physio_s, 'm')
     A=[mmx mmy mmz degx degy degz ones(tdim,1)];
  end
                                                                                      
  % QR-decomposition
  [Q,R] = qr(A,0);
  
  mask_z = mask(:,:,z);
  im_z = squeeze(img(:,:,z,:));
  mask1D = reshape(mask_z, xdim*ydim,1);
  im1Dt = reshape(im_z,xdim*ydim,tdim);
  im1D = im1Dt(logical(mask1D),:);

  tmpim_ca = zeros(size(im1D,1),M); %zeros(length(logical(mask1D)),M);
  tmpim_cb = tmpim_ca; 
  tmpim_ra = tmpim_ca;
  tmpim_rb = tmpim_ca; 
  tmptim_ca = tmpim_ca; 
  tmptim_cb = tmpim_ca; 
  tmptim_ra = tmpim_ca; 
  tmptim_rb = tmpim_ca; 
  tmpretima = zeros(size(im1D,1),tdim);
  fittedim = tmpretima;
  tmp_pim = zeros(size(im1D,1),size(pim,4));
  
  tmp_motionim = zeros(size(im1D,1),6);
  
  p2 = zeros(size(im1D,1),1);
  
  for ii=1:size(im1D,1)
      vox = squeeze(im1D(ii,:)'); 
%   for y=1:ydim
%     for x=1:xdim
%       vox=squeeze(img(x,y,z,:));
%       if (mask(x,y,z)==0)
%         retima(x,y,z,:)=vox;
%         continue;
%       end
      % matrix division to get amplitudes of fits
      p = R\(Q'*vox); % the same operation with p = Q*R\vox;
%      p=regress(vox,A);
      if strcmp(physio_s, 'c') || strcmp(physio_s, 'r')
          tmp_pim(ii,:) = p([1:2*M end]);
          tmp_motionim(ii,:)= p(M*2+1:M*2+6);
      end
      if strcmp(physio_s, 'cr')
          tmp_pim(ii,:) = p([1:2*M*2 end]);
          tmp_motionim(ii,:)= p(2*M*2+1:2*M*2+6);
      end
      
%       im_ca(x,y,z,:)=p(1:M);
%       im_cb(x,y,z,:)=p(M+1:M*2);
      if strcmp(physio_s, 'c') || strcmp(physio_s, 'r')
        tmpim_ca(ii,:)=p(1:M);
        tmpim_cb(ii,:)=p(M+1:M*2);
      end
      
      if strcmp(physio_s, 'cr')
        tmpim_ca(ii,:)=p(1:M);
        tmpim_cb(ii,:)=p(M+1:M*2);
        tmpim_ra(ii,:)=p(M*2+1:M*3);
        tmpim_rb(ii,:)=p(M*3+1:M*4);
      end
      
      % residuals
      res=vox-A*p; % zeros mean
      % mean square error
      mse=res'*res./(dof+1);
      Rinv=pinv(R);
      % error covariance matrix
      covb=(Rinv*Rinv')*mse;
      % get standardized error (variance that is not explained by fit to design matrix)
      stand_err=sqrt(diag(covb));
      % test statistics
%       tim_ca(x,y,z,:)=p(1:M)./stand_err(1:M);
%       tim_cb(x,y,z,:)=p(M+1:M*2)./stand_err(M+1:M*2);
      if strcmp(physio_s, 'c') || strcmp(physio_s, 'r')
         tmptim_ca(ii,:)=p(1:M)./stand_err(1:M);
         tmptim_cb(ii,:)=p(M+1:M*2)./stand_err(M+1:M*2);
      end
      if strcmp(physio_s, 'cr') 
         tmptim_ca(ii,:)=p(1:M)./stand_err(1:M);
         tmptim_cb(ii,:)=p(M+1:M*2)./stand_err(M+1:M*2);
         tmptim_ra(ii,:)=p(M*2+1:M*3)./stand_err(M*2+1:M*3);
         tmptim_rb(ii,:)=p(M*3+1:M*4)./stand_err(M*3+1:M*4);
      end
      
      % residuals
      % don't remove the mean from voxel timeseries
%       p(end)=0;     % p(end) is mean of img. this does not subtract the mean for res
      %retima(x,y,z,:)=vox-A*p;
      %tmpretima(ii,:)= res;   % the vox-A*p has zero mean because p(end) includes mean intensity 
      %tmpretima(ii,:)=vox-A(:,1:end-1)*p(1:end-1); % keep signal intensity because the last comp has mean intensity
      tmpretima(ii,:)=vox-A*p; % remove base intensity (compare with the previous)
      
      %fittedim(ii,:) = A*p;
      % should be remove motion compoent !!!
      if strcmp(physio_s, 'm')
          fittedim(ii,:) = A(:,1:6)*p(1:6);
      else
          fittedim(ii,:) = A(:,[1:4 11])*p([1:4 11]);
      end
      
      %p2(ii) = regress(vox,A*p);
  %  end
  end
  
  if strcmp(physio_s, 'cr') 
     tmp = zeros(xdim*ydim,M);
     tmp(logical(mask1D),:) = tmpim_ca;
     im_ca(:,:,z,:) = reshape(tmp,xdim,ydim,M);
 
     tmp(logical(mask1D),:) = tmpim_cb;
     im_cb(:,:,z,:) = reshape(tmp,xdim,ydim,M);
     
     tmp(logical(mask1D),:) = tmpim_ra;
     im_ra(:,:,z,:) = reshape(tmp,xdim,ydim,M);

     tmp(logical(mask1D),:) = tmpim_rb;
     im_rb(:,:,z,:) = reshape(tmp,xdim,ydim,M);
  
     tmp(logical(mask1D),:) = tmptim_ca;
     tim_ca(:,:,z,:) = reshape(tmp,xdim,ydim,M);
     tmp(logical(mask1D),:) = tmptim_cb;
     tim_cb(:,:,z,:) = reshape(tmp,xdim,ydim,M);
     
     tmp(logical(mask1D),:) = tmptim_ra;
     tim_ra(:,:,z,:) = reshape(tmp,xdim,ydim,M);
     tmp(logical(mask1D),:) = tmptim_rb;
     tim_rb(:,:,z,:) = reshape(tmp,xdim,ydim,M);
  
     tmp2 = zeros(xdim*ydim,tdim);
     tmp2(logical(mask1D),:) = tmpretima;
     retima(:,:,z,:) = reshape(tmp2,xdim,ydim,tdim);
  
     tmp2 = zeros(xdim*ydim,tdim);
     tmp2(logical(mask1D),:) = fittedim;
     fitted_im(:,:,z,:) = reshape(tmp2,xdim,ydim,tdim); 
     % in case of 'cr', fitted_im is useless because of cr-combined map
  
%      tmp3 = zeros(xdim*ydim,1);
%      tmp3(logical(mask1D)) = p2;
%      one_p(:,:,z) = reshape(tmp3,xdim,ydim); % similar as pestica
     
     tmp4 = zeros(xdim*ydim,size(pim,4));
     tmp4(logical(mask1D),:) = tmp_pim;
     pim(:,:,z,:) = reshape(tmp4,xdim,ydim,size(pim,4));
  end
  
  if strcmp(physio_s, 'c') || strcmp(physio_s, 'r') || strcmp(physio_s, 'm')
     tmp = zeros(xdim*ydim,M);
     tmp(logical(mask1D),:) = tmpim_ca;
     im_ca(:,:,z,:) = reshape(tmp,xdim,ydim,M);
 
     tmp(logical(mask1D),:) = tmpim_cb;
     im_cb(:,:,z,:) = reshape(tmp,xdim,ydim,M);
  
     tmp(logical(mask1D),:) = tmptim_ca;
     tim_ca(:,:,z,:) = reshape(tmp,xdim,ydim,M);
     tmp(logical(mask1D),:) = tmptim_cb;
     tim_cb(:,:,z,:) = reshape(tmp,xdim,ydim,M);
  
     tmp2 = zeros(xdim*ydim,tdim);
     tmp2(logical(mask1D),:) = tmpretima;
     retima(:,:,z,:) = reshape(tmp2,xdim,ydim,tdim);
  
     tmp2 = zeros(xdim*ydim,tdim);
     tmp2(logical(mask1D),:) = fittedim;
     fitted_im(:,:,z,:) = reshape(tmp2,xdim,ydim,tdim);
  
%      tmp3 = zeros(xdim*ydim,1);
%      tmp3(logical(mask1D)) = p2;
%      one_p(:,:,z) = reshape(tmp3,xdim,ydim); % similar as pestica
     
     tmp4 = zeros(xdim*ydim,size(pim,4));
     tmp4(logical(mask1D),:) = tmp_pim;
     pim(:,:,z,:) = reshape(tmp4,xdim,ydim,size(pim,4));
     
     tmp5 = zeros(xdim*ydim,6);
     tmp5(logical(mask1D),:) = tmp_motionim;
     pmotionim(:,:,z,:) = reshape(tmp5,xdim,ydim,6);
  end
  
end   % for zdim
disp('Completed RETROICOR');
clear tmp tmp1 tmp2 tmp3 tmp4

% correct for induced negative values
%retima=retima-min(retima(:));
%img_varnom=img_varnom-min(img_varnom(:));

% sum squares of coupling coefficients
im_card=sqrt(sum(im_ca.^2+im_cb.^2,4));
im_resp=sqrt(sum(im_ra.^2+im_rb.^2,4));

% sum squares of t-statistics
tim_card=sqrt(sum(tim_ca.^2+tim_cb.^2,4));
tim_resp=sqrt(sum(tim_ra.^2+tim_rb.^2,4));

% re-introduce the image variance
%retima=retima.*repmat(variance,[1 1 1 tdim]);
% img_varnorm=img_varnom.*repmat(variance, [1 1 1 tdim]);

if display == 1
%    figure, imagesc(tile3d(im_card),[0 500])
%    title('all slice phase for card')
%    figure, imagesc(tile3d(im_resp),[0 500])
%    title('all slice phase for resp')

%end
    figure, imagesc(tile3d(tim_card(:,:,1:2:end)),[0 20])
    title('all slice phase for t-stat for card')
%    figure, imagesc(tile3d(tim_resp),[0 20])
%    title('all slice phase for t-stat for resp')
end

field1 = 'im_ca'; value1 = im_ca;
field2 = 'im_cb'; value2 = im_cb;
field3 = 'im_ra'; value3 = im_ra;  % resp
field4 = 'im_rb'; value4 = im_rb;
field5 = 'tim_ca'; value5 = tim_ca;
field6 = 'tim_cb'; value6 = tim_cb;
field7 = 'tim_ra'; value7 = tim_ra;
field8 = 'tim_rb'; value8 = tim_rb;
field9 = 'im_card'; value9 = im_card;
field10 = 'im_resp'; value10 = im_resp;
field11 = 'tim_card'; value11 = tim_card;
field12 = 'tim_resp'; value12 = tim_resp;
field13 = 'retima'; value13 = retima;
field14 = 'p_im'; value14 = pim;
field15 = 'p_motion_im'; value15 = pmotionim;
field16 = 'fitted_im'; value16 = fitted_im;
results = struct(field1,value1, field2,value2, field3,value3, field4,value4, ...
    field5,value5, field6,value6, field7,value7, field8,value8, ...
    field9,value9, field10,value10, field11,value11, field12,value12, ...
    field13,value13, field14,value14, field15,value15, field16,value16);

% original = fitted_im + retima + p_motion_im
