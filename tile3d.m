function a=tile3d(b,barflag,minmax)

if (ndims(b)<3), a=b; return; end;
c=size(b);
if (nargin<2), 
    x=ceil(sqrt(c(3)));
    a1=zeros(c(1),c(2),x*x);
else
    x = barflag;
    a1 = zeros(c(1),c(2),size(b,3));
end;

a1(:,:,1:c(3))=b;

tmp=sprintf('a=[');
for mm=1:x,
  cc=0;
  for nn=1:x,
    cc=cc+sum(sum(a1(:,:,(mm-1)*x+nn)));
  end;
  %cc,
  if (cc==0), dothisrow=0; else, dothisrow=1; end;
  if (dothisrow),
    for nn=1:x,
      tmp=sprintf('%s a1(:,:,%d)''',tmp,(mm-1)*x+nn);
    end;
  end;
  tmp=sprintf('%s;',tmp);
end;
tmp=sprintf('%s];',tmp);
%disp(tmp);
eval(tmp);
clear a1

%a1=reshape(b(:,:,1:ceil(c(3)/2)),[c(1) c(2)*ceil(c(3)/2)]);
%a2=reshape(b(:,:,ceil(c(3)/2)+1:end),[c(1) c(2)*ceil(c(3)/2)]);
%a=reshape(b,[c(1) c(2)*c(3)]);

if nargout==0,
  if (nargin<3),
    %show(a)
    imagesc(a)
  else,
    %show(a,minmax)
    a=minmax(1)*(a<minmax(1))+minmax(2)*(a>minmax(2))+a.*((a>minmax(1))&(a<minmax(2)));
    imagesc(a)
  end;
  axis('off');
  title(sprintf('min/max = [%f %f]',min(min(a)),max(max(a))));
  if (nargin>2), colorbar, end;
  clear a
end;

