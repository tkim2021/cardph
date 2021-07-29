function [final_angph, s, cardph, respph] = pca_finalangph2(tmp_tc, dims)
% COPYRIGHT 2021 Tae Kim (tak19@pitt.edu)

%[coeff, score, latent, ~, explained, mu] = pca(tmp_tc);
[~, score, ~, ~, ~, ~] = pca(tmp_tc);
s1 = score(:,2);
% s = rescale(s1,-1,1); % not working for 2016b
s = math_scale_values(s1,min(s1),max(s1),-1,1);

ang = wrapTo2Pi(angle(hilbert(s)));
%uang1 = unwrap(ang);   % unwarp to align the first phase to 0
%simpleang = uang1-uang1(1);
%final_angph = wrapTo2Pi(simpleang);

final_angph = ang;

cardph = zeros(dims.z, dims.t); respph = cardph;
for ii=1:dims.z
    cardph(ii,:) = final_angph;
    respph(ii,:) = ones(1,dims.t);
end

