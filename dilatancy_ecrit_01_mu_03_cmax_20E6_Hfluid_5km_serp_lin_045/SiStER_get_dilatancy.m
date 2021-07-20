function [psi]=SiStER_get_dilatancy(im,ep,MAT)
% [fric]=SiStER_get_friction(im,ep,MAT)
% compute friction on markers based on plastic strain, like cohesion
% J.-A. Olive 4/2017


psi = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    psimax=MAT(types(i)).psi;

    epscrit=MAT(types(i)).ecrit;

    % get cohesion
    psi(imInd)=max(psimax+(0-psimax).*ep(imInd)./epscrit,0);
end

return


