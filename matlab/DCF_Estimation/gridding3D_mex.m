function  [kspCart]=gridding3D_mex(kspRad,kmap,kernel,opts)
% [kspCart]=gridding3D_mex(kspRad,kmap,kernel,opts)

% check if kspCart is complex
FlgR=0;
if(isreal(kspRad))
    FlgR = 1;
    kspRad = complex(kspRad,eps*ones(size(kspRad)));
end

% sort to improve rcn time when using cgrid with OMP
[tmp,idx1]=sort(kmap(:,3));
kmap=kmap(idx1,:); kspRad=kspRad(idx1);

kspCart = cGrid(single(kspRad), single(kmap), single(kernel), opts);
if(FlgR)
    kspCart = real(kspCart);
end

end
