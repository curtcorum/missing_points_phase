function [kspRad] = regrid3D_mex(kspCart,kmap,kernel,opts)

% check if kspCart is complex
flgR=0;
if(isreal(kspCart))
    flgR = 1;
    kspCart = complex(kspCart,eps*ones(size(kspCart)));
end

[kspRad]=cInvGridNew(single(kspCart),single(kmap),single(kernel),opts);
if(flgR)
    kspRad = real(kspRad);
end
end