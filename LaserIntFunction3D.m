function I_EL= LaserIntFunction3D(R,R0,phi,rr,W0,E0,k,vec,beta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%              Laer Beam Intensity function                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    clear x1 y1 z1 EL I_EL x1 y1 ELP ELM

%
    x1=rr.*cos(phi);
    y1=rr.*sin(phi);
    z1=R0;
%}
    [EL1, EL1x, EL1y]=LaserFieldFunction3D(x1,y1,z1,W0,E0,k,beta);
    
    if vec==0
    
        EL=EL1;
    
    else
        
        EL=EL1x+EL1y;                 %Total Gaussian field   
        
    end

    I_EL=abs(EL).^2;

