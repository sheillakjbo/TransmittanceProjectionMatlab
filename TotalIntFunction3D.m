function I_E= TotalIntFunction3D(x_r,y_r,z_r,W0,E0,k,sigma1, Gamma_0,R, R0,phi,rr, vec,beta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%              Total Field Intensity Function                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    clear x1 y1 EL1 EL I_E E_Scatt E EScattP EScattM EPtot EMtot

    
    %Laser Field

    %
    x1=rr.*cos(phi);
    y1=rr.*sin(phi);
    z1=R0;
    %}
    [EL1, EL1x, EL1y] = LaserFieldFunction3D(x1,y1,z1,W0,E0,k,beta);
       
    
    
    
    %Scattered Field
    if vec==0        
        
        %Scalar Light
        
        %Laser field for scalar polarization
        EL=EL1; 
        
        %Scattered field for scalar polarization
        EScatt=TotalScatFieldFunction3D(x_r,y_r,z_r,k,sigma1, Gamma_0,R,R0,phi,rr);
        
        %Total field for scalar polarization
        E=EL+1i*Gamma_0/2*EScatt;
        
        %I_E=abs(conj(EL).*E).^2;  
        I_E=abs(E).^2;  
        
    else
        
        %Vector Light        
        
        [Escattx, Escatty, Escattz]= TotalScatFieldFunction3DVec(x_r,y_r,z_r,k,sigma1, Gamma_0,R, R0,phi,rr);
        
        %Total field in the x direction       
        Extot = EL1x + Gamma_0.*Escattx;
        
        %Total field in the y direction    
        Eytot = EL1y + Gamma_0.*Escatty;
    
        %Total field in the z direction    
        Eztot =  1*Gamma_0.*Escattz;        
    
        %Coherent field
        %Total field
        E = Extot + Eytot + Eztot; 

        %I_E=abs(conj(EL1x).*Extot +conj(EL1y).*Eytot).^2;  
        I_E=abs(E).^2;
    
    end
        
 
         

  
