function E_tot= TotalFieldFunction3D(x_r,y_r,z_r,W0,E0,k,sigma1, Gamma_0,R,phi,theta, vec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%              Total Field Intensity Function                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    clear x1 y1 EL1 EL I_E E_Scatt E EScattP EScattM EPtot EMtot

    
    %Laser Field

    %
    x1=R.*cos(phi).*sin(theta);
    y1=R.*sin(phi).*sin(theta);
    z1=R*cos(theta);
    %}
    [EL1, EL1x, EL1y] = LaserFieldFunction3D(x1,y1,z1,W0,E0,k);
       
    
    
    
    %Scattered Field
    if vec==0        
        
        %Scalar Light
        
        %Laser field for scalar polarization
        EL=EL1; 
        
        %Scattered field for scalar polarization
        EScatt=TotalScatFieldFunction3D(x_r,y_r,z_r,k,sigma1, Gamma_0,R,phi,theta);
        
        %Total field for scalar polarization
        E_tot=EL+1i*Gamma_0/2*EScatt;
        
    else
        
        %Vector Light        
        
        [Escattx, Escatty, Escattz]= TotalScatFieldFunction3DVec(x_r,y_r,z_r,k,sigma1, Gamma_0,R,phi,theta);
        
        %Total field in the x direction       
        Extot = EL1x + 1.*Gamma_0.*Escattx;
        
        %Total field in the y direction    
        Eytot = EL1y + 1.*Gamma_0.*Escatty;
    
        %Total field in the z direction    
        Eztot =  1.*Gamma_0.*Escattz;        
    
        %Coherent field
        %Total field
        E_tot = Extot + Eytot + Eztot; 

        

    
    end
        

           

  
