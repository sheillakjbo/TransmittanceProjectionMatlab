function [Escattx, Escatty, Escattz]= TotalScatFieldFunction3DVecV2(x_r1,y_r1,z_r1,k0,sigma1, Gamma_0,R,R0,x1, y1, z1)
  
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%              Scattered Field  - Far Field                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    clear sumE EScatt 
  
    %syms phi theta
    N = numel(x_r1);
    %{
    x1=rr.*cos(phi);
    y1=rr.*sin(phi);
    z1=R0; % + R.*cos(theta);
  %}
  
    Escattx = 0; Escatty = 0; Escattz = 0;
    for ii = 1:N
        clear xoj Roj Got God
        xoj = x1-x_r1(ii); yoj = y1-y_r1(ii); zoj = z1-z_r1(ii);
        Roj = sqrt( xoj.^2+yoj.^2+zoj.^2 );

        God = 3/4*exp(1i*k0.*Roj)./(k0.*Roj).*(1 + 1i./(k0.*Roj) - 1./( (k0.*Roj).^2));
        
        %Green matrix: Second portion of the sum
        Got = 3/4*exp(1i*k0.*Roj)./(k0.*Roj).*(-1 - 3i./(k0.*Roj) + 3./((k0.*Roj).^2))./(Roj.^2);
          

        Escattx = Escattx - ((God+ Got.*xoj.^2).*sigma1(ii)...
            + Got.*xoj.*yoj.*sigma1(N+ii)...
            + Got.*xoj.*zoj*sigma1(2*N+ii));
        
        Escatty = Escatty - (Got.*xoj.*yoj.*sigma1(ii)...
            + (God + Got.*yoj.^2).*sigma1(N+ii)...
            + Got.*zoj.*yoj.*sigma1(2*N+ii));
       
        Escattz = Escattz - (Got.*xoj.*zoj.*sigma1(ii)...
            + Got.*yoj.*zoj.*sigma1(N+ii)...
            + (God + Got.*zoj.^2).*sigma1(2*N+ii));
    end
end
   
