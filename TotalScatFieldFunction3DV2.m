function EScatt= TotalScatFieldFunction3DV2(x_r1,y_r1,z_r1,k,sigma1, Gamma_0,R,R0,x1, y1, z1)
  
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%              Scattered Field  - Far Field                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    clear  EScatt 
  

     
           EScatt=0;
    for ii=1:numel(x_r1)
        %clear Drrj
        Drrj= sqrt( (x1-x_r1(ii)).^2+ (y1-y_r1(ii)).^2 + (z1-z_r1(ii)).^2);
        EScatt=EScatt - exp(1i*k*Drrj)./(1i*k*Drrj).*sigma1(ii);
    end
     
     
     


