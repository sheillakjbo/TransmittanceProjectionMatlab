function [EL, ELx, ELy]= LaserFieldFunction3D(x,y,z,W0,E0,k,beta)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                Laser Beam Field Function                        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    clear alpha EL
    %beta = 1;
    alpha=1+2*1i*beta*z./(k*W0^2);
    
    EL=E0*exp(1i*k*z)./(alpha).*exp(-1*(x.^2+y.^2)./(alpha.*W0^2));
    ELx = EL./sqrt(2);     %Gaussian beam in the x direction
    ELy = 1i*EL./sqrt(2);     %Gaussian beam in the y direction    
    %
    
    %{
    E0 = 1; %sqrt(2*I0/c/epsilon0);

    % Incident beam (Gaussian)
    w0 = 1e-6;
    zr = pi*w0^2/lambda; % Rayleigh length
    waist = @(z) w0*sqrt(1+(z/zr).^2);
    Rcurv = @(z) z.*(1+(zr./z).^2);
    Gouyphase = @(z) atan2(z,zr);
    
    r2 = x.^2+y.^2;
    EL = E0*w0./waist(z).*exp( -r2./waist(z).^2 + ...
    1i*k*z + 1i*k*r2/2./Rcurv(z) - 1i*Gouyphase(z) );
    ELx = EL./sqrt(2);     %Gaussian beam in the x direction
    ELy = EL./sqrt(2);     %Gaussian beam in the y direction  
    %}