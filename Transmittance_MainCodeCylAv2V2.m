tic

clear all

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%




                %%%    General Parameters      %%%
beta = 1;          
zz = 2;
vec = 1;                         %vec~=1 => Scalat light If vec==1 => Vector light         
ro_lambda3 = 0.1;
%N = 1000;                        %Number of particles
steps = 11;
%k0 = 1;                         %wave number
lamb = 780*1e-9;
k0 = 2*pi/lamb;
k = k0;
%lamb = 2*pi/k;
R = 4*lamb;
x0 = 0;                          %x size
y0 = 0;
z0 = zz*1e-6; % 5*lamb; %np.linspace(2,10,steps).reshape((steps,1))                    %z size
N = round(ro_lambda3/lamb^3*(2*pi)^(1/2)*z0*pi*R^2);
Gamma_0 = 2*pi*30.5e6;                     %Decay rate
Gamma_1 = Gamma_0/2;             %Decay rate
Domega = transpose(linspace(-10,10,steps))*Gamma_0;       % 1*np.ones((steps,1))*6%             %Single value of Domega (for tests)
W0 = 3*lamb;
E0 = 1;
Rlz = 1;
rho = N/(pi*R^2*z0);
R0 = 50*R ; % 10*R+ (k*(W0^2))/2 ;  %Distance of the sensors for far field approximation
rm = 8;
                %%%   Fix field parameters       %%%
    %Optical thickness


%%%%%%%%%%%%%%%%%%%%% Adjusted optical thickness %%%%%%%%%%%%%%%%%%%%%%
if vec==1
     %b0 = 6*N./((R*k)^2);
     b0 = 6*N*erf(1/sqrt(2))/((k*R)^2*sqrt(2*pi));%Normal z-axis distribution exact
     delt=2*Domega./Gamma_0;
     b=b0./(1+delt.^2);
     %b=2.*b0./(16.*Domega.^2+1);   %Optical thickness of vector light

else
    %b0 = 4*N./((R*k)^2); %Uniform distribution
    b0 = 4*N*erf(1/sqrt(2))/((k*R)^2*sqrt(2*pi)); %Normal z-axis distribution exact
    delt=2*Domega./Gamma_0;
    b = b0./(1+delt.^2);
    %b=b0./(4*Domega.^2+1);      %Optical thickness of scalar light

end


if vec==1
    NN=3*N;
else
    NN=N;
end

c=rho/k^3;
cc = rho*lamb^3;
stg=string(c);
dens=strrep(stg,'.','');


zr = pi*W0^2/lamb; % Rayleigh length
waist = @(z) W0*sqrt(1+(z/zr).^2);
R1 = 5*waist(R0);

%Funnction of the laser beam intensity for theta%
%I_EL=@(phi, rr) LaserIntFunction3D(R1, R0,phi,rr,W0,E0,k,vec,beta); 


%Integral of the laser beam intensity in theta%
%IL_tot=integral2(I_EL,0,2*pi,0,R1);


%%%%%%%%%%%%%%%%% INCOHERENT fIELD %%%%%%%%%%%%%%%%%%%%%%%%%
I_tot = zeros(steps,1);
for rr=1:Rlz
    
    [x,y,z, DeltaX, DeltaY, DeltaZ, D] = Data(N,z0,R,rm);
    

    for ii=1:steps % Domega(i)=> Detuning (detuning loop)
        Domega1 = Domega(ii);
        sigma1 = Sigma(Gamma_0,N,k0,Domega1,W0,E0,x, y, z, DeltaX, DeltaY, DeltaZ, D,vec,beta);
        %Function of the total intensity for theta and Domega%                 

        I_E=@(phi,rr) TotalIntFunction3D(x,y,z,W0,E0,k,sigma1, Gamma_0,R1,R0,phi,rr, vec,beta);
        
        
        
        HalfRange = R1; 
        nn = 51;
        dl = 2*HalfRange/(nn-1);
        Xobsvec = linspace(-1,1,nn)*HalfRange;
        Yobsvec = Xobsvec;
        [Xobs,Yobs] = meshgrid(Xobsvec,Yobsvec);
        Zobs = R0*ones(size(Xobs));
        %NRange = length(Xobs);
        
        El = LaserFieldFunction3D(Xobs,Yobs,Zobs,W0,E0,k,beta);
        

      
        %Scalar
        %{        
        Etot = El+1i*Gamma_0/pi*TotalScatFieldFunction3DV2(x,y,z,k,sigma1, Gamma_0,R1,R0,Xobs, Yobs, Zobs);

        I_tot(ii)= I_tot(ii) + (abs(sum(sum(conj(El).*Etot*dl^2)))/(pi/2*W0^2)).^2/Rlz;
        %}


        %
        %vector
        [Escattx, Escatty, Escattz]= TotalScatFieldFunction3DVecV2(x,y,z,k,sigma1, Gamma_0,R1,R0,Xobs, Yobs, Zobs);
        Etotx = El/sqrt(2)+1i/2*Gamma_0*Escattx;
        Etoty = 1i*El/sqrt(2)+1i/2*Gamma_0*Escatty;
        %Etotz = 1i*Gamma_0/2*Escat3z;        
        
        I_tot(ii)= I_tot(ii) + (abs(sum(sum((conj(El/sqrt(2)).*Etotx + conj(1i*El/sqrt(2)).*Etoty)*dl^2) ))/(pi/2*W0^2)).^2/Rlz;
        %}
        
    end

end

    


  

%%%   Transmittance for incoherent field   %%
T_c=I_tot;%./IL_tot;   










%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%Paht where the figures will be saved
if vec==1
    
    pth='/home/sheilla/Documents/MATLAB/Transmittance_figs/bVary/Vector';
    
else
    
   pth='/home/sheilla/Documents/MATLAB/Transmittance_figs/bVary/Scalar'; 
    
end    


            %Normalized density


%{
figure
fig1=plot(b,exp(-b),'k');
hold on
plot(b((steps-1)/2+1:steps),T_c((steps-1)/2+1:steps),'r-o')
hold on
plot(b(1:(steps-1)/2+1),T_c(1:(steps-1)/2+1),'b-o')
%plot(b((steps-1)/2+1:steps),T((steps-1)/2+1:steps),'b-*')
xlabel('b');
ylabel('Transmission');
title(['N=',num2str(N), '; b_0=',num2str(b0), '; \rho/k²=',num2str(c)])
legend('Beer-Lambert Law','Delta omega<1','Delta omega >1')

figure
fig2=plot(b,exp(-b),'k');
set(gca,'yscale','log','xscale','log')
hold on
plot(b((steps-1)/2+1:steps),T_c((steps-1)/2+1:steps),'r-o')
set(gca,'yscale','log','xscale','log')
hold on
plot(b(1:(steps-1)/2+1),T_c(1:(steps-1)/2+1),'b-o')
%plot(b((steps-1)/2+1:steps),T((steps-1)/2+1:steps),'b-*')
set(gca,'yscale','log','xscale','log')
xlabel('b');
ylabel('Transmission');
title(['N=',num2str(N), '; b_0=',num2str(b0),  '; \rho/k²=',num2str(c)])
%legend('Coherent: |<E>|^2','Incoherent: <|E|^2>')
%}
figure
fig3=plot(Domega,exp(-b),'k');
hold on
plot(Domega,T_c,'r-o');
%hold on
%plot(Domega,T,'b-*')
xlabel('\delta/\Gamma_0');
ylabel('Transmission');
title(['N=',num2str(N), '; b_0=',num2str(b0),  '; \rho/k²=',num2str(c)])
legend('B-L','Coherent: |<E>|^2','Incoherent: <|E|^2>')

figure
fig4=plot(Domega,T_c,'r-o');
set(gca,'yscale','log')
%hold on
%plot(Domega,T,'b-*')
%set(gca,'yscale','log')
xlabel('\delta/\Gamma_0');
ylabel('Transmission');
title(['N=',num2str(N), '; b_0=',num2str(b0),  '; \rho/k²=',num2str(c)])
%legend('Coherent: |<E>|^2','Incoherent: <|E|^2>')

%Saving data
Data1= [T_c T_c b Domega];
%save(['r_min' num2str(rm) '.mat'],'Data1')
%save b048.txt T T_c b Domega -ascii

%Loading saved data
%filename = 'rlz100.txt';
%[A6,delimiterOut]=importdata(filename);


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Saving figures       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
filename1 = ['N',num2str(N),'c',num2str(dens)];
saveas(fig1,fullfile(pth, filename1), 'fig');
saveas(fig1,fullfile(pth, filename1), 'svg');

filename2 = ['LogLog_N',num2str(N),'c',num2str(dens)];
saveas(fig2,fullfile(pth, filename2), 'fig');
saveas(fig2,fullfile(pth, filename2), 'svg');

filename3 = ['Delta_N',num2str(N),'c',num2str(dens)];
saveas(fig3,fullfile(pth, filename3), 'fig');
saveas(fig3,fullfile(pth, filename3), 'svg');

filename4 = ['Logy_Delta_N',num2str(N),'c',num2str(dens)];
saveas(fig4,fullfile(pth, filename4), 'fig');
saveas(fig4,fullfile(pth, filename4), 'svg');


%}

%}
skn = skewness(T_c);

toc  
%{
c=rho/k^2;  
figure
scatter(Gamma_n, omega_n)
%set(gca,'xscale', 'log','yscale','log')
xlabel('\omega_n');
ylabel('\Gamma_n');
title(['\rho=',num2str(c),'; N=',num2str(N), '; b_0=',num2str(b0)])


xy=[x y];
save(['sigma_N' num2str(N) '_R' num2str(rho) '.mat'],'sigma')
save(['xy_N' num2str(N) '_R' num2str(rho) '.mat'],'xy')
save(['rho' num2str(rho) '.mat'],'Data')
%}

%asimetria = [T(Delta)-T(-Delta)]/[T(Delta)+T(-Delta)]

function [x,y,z, DeltaX, DeltaY, DeltaZ, D] = Data(N,z0,R,rm)



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Fill for more than one disorder realization                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %Defining sigma as a function of
        %Rlz: Disorder realizations
        %steps: Detuning
        %NN: number of modes


         %for j=1:Rlz % Disorder realization loop
            clear x y z 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%Random distrution of points in a retangle%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
%uniform distribution
%{
    V = pi*R^2*z0;
    rho = N/(V); 
    r_min = rho^(-1/3)/(.1*rm*pi);         %Minimal distance between two particles

    angle1 = 0;
    angle2 = 2*pi ;           

    x = zeros(N,1);
    y = zeros(N,1);
    z = zeros(N,1);

    theta = (angle2 - angle1) * rand + angle1;
    r = R*sqrt(rand); 
    x(1) = r*cos(theta);
    y(1) = r*sin(theta);
    z(1) = z0*rand - z0/2;


     nn=2;

            for rr = 1:100*N
                theta = (angle2 - angle1) * rand + angle1;
                r = R*sqrt(rand);
                x1 = r*cos(theta);
                y1 = r*sin(theta);
                z1 = z0*rand - z0/2;

                r1 = sqrt((x1-x).^ 2 + (y1-y).^2 + (z1 -z).^2);     %Distance between the test poin and the other points
                test1 = r1 < r_min;                          %Testing if the distance bewteen the points are large enough
                test = any(test1);


                if test == 0
                    x(nn) = x1;
                    y(nn) = y1;
                    z(nn) = z1;
                    nn = nn +1;
                end

                if nn == N +1
                    break
                end

            end
%}

%Gaussian 

rhosq = unifrnd(0,R^2,N,1);
rho1 = sqrt(rhosq);
theta = unifrnd(0,2*pi,N,1);
x = rho1.*cos(theta);
y = rho1.*sin(theta);
z = normrnd(-z0/2,z0/2,N,1);


%rf = (std(x)*std(y))^(1/2); %Raio efetivo
%1/(exp(1)^2)*1/6*rf

clear Repx Repx1 Repy Repy1 Repz Repz1 DeltaX DeltaY Deltaz D

    Repx = repmat(x,[1,N]);
    Repx1 = repmat(transpose(x),[N,1]);

    DeltaX = Repx - Repx1;

    Repy = repmat(y,[1,N]);
    Repy1 = repmat(transpose(y),[N,1]);

    DeltaY = Repy - Repy1;

    Repz = repmat(z,[1,N]);
    Repz1 = repmat(transpose(z),[N,1]);

    DeltaZ = Repz - Repz1;

    D = sqrt(DeltaX.^2 + DeltaY.^2 + DeltaZ.^2);
    %np.fill_diagonal(D, 1)
    
end

function sigma = Sigma(Gamma_0,N,k0,Domega1,W0,E0,x, y, z, DeltaX, DeltaY, DeltaZ, D,vec,beta)

        k=k0;
         %Laser beam fiend in the Particle positions:
         
         [EL1, EL1x, EL1y] = LaserFieldFunction3D(x,y,z,W0,E0,k,beta);         

         EL1z = zeros(N,1);
         


    
                %%%%%%%%%%%%%%%%%%%%%%%%
                %%%     Constructing the Green Matrix      %%%
                %%%%%%%%%%%%%%%%%%%%%%%%



                
                clear A H0 HI Hp Hm 

                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%     Vector Light     %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 if vec == 1
                    %Green matrix: First portion of the sum
                    Gd = 1i*(3*Gamma_0/4).*exp(1i*k0.*D)./(k0.*D).*(1 + 1i./(k0.*D) - 1./( (k0.*D).^2));
                    Gd(isinf(Gd))=0;
                    %Green matrix: Second portion of the sum
                    Gt = 1i*(3*Gamma_0/4).*exp(1i*k0.*D)./(k0.*D).*(-1 - 3i./(k0.*D) + 3./((k0.*D).^2))./(D.^2);
                    Gt(isinf(Gt))=0;
                    %Green Matrix: Positions product (ra*ra')
                    XX = DeltaX.^2;
                    YY = DeltaY.^2;
                    ZZ = DeltaZ.^2;
                    XY = DeltaX.*DeltaY;
                    XZ = DeltaX.*DeltaZ;
                    YZ = DeltaY.*DeltaZ;
                    %Green Matrix: Diagonal
                    Diag = eye(3*N);
                    HI = (1i * Domega1 - Gamma_0/2) .* Diag;
                    %Green Matrix: Assemblong Blocks
                    Zr = zeros(N,N);
                    GGd = [Gd Zr Zr ; Zr Gd Zr ; Zr Zr Gd];
                    GGt = [Gt Gt Gt ; Gt Gt Gt ; Gt Gt Gt];
                    XYZ = [XX XY XZ ; XY YY YZ ; XZ YZ ZZ];

                    A = HI + (GGd + GGt.*XYZ); 
                    E1= 1/2.*[EL1x; EL1y; EL1z];


                 else

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%      Scalar Light     %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  




                    H0 = (-Gamma_0/2).*(exp(1i*k0.*D)./(1i.*k0*D));
                    H0(isnan(H0))=0;
                    H0(isinf(H0))=0;
                    Diag=eye(N);
                    HI=(1i*Domega1-Gamma_0./2).*Diag;
                    H0=H0+HI;


                    %Matrix A
                    %Linear system:
                    %A*sigma=E1

                    A= H0;
                    E1= (1i/2).*EL1;
                 end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     Stationary  solutuion                       %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sigma=linsolve(A,E1);

           


end
