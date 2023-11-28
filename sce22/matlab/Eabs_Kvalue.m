%RI, wavelength, and BC_den
RI_BC= complex(1.95,0.96);  
RI_shell= complex(1.53,0);
dens_BC= 1800; 
wavelength = 633*1e-9;
sigma= 1.3;
CMD=70;
n_BC=@(D)1/sqrt(2*pi)./D/log(sigma).*...
    exp(-0.5*((log(D)-log(CMD))/log(sigma)).^2); % D unit:nm
%set Dc list
Dc = 1:1000; %nm
Dc_m = Dc*1e-9;
Dp_m = (Dc+49)*1e-9;
n_sum=0;
for j=1:length(Dc)
    xcor = pi*Dc_m(j)/wavelength;
    xman = pi*Dp_m(j)/wavelength;
    result=Miecoated(RI_BC,RI_shell,xcor,xman,1);
    result_ex=Miecoated(RI_BC,RI_shell,xcor,xcor,1);
    qext=result(1);
    qabs=result(3);
    qext_ex=result_ex(1);
    qabs_ex=result_ex(3);
    cabs_shell(j)=1/4*pi*Dp_m(j).^2 * qabs * n_BC(Dc(j));
    cabs_core(j)=1/4*pi*Dc_m(j).^2 * qabs_ex * n_BC(Dc(j));
    mass_bc(j)=1/6*pi*Dc_m(j).^3 * dens_BC * n_BC(Dc(j));
    n_sum=n_sum+n_BC(Dc(j));
end
accum_cabs_shell=trapz(Dp_m,cabs_shell);
accum_cabs_core=trapz(Dc_m,cabs_core);
accum_massbc1=trapz(Dc_m,mass_bc);

% end
mac_shell = accum_cabs_shell./accum_massbc1 .* 1e-3; %m2/g;
mac_core = accum_cabs_core./accum_massbc1 .* 1e-3;
Eabs = mac_shell/mac_core;
