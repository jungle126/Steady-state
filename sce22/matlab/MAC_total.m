function MAC = MAC_total(Dp_list,Dc_list,Conc_list)
%MAC_TOTAL calculate MAC of a time

%   set RI BC_density and wavelength
RI_BC= complex(1.95,0.96);  
RI_shell= complex(1.53,0);
BC_den= 1800; 
wavelength = 633*1e-9;
Dp_m = Dp_list*1e-9;  % nm to m
Dc_m = Dc_list*1e-9;   % nm to m
for j=1:length(Dc_m)
    xcor = pi*Dc_m(j)/wavelength;
    xman = pi*Dp_m(j)/wavelength;
    result=Miecoated(RI_BC,RI_shell,xcor,xman,1);
    qabs=result(3);
    cabs_shell(j)=1/4*pi*Dp_m(j).^2 * qabs * Conc_list(j);
    mass_bc(j)=1/6*pi*Dc_m(j).^3 * BC_den * Conc_list(j);
end
% calculate MAC
MAC = sum(cabs_shell)/sum(mass_bc)*1e-3;  %m^2/g


