clc
clear

% set Dc distribution paramter
sigma = 1.3; 
CMD = 70;  %70 nm
% set Dc distribution function
n_BC=@(D)1/sqrt(2*pi)./D/log(sigma).*...
    exp(-0.5*((log(D)-log(CMD))/log(sigma)).^2); % D unit:nm

% set Dc distribution paramter k
k = 0.021;
% set CT distribution function
n_CTF =@(CT)k*exp(-k*CT);

% set bin paramter
bin = 1; % set bin of Dc and CT
Dc_min=0;  % set the range of Dc
Dc_max=300;
CT_min=0;  % set the range of CT
CT_max=400;
% set Dc dist
Dc_mid = Dc_min+bin/2:bin:Dc_max;
n_Dc = zeros(1,length(Dc_mid));

for i=1:length(Dc_mid)
    n_Dc(i) = n_BC(Dc_mid(i))*bin;
end
N_Dc = sum(n_Dc);

% set CT dist
CT_mid = CT_max-bin/2:-bin:CT_min;
n_CT = zeros(1,length(CT_mid));

for i=1:length(CT_mid)
    n_CT(i) = n_CTF(CT_mid(i))*bin;
end
N_CT = sum(n_CT);

% get the n(Dc. CT)
n_matrix = repmat(n_Dc,length(n_CT),1).*repmat(transpose(n_CT),1,length(n_Dc));
N=sum(sum(n_matrix));

VF_DcCT = zeros(length(CT_mid),length(Dc_mid));

% get the VF(Dc,CT)
for i =1:length(Dc_mid)
    Dc  = Dc_mid(i);
    for j =1:length(CT_mid)
        CT = CT_mid(j);
        VF_DcCT(j,i)= Dc^3/(Dc+CT)^3;
    end
end

% get the VF(Dc)
V_core = zeros(1,length(Dc_mid));
V_P = zeros(1,length(Dc_mid));
VF_Dc = zeros(1,length(Dc_mid));
for i =1:length(Dc_mid)
    Dc  = Dc_mid(i);
    V_core(i) = Dc^3;
    V_P(i) = 0;
    for j =1:length(CT_mid)
        CT = CT_mid(j);
        V_P_j = (Dc+CT)^3 * n_CT(j);
        V_P(i) = V_P(i)+V_P_j;
    end
    VF_Dc(i) = V_core(i)/V_P(i);
end

%Dc(nm) to mass(pg)
rho_BC = 1800;  %kg/m3
rho_NBC = 1500;  %kg/m3
mass_pg =@(D)D^3*pi/6*10^(-12)* rho_BC;
Dc_mass_mid = zeros(1,length(Dc_mid));
for i=1:length(Dc_mid)
    Dc_mass_mid(i) = mass_pg(Dc_mid(i));
end
% get mass
V_BC = zeros(1,length(Dc_mass_mid));
mass_BC = zeros(1,length(Dc_mass_mid));
V_NBC = zeros(1,length(Dc_mass_mid));
mass_NBC = zeros(1,length(Dc_mass_mid));
mass_P = zeros(1,length(Dc_mass_mid));
for i =1:length(Dc_mass_mid) 
    V_BC(i) = Dc_mid(i)^3*pi/6*n_Dc(i);
    mass_BC(i) = V_BC(i)*rho_BC;
    V_P(i) = V_BC(i)/VF_Dc(i);
    V_NBC(i) =  V_P(i)- V_BC(i);
    mass_NBC(i) = V_NBC(i)*rho_NBC;
    mass_P(i) = mass_BC(i) + mass_NBC(i);
end
to1MAX = max(mass_P);
mass_BC = mass_BC./to1MAX;
mass_P = mass_P./to1MAX;
mass_P_uniform = mass_BC.*sum(mass_P)./sum(mass_BC);


figure;
semilogx(Dc_mass_mid, mass_BC, 'k', 'LineWidth', 2);
hold on;
semilogx(Dc_mass_mid, mass_P, 'r', 'LineWidth', 2);
hold on;
semilogx(Dc_mass_mid, mass_P_uniform, 'r--', 'LineWidth', 2);
set(gca,'LineWidth',0.9,'FontName','Arial','FontSize', 10);
% set legend
set(gca,'xminortick','off');
tickLength = [0.01, 0.02]; 
set(gca, 'TickLength', tickLength);
legend('BC', 'Particle', 'Particle (uniform)');
legend('boxoff')
xlabel('Mass of BC contained in each particle(pg)');
ylabel(sprintf('The distribution of aerosol components\nwith respect to per-particle BC mass'));
% Set the boundary values for the x-axis
x_min = 1e-5;  % Minimum value for the x-axis
x_max = 1e-2;  % Maximum value for the x-axis
xlim([x_min, x_max]);
% Set the boundary values for the y-axis
y_min = 0;     % Minimum value for the y-axis
y_max = 1.2;   % Maximum value for the y-axis
ylim([y_min, y_max]);
% set subtitile
title('(c)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left', 'FontSize', 12);

set(gcf, 'unit', 'centimeters', 'position', [0 0 17.78 8])
print('-r1000','-dpng','../Figure/Figure4_c.png');


csv_file_path = '../Data/VF_massBC.csv';
% read csv
data_matrix = csvread(csv_file_path);
Dc_mass_mid_P = data_matrix(2,:);
VF_Dc_P = data_matrix(1,:);
idx = VF_Dc_P ~= 0;
VF_Dc_P = VF_Dc_P(idx);
Dc_mass_mid_P = Dc_mass_mid_P(idx);

figure;
semilogx( Dc_mass_mid_P, VF_Dc_P , 'k',  'LineWidth', 1,'Marker','*');
hold on;
semilogx( Dc_mass_mid, VF_Dc , 'k','LineWidth', 2);
set(gca,'LineWidth',1,'FontName','Arial','FontSize', 10);
set(gca,'xminortick','off');
tickLength = [0.02, 0.03]; 
set(gca, 'TickLength', tickLength);
% fill
threshold = 0.5;
x_fill = [Dc_mass_mid_P, fliplr(Dc_mass_mid_P)]; % 
y_fill = [VF_Dc_P, zeros(1, length(VF_Dc_P))]; % 
y_fill2 = [VF_Dc_P, ones(1, length(VF_Dc_P))];
%y_fill(VF_Dc_P <= threshold) = threshold; % 
fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3); %
fill(x_fill, y_fill2, 'r', 'FaceAlpha', 0.3); %
legend('Particle resolved model','Theoretical model' ,'BC','non-BC');
legend('boxoff')
% set label
xlabel('Mass of BC contained in each particle(pg)');
ylabel(sprintf('Volume fraction'));
% Set the boundary values for the x-axis
x_min = min(Dc_mass_mid_P);  % Minimum value for the x-axis
x_max = max(Dc_mass_mid_P);  % Maximum value for the x-axis
xlim([x_min, x_max]);
% Set the boundary values for the y-axis
y_min = 0;     % Minimum value for the y-axis
y_max = 1;   % Maximum value for the y-axis
ylim([y_min, y_max]);
% set subtitile
title('(b)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left', 'FontSize', 12);
set(gcf, 'unit', 'centimeters', 'position', [0.5 0 8.5 7])
print('-r1000','-dpng','../Figure/Figure4_b.png');

