clc
clear
% set Dc distribution paramter
sigma = 1.3;
CMD = 70;
% set Dc distribution function
n_BC=@(D)1/sqrt(2*pi)./D/log(sigma).*...
    exp(-0.5*((log(D)-log(CMD))/log(sigma)).^2); % D unit:nm
% set Dc distribution function
n_BC=@(D)1/sqrt(2*pi)./D/log(sigma).*...
    exp(-0.5*((log(D)-log(CMD))/log(sigma)).^2); % D unit:nm

% set Dc distribution paramter
k = 0.03;
% set CT distribution function
n_CTF =@(CT)k*exp(-k*CT);

% set bin paramter
bin = 1;
Dc_min=0;
Dc_max=200;
CT_min=0;
CT_max=400;
% set Dc dist
Dc_mid = Dc_min+bin/2:bin:Dc_max;
for i=1:length(Dc_mid)
    n_Dc(i) = n_BC(Dc_mid(i))*bin;
end
N_Dc = sum(n_Dc);

% set CT dist
CT_mid = CT_max-bin/2:-bin:CT_min;
for i=1:length(CT_mid)
    n_CT(i) = n_CTF(CT_mid(i))*bin;
end
N_CT = sum(n_CT);

% get the n(Dc. CT)
n_matrix = repmat(n_Dc,length(n_CT),1).*repmat(transpose(n_CT),1,length(n_Dc));
N=sum(sum(n_matrix));
% plot Figure1
figure;
imagesc(n_matrix);
% add colorbar
cb = colorbar;
cb.Label.String = 'n(Dc,CT)'; % 设置颜色条标签

%colormap(jet);
colormap(parula); % hot还可以 bone还可以
%set(gca,'FontName','Arial','Fontsize',12);
% set label
xlabel('Dc(nm)');
ylabel('CT(nm)');
% set tick
xticks(Dc_min/bin+50/bin:50/bin:Dc_max/bin);
yticks(CT_min/bin:100/bin:CT_max/bin);

% set tickslabel
set(gca, 'XTickLabels', {Dc_min+50:50:Dc_max}, ...
         'FontName', 'Arial', 'FontSize', 12)
%xticklabels({Dc_min+50:50:Dc_max},'FontName','Arial','Fontsize',12);
yticklabels({CT_max:-100:CT_min});

% set subtitle
title('(a)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left', 'FontSize', 12);

% keep squre
%axis equal;

set(gcf, 'unit', 'centimeters', 'position', [0 0 8.5 7])
print('-r1000','-dpng','n(Dc,CT).png');
