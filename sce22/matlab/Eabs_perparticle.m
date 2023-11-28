clc;
clear;
% set path
%csv_file_path = './Dp_Dc_Conc.csv'; % window
csv_file_path = '../Data/Dp_Dc_Conc.csv'; % project
% read data
data_matrix = csvread(csv_file_path);
% introdutction£º
% input matrix: 3*240£¬3 mean: Dp, Dc and Conc£¬240 means 240 time
% output matrix: MAC_core, MAC_shell and Eabs

% set n=time=240
n=240;
% MAC_shell
MAC_A = zeros(1,n);
% MAC_core
MAC_EA = zeros(1,n);
% Eabs
Eabs_A = zeros(1,n);

for i=1:240
    MAC_A(i)=MAC_total(data_matrix(i,:),data_matrix(i+240,:),data_matrix(i+480,:));
    MAC_EA(i)=MAC_total(data_matrix(i+240,:),data_matrix(i+240,:),data_matrix(i+480,:)); %Dc=Dp
end
% calculate Eabs
Eabs_A = MAC_A./MAC_EA;

MAC_ABC_matrix = [MAC_A;MAC_EA;Eabs_A];
%csvwrite('./MAC_MACE_Rabs_MACEX_ABC_time.csv', MAC_ABC_matrix);%window
csvwrite('../Data/MAC_MACE_Rabs_MACEX_ABC_time.csv', MAC_ABC_matrix);%project
