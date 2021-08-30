clc
clear all

A=xlsread('position_1hr_1.xlsx');
Nodes=11; % number of nodes
Node_Corr=zeros(Nodes,Nodes);

for n1=1:Nodes
    for n2=1:n1
        n1_xyz=A(4:end,3*(n1-1)+2:3*(n1-1)+4);     
        n2_xyz=A(4:end,3*(n2-1)+2:3*(n2-1)+4);
        Node_Corr(n1,n2)=corr3D_ekta(n1_xyz,n2_xyz);
    end
end

figure
imagesc(Node_Corr)
caxis([-1 1])
colorbar
%axis off

%xlswrite ('position_1hr_1.xlsx',Node_Corr)
