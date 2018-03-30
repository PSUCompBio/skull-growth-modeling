clear
clc

data=load('XDX_E13_E15_outer_wall');
X=0.001*[data(:,1),data(:,2),data(:,3)];
DX=0.001*172.8/172800*[data(:,4),data(:,5),data(:,6)];
Xq = load ('Xq_E13_32_outer_wall');
%DXq = griddata(XDX(:,1),XDX(:,2),XDX(:,3),XDX(:,4),Xq(:,1),Xq(:,2),Xq(:,3));
Fx = scatteredInterpolant(X(:,1),X(:,2),X(:,3),DX(:,1), 'linear', 'nearest');
Fy = scatteredInterpolant(X(:,1),X(:,2),X(:,3),DX(:,2), 'linear', 'nearest');
Fz = scatteredInterpolant(X(:,1),X(:,2),X(:,3),DX(:,3), 'linear', 'nearest');
DXq = Fx(Xq(:,1),Xq(:,2),Xq(:,3));
DYq = Fy(Xq(:,1),Xq(:,2),Xq(:,3));
DZq = Fz(Xq(:,1),Xq(:,2),Xq(:,3));
Dq = [DXq';DYq';DZq'];

fid1=fopen('Dx_E13_32_E15_outer_wall','w');
fprintf(fid1,'( %d %d %d) \n',Dq)
fclose(fid1)

figure(2)
quiver3(Xq(:,1),Xq(:,2),Xq(:,3),DXq, DYq, DZq,2)