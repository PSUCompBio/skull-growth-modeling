clear
clc

fid = fopen('E15_inner.surf','r')
for k=1:18
    tline = fgets(fid);
end
A = fscanf(fid,'%f %f %f',[3,8816]);

for i=1:12
    tline = fgets(fid);
end
B = fscanf(fid,'%f %f %f',[3,17628]);

fclose(fid)


fid3 = fopen('E13_inner.surf','r')
for k=1:18
    tline = fgets(fid3);
end
D = fscanf(fid3,'%f %f %f',[3,7814]);

for i=1:12
    tline = fgets(fid3);
end
E = fscanf(fid3,'%f %f %f',[3,15624]);

fclose(fid3)

FV.vertices = D';
FV.faces = E'

points = A';

%      FV.faces    = [5 3 1; 3 2 1; 3 4 2; 4 6 2];
%      FV.vertices = [2.5 8.0 1; 6.5 8.0 2; 2.5 5.0 1; 6.5 5.0 0; 1.0 6.5 1; 8.0 6.5 1.5];
%      points      = [2 4 2; 4 6 2; 5 6 2];
      [distances,surface_points] = point2trimesh(FV, 'QueryPoints', points); 
      patch(FV,'FaceAlpha',.5); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
      plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
      plot3M(points,'*r')
      plot3M(surface_points,'*k')
      plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1);shiftdim(points,-1)*NaN],[],3),'g')

v = points - surface_points;
rev_v =-v;
final = [surface_points v];
final_rev = [points rev_v];
final_data = final';

fid4=fopen('XDX_E13_E15_inner','w');
fprintf(fid4,'%d %d %d %d %d %d \n',final_data)
fclose(fid4)
