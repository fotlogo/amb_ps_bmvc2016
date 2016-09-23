function [ N ] = visualise_reconstruction(X,Y,Z,C,mask,f,cc,S,Sd,Phi,mu,epsil,mm_to_px, title_str )
%UNTITLED4 Summary of this function goes here
[nrows,ncols]=size(Z);
% [y,x] = meshgrid(1:ncols,1:nrows);
% x=x-cc(1);
% y=y-cc(2);
[x,y] = meshgrid(1:ncols,1:nrows);
x=x-cc(2);
y=y-cc(1);
%for quick check cause of too much confusion :P
xmin=min(x(:)); xmax=max(x(:));
ymin=min(y(:)); ymax=max(y(:));
fprintf(1,'image plane coord lims [%.0f %.0f %.0f %.0f] \n', xmin, xmax, ymin, ymax);
% Normals
% [zy,zx] = gradient(Z);
[zx,zy] = gradient(Z);
N = zeros(nrows,ncols,3);
Nz = 1+Z/f + x.*zx/f + y.*zy/f;
normalisation = sqrt(zx.^2+zy.^2+Nz.^2);
normalisation(mask==0)=NaN;
N(:,:,1) = zy./normalisation;
N(:,:,2) = zx./normalisation;
N(:,:,3) = -Nz./normalisation; 

[ A, H ] = calculate_A_H_fields(X,Y,Z,C,(mask>0),S,Sd,Phi,mu,epsil );

% % Px ==> mm
X= X/mm_to_px;
Y= Y/mm_to_px;
Z= Z/mm_to_px;

C_affichage = 0.6;
figure()
num_image = 2;
h=surf(X,Y,Z); %FOT it was initialy -Z giving negatives
shading flat
colormap gray
 view(0,-45);
% view(-90,0);
Ireproj = ((A(:,:,num_image)).*(max(0,sum(N.*H(:,:,:,num_image),3))).^(1./C_affichage));%
set(h,'CData',(Ireproj))
axis equal  

% Sd=[0.930340420000000,-0.195723600000000,0.305429930000000;
%     0.923262530000000,0.248279910000000,0.289477700000000;
%     0.755365550000000,-0.517115260000000,0.398578590000000;
%     0.0944696210000000,0.494408530000000,0.860295440000000;
%     0.195812240000000,-0.790059490000000,0.577170620000000;
%     -0.465786120000000,0.525195070000000,0.708990210000000;
%     -0.756029390000000,-0.472338480000000,0.448976960000000;
%     -0.800318460000000,-0.154047240000000,0.574830790000000];
% 
% Sd=Sd';

hold on; %%uncoment to debug light source positions 

%   plot3(0,0,0,'wo'); 
%   plot3(0,0,-f/mm_to_px,'k*'); 
%  for ii=0: floor(size(S,2)/6)
%     ss=6*ii+1;
%     if(ss<=size(S,2))
%         plot3( S(1,ss)/mm_to_px, S(2,ss)/mm_to_px, S(3,ss)/mm_to_px, 'r*', 'markers',5*ii+10);
%     end
%     ss=6*ii+2;
% 	if(ss<=size(S,2))
%         plot3( S(1,ss)/mm_to_px, S(2,ss)/mm_to_px, S(3,ss)/mm_to_px, 'g*', 'markers',5*ii+10);
% 	end
% 	ss=6*ii+3;
% 	if(ss<=size(S,2))
%         plot3( S(1,ss)/mm_to_px, S(2,ss)/mm_to_px, S(3,ss)/mm_to_px, 'b*', 'markers',5*ii+10);
% 	end
% 	ss=6*ii+4;
% 	if(ss<=size(S,2))
%         plot3( S(1,ss)/mm_to_px, S(2,ss)/mm_to_px, S(3,ss)/mm_to_px, 'c*', 'markers',5*ii+10);
% 	end
%    ss=6*ii+5;
% 	if(ss<=size(S,2))
%         plot3( S(1,ss)/mm_to_px, S(2,ss)/mm_to_px, S(3,ss)/mm_to_px, 'y*', 'markers',5*ii+10);
% 	end
% 	ss=6*ii+6;
% 	if(ss<=size(S,2))
%         plot3( S(1,ss)/mm_to_px, S(2,ss)/mm_to_px, S(3,ss)/mm_to_px, 'm*', 'markers',5*ii+10);
% 	end
%  end
% for ii=1: size(S,2)
%       plot3( S(1,ii)/mm_to_px, S(2,ii)/mm_to_px, S(3,ii)/mm_to_px, 'b*', 'markers',10);   
%       
% %        for jj=1:1
% %       
% %       plot3( S(1,ii)/mm_to_px+3*jj^2*Sd(1,ii), S(2,ii)/mm_to_px+3*jj^2*Sd(2,ii), S(3,ii)/mm_to_px+3*jj^2*Sd(3,ii), 'g*', 'markers',10);   
% %       
% %        end
% end
% 
%  for ii=1:5
% 	plot3( ii^2, 0,0,'r*');
%     plot3( 0, ii^2, 0, 'g*');
% 	plot3( 0,0, ii^2, 'b*');
%  end
hold off;

title(title_str);

% zu = Z([2:end end],:)-Z;
% zu(mask ==0) = NaN;
% zv = Z(:,[2:end end])-Z;
% zv(mask ==0) = NaN;
% N = zeros(nrows,ncols,3);
% N(:,:,1) = zu;
% N(:,:,2) = zv;
% N(:,:,3) = ((f+Z)./f + (x.*zu+y.*zv)./f);
end

