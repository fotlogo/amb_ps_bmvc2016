clear all %TODO improve this and reuse into others
close all

results_dir='./results/';
data_dir='./data/';

data_f=[data_dir,'buddha.mat'];      
%  data_f=[data_dir,'head_far.mat'];

load(data_f)   

C =1*ones(size(I,1),size(I,2)); 
mu = 0.5*ones(24,1); 
% mu = [2*ones(12,1); 0.2*ones(12,1); ];
Sd=zeros(24,3);
Sd(:,3)=1;

%% INCREASE AMB FOR statue_24_ambient
I=I+1.5*repmat(AMB,1,1,size(S,2));
AMB=2.5*AMB;

I = max(0,I);
%% SELECT ONLY A FEW IMAGES
% % % %   images=[1;12;16;19;16;12;1]; %4 images that work
%   images=[1;3;6;9]; %for iso-depth
% % % % %      
% I=I(:,:,images);
% S=S(:,images);
% Phi=Phi(images);
% mu=mu(images);
% Sd=Sd(images,:);
%% RESIZE
[nrows,ncols,nimages] = size(I);
ratio =8;
I=imresize(I,[nrows,ncols]/ratio); 
mask=imresize(mask,[nrows,ncols]/ratio);
C=imresize(C,[nrows,ncols]/ratio);
AMB=imresize(AMB,[nrows,ncols]/ratio);
% AMB=AMB/(256*256); %FIX DATA
S=S/ratio;

f = f/ratio;
cc = cc/ratio;
mean_distance =mean_distance/ratio;
mm_to_px =mm_to_px/ratio;

[nrows,ncols,nb_images] = size(I);
I = max(0.01,I);
%
figure(1);
imshow(AMB);
title('Ambient');

% IM=max(I,[],3);
% Im=mean(I,3);
% 
% %   C=(0.98-2*((IM-Im).^2));
% % % C(mask==0)=NaN;
% % 
%  C=imresize(C,[nrows,ncols]/4);
%  C=imresize(C,[nrows,ncols]);
%  Ia=I-1*repmat(AMB,1,1,size(S,2));
%  Ia=max(Ia,0);
 
epsil = 0.01;   % Lambertian / specular weight
refine_C=1;
ambient=1; 
%%
cam.f=f;
cam.cc=cc;

S_struct.S=S;
S_struct.Sd=Sd;
S_struct.Phi=Phi;
S_struct.mu=mu;

shadow_threshold = 0.01; 
saturation_thress=0.99;
thresholds=[shadow_threshold,saturation_thress];

[XA,YA,ZA, C_refined] = ambient_ps(I, mask, mean_distance,cam, C,S_struct, epsil,thresholds,refine_C,ambient);

C=C_refined;
mask_out=mask;
mask_out(isnan(ZA))=0;% 

% 
title_str=sprintf('Ambient perspective PS(new method) with %d images', size(S,2));
[ N ] = visualise_reconstruction(XA,YA,ZA,C,mask_out,f,cc,S,Sd,Phi,mu,epsil,mm_to_px,title_str ); 

figure;
imshow(C_refined);
title('refined c map');

XYZ = cat(3,XA,YA,ZA)/mm_to_px;
export_ply(XYZ,mask_out,[results_dir,'buddha.ply']);

