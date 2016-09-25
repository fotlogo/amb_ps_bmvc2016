%% Main test code for the paper 'Near-Field Photometric Stereo in Ambient Light'. from BMVC 2016
% Author Fotios Logothetis fl302@cam.ac.uk
clear all 
close all

results_dir='./results/';
data_dir='./data/';

data_f=[data_dir,'buddha.mat'];      
%  data_f=[data_dir,'head_far.mat'];

%% Load data
%Data MUST contain: 
%--I: nrows x ncols x nimages . image raw data. values assumed to be
% NORMALISED in [0,1]
%--mask: nrows x ncols. forground mask
%--S: 3 x nimages vector containing light source position in pixels. It is
%VERY IMPORTANT to get this right otherwise the reconstruction will fail.
%run the test_axis script to see a demonstration of the convention assumed
%--f: camera focal length in pixels
%--cc: [y0,x0]. camera principal point
%--Phi: nimages x 1 vector containing light source luminances (e.g.
%measured with lux meter). Scale does not matter due to image ratios so
%keep them normalised (in (0,1]) for better numerical stability
%--mean_distance: mean distance between camera plane and object (in
%pixels). Just a rough estimate (e.g. measure with a ruler) should be ok
%--mu: nimages x 1 vector containing light source radial attenuation
%params. If not sure just set the zeros(nimages,1)
%--Sd: nimages x 3 vector containing light source maximum illumination
%directions. If mu=0 it does not matter (as long as it non-zero). if not
%sure use [0;0;1] for each source
load(data_f)   

% mu = 0.5*ones(24,1); 
% Sd=zeros(3,size(I,3));
% Sd(3,:)=1;
% save(data_f,'I','mask','mean_distance','S','Phi','f','cc', 'mm_to_px', 'mu', 'Sd');
% return;  
C =1*ones(size(I,1),size(I,2)); 
%% SELECT ONLY A FEW IMAGES
%    images=[1;3;5;7;10;11;13;16;17;20;22;23]; 
% %    images=1:2:24;
% % % % %      
% I=I(:,:,images);
% S=S(:,images);
% Phi=Phi(images);
% mu=mu(images);
% Sd=Sd(:,images);
%% RESIZE
%This helps reducing computation time and RAM
[nrows,ncols,nimages] = size(I);
ratio =4;
I=imresize(I,[nrows,ncols]/ratio); 
mask=imresize(mask,[nrows,ncols]/ratio);
C=imresize(C,[nrows,ncols]/ratio);
S=S/ratio;

f = f/ratio;
cc = cc/ratio;
mean_distance =mean_distance/ratio;
mm_to_px =mm_to_px/ratio;

[nrows,ncols,nb_images] = size(I);
I = max(0.01,I);
 
epsil = 0.01;   % Lambertian / specular weight
refine_C=1;
%% group vars
cam.f=f;
cam.cc=cc;

S_struct.S=S;
S_struct.Sd=Sd;
S_struct.Phi=Phi;
S_struct.mu=mu;

shadow_threshold = 0.01; 
saturation_thress=0.99;
thresholds=[shadow_threshold,saturation_thress];
%% run main function
[XA,YA,ZA, C_refined] = ambient_ps(I, mask, mean_distance,cam, C,S_struct, epsil,thresholds,refine_C);

C=C_refined;
mask_out=mask;
mask_out(isnan(ZA))=0;% 

%% Display 
title_str=sprintf('Ambient perspective PS with %d images', size(S,2));
[ N ] = visualise_reconstruction(XA,YA,ZA,C,mask_out,f,cc,S,Sd,Phi,mu,epsil,mm_to_px,title_str ); 

figure;
imshow(C_refined);
title('estimated c map');

XYZ = cat(3,XA,YA,ZA)/mm_to_px;
export_ply(XYZ,mask_out,[results_dir,'buddha.ply']);

