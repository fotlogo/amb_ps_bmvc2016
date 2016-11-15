%% Main test code for the paper 'Near-Field Photometric Stereo in Ambient Light'. from BMVC 2016
%% Code for 'A Single Lobe Photometric Stereo Approach for Heterogeneous Materia' from SIAM 2016 also included
% Author Fotios Logothetis fl302@cam.ac.uk
clear all 
close all

results_dir='./results/';
data_dir='./data/';

  name='buddha'; 
%  name='ArlequinMask';

data_f=[data_dir,name,'.mat'];  
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
%--Sd: 3 x nimages vector containing light source maximum illumination
%directions. If mu=0 it does not matter (as long as it non-zero). if not
%sure use [0;0;1] for each source
%--AMB:(OPTIONAL).used to compare with the old, non-ambient version

load(data_f) 
%% SOME TESTS
[nrows,ncols,nimages] = size(I);
assert(size(S,1)==3);
assert(size(S,2)==nimages);
assert(size(Sd,1)==3);
assert(size(Sd,2)==nimages);
assert(size(Phi,1)==nimages);
assert(size(mu,1)==nimages);
%% RESIZE
%This helps reducing computation time and RAM
%The L1 optimiser will need 10+GB of RAM if running at full resolution
ratio =8;
I=imresize(I,[nrows,ncols]/ratio); 
mask=imresize(mask,[nrows,ncols]/ratio);
AMB=imresize(AMB,[nrows,ncols]/ratio);
S=S/ratio;

f = f/ratio;
cc = cc/ratio;
mean_distance =mean_distance/ratio;
mm_to_px =mm_to_px/ratio;

[nrows,ncols,nb_images] = size(I);
I = max(0.01,I);
%% group vars
cam.f=f;
cam.cc=cc;

S_struct.S=S;
S_struct.Sd=Sd;
S_struct.Phi=Phi;
S_struct.mu=mu;
%% Misk opts
use_L2=0; %use (much) faster and less memory consuming L2 optimiser (instead of L1). The L1 optimiser 'should' be more accurate
C =1*ones(size(I,1),size(I,2));  %initialise C as being Lambertian.
refine_C=1;
shadow_threshold = 0.03; 
saturation_thress=0.99;
thresholds=[shadow_threshold,saturation_thress];

ambient=1;
%% run main function-AMBIENT
[XA,YA,ZA, C_refined] = perform_ps(I, mask, mean_distance,cam, C,S_struct,thresholds,use_L2,refine_C,ambient);

mask_out=mask;
mask_out(isnan(ZA))=0;% 
%% Display 
title_str=sprintf('Ambient perspective PS with %d images', size(S,2));
[ ~ ] = visualise_reconstruction(XA,YA,ZA,C_refined,mask_out,f,cc,S,Sd,Phi,mu,mm_to_px,title_str ); 

% figure;
% imshow(C_refined);
% title('estimated c map');

XYZ = cat(3,XA,YA,ZA)/mm_to_px;
export_ply(XYZ,mask_out,[results_dir,name,'.ply']);
%%
%% run main function-DARK to compare
if exist('AMB','var')
%remove ambient else reconstruction will be very flat
%of course if data was captured under no ambient just run directly
Ia=I-repmat(AMB,1,1,size(I,3));
Ia=max(Ia,0.01);
ambient=0;
[XD,YD,ZD, C_refined] = perform_ps(Ia, mask, mean_distance,cam, C,S_struct,thresholds,use_L2,refine_C,ambient);

mask_out=mask;
mask_out(isnan(ZA))=0;% 
%% Display 
title_str=sprintf('Perspective PS (SIAM) with %d images', size(S,2));
[ ~ ] = visualise_reconstruction(XD,YD,ZD,C_refined,mask_out,f,cc,S,Sd,Phi,mu,mm_to_px,title_str ); 
end