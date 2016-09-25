function [X,Y,Z, C_refined] = ambient_ps(I, mask, mean_distance, cam, C, S_struct, epsil,thresholds,refine_C)
%ambient_ps main calculation function
% Author Fotios Logothetis fl302@cam.ac.uk
[nrows,ncols,~] = size(I);

%uncompress parameters
f=cam.f;
cc=cam.cc;

S=S_struct.S;
Sd=S_struct.Sd;
Phi=S_struct.Phi;
Phi=Phi(:);
mu=S_struct.mu;
%%% 
%optimiser params. leave them as they are   
alpha = 50;
minit = 5;
maxit = 20;

tol = 5e-4;
tol_global = 1e-4;
lambda = 1e-9;
nGS = 100;
nloops = 10;

use_L2=1; %use (much) faster L2 optimiser (instead of L1). The L1 optimiser 'should'(in theory) be more accurate

%%%
indices_mask = (mask>0);
%% PREPROCESING
% Pixel coordinates
[x,y] = meshgrid(1:ncols,1:nrows);
x=x-cc(2);
y=y-cc(1);
%for quick check cause of too much confusion :P
xmin=min(x(:)); xmax=max(x(:));
ymin=min(y(:)); ymax=max(y(:));
%it is very common to mess x0,y0 or axis conventions etc. so display this
%and check that it is ok
fprintf(1,'image plane coord lims [%.0f %.0f %.0f %.0f] Do they look ok?\n', xmin, xmax, ymin, ymax);

% Prior shape as a plane at distance mean_distance
Z0 = (mean_distance)*(ones(nrows,ncols));
Z = Z0;
Z(mask==0) = NaN;
res_old = Inf;

for loop = 1:nloops
    Zbefore = Z;
%X,Y,Z in Papadimitri, Favaro Convention    
    X = x.*(f+Z)/f;
    Y = y.*(f+Z)/f;      
   
	[ A, H ] = calculate_A_H_fields(X,Y,Z,C,indices_mask,S,Sd,Phi,mu,epsil);   
    % Create b,s fields
    [b,s]=calculate_b_s_fields(I,indices_mask,x,y,Z,f,C,A,H,thresholds);    
    
    clear A H  %save memory  
%     L1 Opti
    order = 1; % 1 for order-1 finite diff (thin details, less robust), 2 for order-2 (smoother, more robust)    
    %     whos    
	if ((use_L2)||(loop==1))
       
        disp('running L2 optimiser...');
       
        [A_system,mapping_matrix,Omega] = make_Matrix_freeboundary_ICCV(b,lambda,mask,nrows,ncols);        
        b_system = make_SecondMember_freeboundary_ICCV(b,s,lambda,Z,mask,mapping_matrix,Omega);
        clear mapping_matrix Omega
        
        Z(indices_mask) =A_system\b_system;    
        clear A_system b_system    
	else 
            % L1 optim
            disp('SB refinement');
            [Z,~,~,~] = minimisation_L1_bis(b,s,mask,Z0,tol,minit,maxit,lambda,order,alpha,nGS,1,Z);%       
	end      
	Z(mask==0) = NaN;
  %%%%stopping condition on residual   
    residual = norm(Zbefore(indices_mask)-Z(indices_mask))/norm(Z(indices_mask));
    
    disp('=============')
    fprintf(1,'Loop %d : residual = %.08f\n',loop,residual);	
    disp('=============')
    
    if(residual<tol_global)
        break;
    end
    if((residual>res_old))
        break;
    end
    res_old = residual;      
	%%% START C estimation
	if(refine_C)
	[ A, H ] = calculate_A_H_fields(X,Y,Z,C,indices_mask,S,Sd,Phi,mu,epsil);
	%% Estimate N
	zu = Z([2:end end],:)-Z;
	zu(mask ==0) = NaN;
	zv = Z(:,[2:end end])-Z;
	zv(mask ==0) = NaN;
	N = zeros(nrows,ncols,3);
	N(:,:,1) = zu;
	N(:,:,2) = zv;
	N(:,:,3) = -((f+Z)./f + (x.*zu+y.*zv)./f);
	N = N./repmat(sqrt(N(:,:,1).^2+N(:,:,2).^2+N(:,:,3).^2),[1 1 3]);      
	%%
	C=re_estimate_c(I,C,N,H,A,mask,thresholds(1)); 
	end 		
end

%% Kill point at boundaries that are badly constrained
zx1 = abs(Z([3:end end end],:)-Z);
zx2 = abs(Z([1 1 1:end-2],:)-Z);
zy1 = abs(Z(:,[3:end end end])-Z); 
zy2 = abs(Z(:,[1 1 1:end-2])-Z); 

boundaries= isnan(zx1) | isnan(zx2) | isnan(zy1) | isnan(zy2);
% no_change= (abs(Z0-Z)<2);
% Z(boundaries & no_change)=NaN;
Z(boundaries)=NaN;

C_refined=C;
end



