function [X,Y,Z, C_refined] = ambient_ps(I, mask, mean_distance, cam, C, S_struct, epsil,thresholds,refine_C,ambient)
%perspective_ps Summary of this function goes here
%   Detailed explanation goes here
[nrows,ncols,nb_images] = size(I);

%uncompress parameters
f=cam.f;
cc=cam.cc;

S=S_struct.S;
Sd=S_struct.Sd;
Phi=S_struct.Phi;
Phi=Phi(:);
mu=S_struct.mu;
%%% 
%optimiser    
alpha = 50;
minit = 5;
maxit = 20;
display = 1; % 0 = no display (fast), 1 = display CV residual, 3 = display maps and CPU times (VERY SLOW), 4 = display 3D (VERY VERY SLOW)
tol = 5e-4;
tol_global = 1e-4;
lambda = 1e-9;
nGS = 100;
nloops = 10;
alpha_admm = 1;
% maxit_admm = 10;
% tol_admm = 1e-2;
% Lp = 1;
fast=0;
silent=0;
%%%
indices_mask = (mask>0);
%%
% Pixel coordinates
[x,y] = meshgrid(1:ncols,1:nrows);
x=x-cc(2);
y=y-cc(1);
%for quick check cause of too much confusion :P
xmin=min(x(:)); xmax=max(x(:));
ymin=min(y(:)); ymax=max(y(:));
fprintf(1,'image plane coord lims [%.0f %.0f %.0f %.0f] \n', xmin, xmax, ymin, ymax);

% Prior shape as a plane at distance mean_distance
% rz=0.7*randn(nrows/10,ncols/10);
% rz=imresize(rz,[nrows,ncols]);
Z0 = (mean_distance)*(ones(nrows,ncols));
Z = Z0;
Z(mask==0) = NaN;
res_old = Inf;

 use_A_table=0;
% % error_old=Inf;
for loop = 1:nloops
    Zbefore = Z;
%X,Y,Z in Papadimitri, Favaro Convention    
    X = x.*(f+Z)/f;
    Y = y.*(f+Z)/f;  
    
    if use_A_table
        [H ] = calculate_H_field(X,Y,Z,C,indices_mask,S,epsil );
        A=getAttenuation(X/mm_to_px,Y/mm_to_px,Z/mm_to_px,1:nb_images);         

%         figure(10);
%         imshow(A(:,:,1));
%         colormap(10,'jet');
%         
%         A=imresize(A,[nrows,ncols]/1.1);
%         A=imresize(A,[nrows,ncols]);
%         
%         figure(20);
%         imshow(A(:,:,1));
%         colormap(20,'jet');
%         pause
        
    else
        [ A, H ] = calculate_A_H_fields(X,Y,Z,C,indices_mask,S,Sd,Phi,mu,epsil);  
    end   
    
%     A(:,:,16:23)=0.001*A(:,:,16:23);
%      [ A, H ] = calculate_A_H_fields(X,Y,Z,C,indices_mask,S,Sd,Phi,mu,epsil);  
%       A=0.5*A/nanmean(A(:));
% %      A=0.5*A/nanmax(A(:));
%      AA=getAttenuation(X/mm_to_px,Y/mm_to_px,Z/mm_to_px,1:nb_images);   
%      AA=0.5*AA/nanmean(AA(:));
%      
%  for k=1:24
%     
% %     figure(3);
%      I1=AA(:,:,k);
% %     imshow(I1);
% %     title('table');
% %     colormap(3, jet);
% % %     
%     I2=A(:,:,k);
% % %     figure(4);
% % %     imshow(I2);
% % %     title('model');
% % %     colormap(4, jet);
% % %     
%      figure(k);
%      imshow(0.5*I1./I2);   
%      colormap(k, jet);
%  end   
%      pause
%        A=getAttenuation(X/mm_to_px,Y/mm_to_px,Z/mm_to_px,1:nb_images);    
    % Create b,s fields
    [b,s]=calculate_b_s_fields(I,indices_mask,x,y,Z,f,C,A,H,thresholds,ambient);     
%     mean(b(:))
%     mean(s(:))
%     pause
    
%     nxym=zeros(nrows,ncols,3);    
%    for ii=1:nrows
%         for jj=1:ncols            
%             if(mask(ii,jj)==0)
%                 continue;
%             end
%        AA=squeeze(b(ii,jj,:,:));
%        BB=squeeze(s(ii,jj,:));
% %        size(AA)
% %        size(BB)
%        nxy=BB'/AA';
%        nn=[nxy';-1];
%        nn=nn/norm(nn);
%         nxym(ii,jj,:)=(nn+1)/2;        
%         end
%    end    
%    figure;
%     imshow(nxym);     
%     pause      
    
    clear A H    
%     L1 Opti
    order = 1; % 1 for order-1 finite diff (thin details, less robust), 2 for order-2 (smoother, more robust)    
    %     whos    
    if (fast==0)
        if silent==0
        disp('constructing system part 1...');
        end
        [A_system,mapping_matrix,Omega] = make_Matrix_freeboundary_ICCV(b,lambda/alpha_admm,mask,nrows,ncols);
        if silent==0
        disp('constructing system part 2...');
        end
        b_system = make_SecondMember_freeboundary_ICCV(b,s,lambda/alpha_admm,Z,mask,mapping_matrix,Omega);
        clear mapping_matrix Omega
        if silent==0
        disp('solving system ...');
        end       
        Z(indices_mask) =A_system\b_system;    
        clear A_system b_system
    elseif (fast==1)
        disp('LS initialization');     
        [Z,~,~] = minimisation_L2(b,s,mask,Z,lambda,order);
    else    
        % L2 initialization
        if(loop==1)
            disp('LS initialization');
            [Z,~,~] = minimisation_L2(b,s,mask,Z,lambda,order);
        else
            % L1 optim
            disp('SB refinement');
            [Z,~,~,~] = minimisation_L1_bis(b,s,mask,Z0,tol,minit,maxit,lambda,order,alpha,nGS,display,Z);%
       end 
    end  
    
	Z(mask==0) = NaN;
% 	zu = Z([2:end end],:)-Z;
% 	zu(mask ==0) = NaN;
% 	zv = Z(:,[2:end end])-Z;
% 	zv(mask ==0) = NaN;
% 	error_current = sum((bsxfun(@times,zu,b(:,:,:,1))+bsxfun(@times,zv,b(:,:,:,2))-s).^2,3);
% 	error_current = error_current/nb_images;
% % 	erreur = erreur+lambda*((Z-Z0).^2);
% 	error_current = nanmean(error_current(:));	   
%     
%      fprintf(1,'Loop %d : error = %f\n',loop, error_current);	    
%     if error_current> error_old
%          Z=Zbefore;
%         break;
%     end
%     
%     if (error_old-error_current)/error_old < tol_global
%        break;
%     end
%     error_old=error_current;
%  %%%%stopping condition on residual   
    residual = norm(Zbefore(indices_mask)-Z(indices_mask))/norm(Z(indices_mask));
    if silent==0
    disp('=============')
    fprintf(1,'Loop %d : residual = %.08f\n',loop,residual);	
    disp('=============')
    end
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
            C=re_estimate_c(I,C,N,H,A,mask,thresholds(1),ambient); 
        end        
		%%% END C estimation
end

%% Kill point at boundaries that remain unchanged
% indx_j=find_jumps(Z);
% fprintf(1,'removing %d points because of z jump \n', sum(indx_j(:)));
% Z(indx_j)=NaN;
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



