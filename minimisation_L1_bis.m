function [zm,d,e,it_breg] = minimisation_L1_bis(b,s,mask,z0,tol,minit,maxit,lambda,alpha,display,zinit)
% Input : 
% b : nrows x ncols x m x 2 -- REQUIRED
% s : nrows x ncols x m -- REQUIRED
% (mask : nrows x ncols) -- DEFAULT : ones(nrows,ncols)
% (z0 : nmaskx1) -- DEFAULT : zeros(nrows,ncols)
% (tol : 1x1) -- DEFAULT : 1e-5
% (maxit : 1x1) -- DEFAULT : 100
% (lambda : 1x1) -- DEFAULT : 1e-9
% (alpha : 1x1) -- DEFAULT : 1e-1
% (nGS : 1x1) -- DEFAULT : 50
% (display : 1x1) -- DEFAULT : 0
% Output : 
% z : nmaskx1
%% TODO vectorise everything else
    [nrows,ncols,m,~] = size(b);
	% Create masks
	disp('Setting variables');
    
	indices_mask = find(mask>0);    
    
% 	indices_nmask = mask==0;
	nb_pixels_mask = length(indices_mask);
    %to know when are we doing fw or backwards differences
    b1f=[-1,-1,+1,+1];
    b2f=[-1,+1,+1,-1];  
    %%
	b1=b(:,:,:,1);
	b2=b(:,:,:,2);
    
    sum_b1_2 = sum(b1.^2,3);
	sum_b2_2 = sum(b2.^2,3); 
    sum_b1_b2 = sum(b1.*b2,3);
 
	[I, J, indices_centre, I_xy, I_x1_y, I_x_y1,Gx,Gy] = make_system_indices_v2(mask);
    
    K=zeros(nb_pixels_mask+9*(size(I_xy{1},1)+size(I_xy{2},1)+size(I_xy{3},1)+size(I_xy{4},1)),1);
    K(1:nb_pixels_mask)=2*lambda*ones(nb_pixels_mask,1)/alpha;   
  
    curr_indx=nb_pixels_mask+1;
    
    for indx=1:4
        N_omega_current=size(I_xy{indx},1);
        % Line of the central pixels
        values = 2*(sum_b1_2 + 2*sum_b1_b2*b1f(indx)*b2f(indx) + sum_b2_2);
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        values =-2*(sum_b1_2 + sum_b1_b2*b1f(indx)*b2f(indx));
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        values =  -2*(sum_b2_2 + sum_b1_b2*b1f(indx)*b2f(indx));
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;	
        % Line of the x-neighbor pixels
        values = -2*(sum_b1_2 + sum_b1_b2*b1f(indx)*b2f(indx));
	    K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        values = 2*sum_b1_2;
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        values = 2*sum_b1_b2*b1f(indx)*b2f(indx);
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        % Line of the y-neighbor pixels
        values = -2*(sum_b2_2 + sum_b1_b2*b1f(indx)*b2f(indx));
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        values = 2*sum_b1_b2*b1f(indx)*b2f(indx);
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
	
        values = 2*sum_b2_2;
        K(curr_indx:curr_indx+N_omega_current-1) = values(indices_centre{indx});
        curr_indx=curr_indx+N_omega_current;
    end
    
	B = zeros(nb_pixels_mask,1);

	% Construction de A : 
	A = sparse(I,J,K);		
	clear I J K	
	
	disp('Solving...');	
	% Some variables
	two_alpha = 2*alpha;
	
	% Initialisation
	zm = zinit;%z(indices_nmask)=NaN;
	d = zeros(nrows,ncols,m);
	e = zeros(nrows,ncols,m);
	
	residual_mat = zeros(nrows,ncols,m);
    reskfull=zeros(nrows,ncols); %TODO vectorise before, at main file
    
    s_breg=s; 
	disp('Cholesky factorization...');	
    R = chol(A); 
	%%% Split-Bregman
	for it_breg = 1:maxit
		z_before = zm;	
		%%% z update
		% Compute System update		
		B(:) = 2*lambda*z0/alpha;       
        
        sum_b1s = sum(b1.*s_breg,3);
        sum_b2s = sum(b2.*s_breg,3);
        %%
        for indx=1:4 
            values =2*(sum_b1s*b1f(indx) + sum_b2s*b2f(indx));
            B(I_xy{indx}) = B(I_xy{indx})+values(indices_centre{indx});
            
            values = -2*sum_b1s*b1f(indx);
            B(I_x1_y{indx}) = B(I_x1_y{indx})+values(indices_centre{indx});
            
            values = -2*sum_b2s*b2f(indx);  
            B(I_x_y1{indx}) = B(I_x_y1{indx})+values(indices_centre{indx});  
        end 
        %% z update through L2 minimisation. use the cholesky factors     
        z1=R'\B;
        zm=R\z1;
% 		%%% Residual
        zx=Gx*zm;
        zy=Gy*zm;  
%         
        for k = 1:m %TODO do the vectorisation before
            b1k=b1(:,:,k);
            b1k=b1k(indices_mask);
            b2k=b2(:,:,k);
            b2k=b2k(indices_mask);
            sk=s(:,:,k);
            sk=sk(indices_mask);
            resk=zx.*b1k+zy.*b2k-sk;
            reskfull(indices_mask)=resk;
            residual_mat(:,:,k) =reskfull;
        end	
        
		%%% d update
		if(display>1)
			toc
			disp('Bregman variables updates')
			tic
		end		
		res_plus_e = (residual_mat+e);
		abs_res_plus_e = abs(res_plus_e);		
		d = res_plus_e.*max(abs_res_plus_e-two_alpha,0)./(abs_res_plus_e);
		d(abs_res_plus_e==0)=0;		
		%%% e update
		e = res_plus_e-d;			
			
		%%% CV test			
		residual = norm(z_before-zm)/norm(zm);
		if(display>0)
			fprintf(1,'Iteration %d : residual = %.08f\n',it_breg,residual);			
		end	
		
		if(residual<tol && it_breg>minit)
			break;
        end
        
         s_breg=s+d-e;

	end
	
end
