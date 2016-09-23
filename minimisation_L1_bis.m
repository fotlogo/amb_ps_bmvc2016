function [z,d,e,it_breg] = minimisation_L1_bis(b,s,mask,z0,tol,minit,maxit,lambda,order,alpha,nGS,display,zinit)
% Input : 
% b : nrows x ncols x m x 2 -- REQUIRED
% s : nrows x ncols x m -- REQUIRED
% (mask : nrows x ncols) -- DEFAULT : ones(nrows,ncols)
% (z0 : nrows x ncols) -- DEFAULT : zeros(nrows,ncols)
% (tol : 1x1) -- DEFAULT : 1e-5
% (maxit : 1x1) -- DEFAULT : 100
% (lambda : 1x1) -- DEFAULT : 1e-9
% (order : 1x1) -- DEFAULT : 1
% (alpha : 1x1) -- DEFAULT : 1e-1
% (nGS : 1x1) -- DEFAULT : 50
% (display : 1x1) -- DEFAULT : 0
% Output : 
% z : nrows x ncols
	
	% Check arguments
	if(nargin<2)
		disp('At least two arguments required')
		return;
	end
	[nrows,ncols,m,~] = size(b);
	if (~exist('mask','var')||isempty(mask)) 
    	mask=ones(nrows,ncols);
	end
	if (~exist('z0','var')||isempty(z0))
        z0=zeros(nrows,ncols);
	end
	if (~exist('tol','var')||isempty(tol)) 
        tol=1e-5;
    end
	if (~exist('maxit','var')||isempty(maxit))
        maxit=100;
    end
	if (~exist('lambda','var')||isempty(lambda))
        lambda=1e-9;
    end
	if (~exist('order','var')||isempty(order)) 
        order=1;
    end
	if (~exist('alpha','var')||isempty(alpha)) 
        alpha=1e-1;
    end
	if (~exist('display','var')||isempty(display))
        display=0;
    end
    if (~exist('zinit','var')||isempty(zinit)) 
        zinit=z0;
    end

	% Create masks
	disp('Setting variables');
	Omega = zeros(size(mask,1),size(mask,2),4);
	Omega_padded = padarray(mask,[1 1],0);
	if(order>1)
		Omega(:,:,1) = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
		Omega(:,:,2) = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
		Omega(:,:,3) = Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
		Omega(:,:,4) = Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
	elseif(order==1)
		Omega(:,:,1) = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
        NOmega(:,:,1) = ~Omega(:,:,1);
        Omega(:,:,2) = NOmega(:,:,1).*Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
        NOmega(:,:,2) = ~Omega(:,:,2);
        Omega(:,:,3) = NOmega(:,:,1).*NOmega(:,:,2).*Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,1:end-2).*mask;
        NOmega(:,:,3) = ~Omega(:,:,3);
        Omega(:,:,4) = NOmega(:,:,1).*NOmega(:,:,2).*NOmega(:,:,3).*Omega_padded(1:end-2,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
	end
	clear Omega_padded  NOmega
	Omega_rep = repmat(Omega,[1 1 1 m]);	
	
	indices_mask = find(mask>0);
	indices_nmask = mask==0;
	nb_pixels_mask = length(indices_mask);
	mapping_matrix = zeros([nrows ncols]);
	mapping_matrix(indices_mask)=1:nb_pixels_mask; 

	% Auxiliary variables
	gamma = zeros(nrows,ncols,m,4,2);
	gamma(:,:,:,1,1) = -b(:,:,:,1);
	gamma(:,:,:,2,1) = -b(:,:,:,1);
	gamma(:,:,:,3,1) = b(:,:,:,1);
	gamma(:,:,:,4,1) = b(:,:,:,1);
	gamma(:,:,:,1,2) = -b(:,:,:,2);
	gamma(:,:,:,2,2) = b(:,:,:,2);
	gamma(:,:,:,3,2) = b(:,:,:,2);
	gamma(:,:,:,4,2) = -b(:,:,:,2);	
	delta = gamma(:,:,:,:,1)+gamma(:,:,:,:,2);
	
	epsilon = zeros(nrows,ncols,4,9);
	epsilon(:,:,:,1) = -2*squeeze(sum(gamma(:,:,:,:,1).*delta(:,:,:,:),3));
	epsilon(:,:,:,2) = -2*squeeze(sum(gamma(:,:,:,:,2).*delta(:,:,:,:),3));
	epsilon(:,:,:,3) = 2*squeeze(sum(gamma(:,:,:,:,1).*gamma(:,:,:,:,2),3));
	epsilon(:,:,:,4) = 2*squeeze(sum(delta(:,:,:,:).^2,3));
	epsilon(:,:,:,5) = 2*squeeze(sum(gamma(:,:,:,:,1).^2,3));
	epsilon(:,:,:,6) = 2*squeeze(sum(gamma(:,:,:,:,2).^2,3));
	s_rep = s(:,:,:,[1 1 1 1]);	
	epsilon(:,:,:,7) = -2*squeeze(sum(delta(:,:,:,:).*s_rep,3));
	epsilon(:,:,:,8) = 2*squeeze(sum(gamma(:,:,:,:,1).*s_rep,3));
	epsilon(:,:,:,9) = 2*squeeze(sum(gamma(:,:,:,:,2).*s_rep,3));
	clear s_rep


	disp('Constructing system');
	% Initialize system
	I = [];
	J = [];
	K = [];
	B = zeros(nb_pixels_mask,1);
	
	% Regularization 	
	[X,Y]=find(mask>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);
	A_centre = 2*lambda*ones(length(I_centre),1)/alpha;	
	I=[I;I_centre(:)]; % Which lines are we talking about ? 
	J=[J;I_centre(:)]; % Which columns ?
	K=[K;A_centre(:)]; % Which values ? 	
	clear I_centre indices_centre
	
	% Terms in Omega(:,:,1)
	[X,Y]=find(Omega(:,:,1)>0);
	indices_centre1=sub2ind(size(mask),X,Y);
	I_centre1 = mapping_matrix(indices_centre1);
	indices_voisins_x1=sub2ind(size(mask),X+1,Y);
	I_voisins_x1 = mapping_matrix(indices_voisins_x1);
	indices_voisins_y1=sub2ind(size(mask),X,Y+1);
	I_voisins_y1 = mapping_matrix(indices_voisins_y1);
	% Line of the central pixels
	values = epsilon(:,:,1,4);
	K = [K;values(indices_centre1)];
	I=[I;I_centre1(:)];
	J=[J;I_centre1(:)];
	values = epsilon(:,:,1,1);
	K = [K;values(indices_centre1)];
	I=[I;I_centre1(:)];
	J=[J;I_voisins_x1(:)];
	values = epsilon(:,:,1,2);
	K = [K;values(indices_centre1)];
	I=[I;I_centre1(:)];
	J=[J;I_voisins_y1(:)];
	% Line of the x-neighbor pixels
	values = epsilon(:,:,1,1);
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_x1(:)];
	J=[J;I_centre1(:)];
	values = epsilon(:,:,1,5);
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_x1(:)];
	J=[J;I_voisins_x1(:)];
	values = epsilon(:,:,1,3);
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_x1(:)];
	J=[J;I_voisins_y1(:)];
	% Line of the y-neighbor pixels
	values = epsilon(:,:,1,2);
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_y1(:)];
	J=[J;I_centre1(:)];
	values = epsilon(:,:,1,3);
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_y1(:)];
	J=[J;I_voisins_x1(:)];
	values = epsilon(:,:,1,6);
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_y1(:)];
	J=[J;I_voisins_y1(:)];
	
	% Terms in Omega(:,:,2)
	[X,Y]=find(Omega(:,:,2)>0);
	indices_centre2=sub2ind(size(mask),X,Y);
	I_centre2 = mapping_matrix(indices_centre2);
	indices_voisins_x2=sub2ind(size(mask),X+1,Y);
	I_voisins_x2 = mapping_matrix(indices_voisins_x2);
	indices_voisins_y2=sub2ind(size(mask),X,Y-1);
	I_voisins_y2 = mapping_matrix(indices_voisins_y2);
	% Line of the central pixels
	values = epsilon(:,:,2,4);
	K = [K;values(indices_centre2)];
	I=[I;I_centre2(:)];
	J=[J;I_centre2(:)];
	values = epsilon(:,:,2,1);
	K = [K;values(indices_centre2)];
	I=[I;I_centre2(:)];
	J=[J;I_voisins_x2(:)];
	values = epsilon(:,:,2,2);
	K = [K;values(indices_centre2)];
	I=[I;I_centre2(:)];
	J=[J;I_voisins_y2(:)];
	% Line of the x-neighbor pixels
	values = epsilon(:,:,2,1);
	K = [K;values(indices_centre2)];
	I=[I;I_voisins_x2(:)];
	J=[J;I_centre2(:)];
	values = epsilon(:,:,2,5);
	K = [K;values(indices_centre2)];
	I=[I;I_voisins_x2(:)];
	J=[J;I_voisins_x2(:)];
	values = epsilon(:,:,2,3);
	K = [K;values(indices_centre2)];
	I=[I;I_voisins_x2(:)];
	J=[J;I_voisins_y2(:)];	
	% Line of the y-neighbor pixels
	values = epsilon(:,:,2,2);
	K = [K;values(indices_centre2)];
	I=[I;I_voisins_y2(:)];
	J=[J;I_centre2(:)];
	values = epsilon(:,:,2,3);
	K = [K;values(indices_centre2)];
	I=[I;I_voisins_y2(:)];
	J=[J;I_voisins_x2(:)];
	values = epsilon(:,:,2,6);
	K = [K;values(indices_centre2)];
	I=[I;I_voisins_y2(:)];
	J=[J;I_voisins_y2(:)];
	
	% Terms in Omega(:,:,3)
	[X,Y]=find(Omega(:,:,3)>0);
	indices_centre3=sub2ind(size(mask),X,Y);
	I_centre3 = mapping_matrix(indices_centre3);
	indices_voisins_x3=sub2ind(size(mask),X-1,Y);
	I_voisins_x3 = mapping_matrix(indices_voisins_x3);
	indices_voisins_y3=sub2ind(size(mask),X,Y-1);
	I_voisins_y3 = mapping_matrix(indices_voisins_y3);
	% Line of the central pixels
	values = epsilon(:,:,3,4);
	K = [K;values(indices_centre3)];
	I=[I;I_centre3(:)];
	J=[J;I_centre3(:)];
	values = epsilon(:,:,3,1);
	K = [K;values(indices_centre3)];
	I=[I;I_centre3(:)];
	J=[J;I_voisins_x3(:)];
	values = epsilon(:,:,3,2);
	K = [K;values(indices_centre3)];
	I=[I;I_centre3(:)];
	J=[J;I_voisins_y3(:)];
	% Line of the x-neighbor pixels
	values = epsilon(:,:,3,1);
	K = [K;values(indices_centre3)];
	I=[I;I_voisins_x3(:)];
	J=[J;I_centre3(:)];
	values = epsilon(:,:,3,5);
	K = [K;values(indices_centre3)];
	I=[I;I_voisins_x3(:)];
	J=[J;I_voisins_x3(:)];
	values = epsilon(:,:,3,3);
	K = [K;values(indices_centre3)];
	I=[I;I_voisins_x3(:)];
	J=[J;I_voisins_y3(:)];	
	% Line of the y-neighbor pixels
	values = epsilon(:,:,3,2);
	K = [K;values(indices_centre3)];
	I=[I;I_voisins_y3(:)];
	J=[J;I_centre3(:)];
	values = epsilon(:,:,3,3);
	K = [K;values(indices_centre3)];
	I=[I;I_voisins_y3(:)];
	J=[J;I_voisins_x3(:)];
	values = epsilon(:,:,3,6);
	K = [K;values(indices_centre3)];
	I=[I;I_voisins_y3(:)];
	J=[J;I_voisins_y3(:)];
	
	% Terms in Omega(:,:,4)
	[X,Y]=find(Omega(:,:,4)>0);
	indices_centre4=sub2ind(size(mask),X,Y);
	I_centre4 = mapping_matrix(indices_centre4);
	indices_voisins_x4=sub2ind(size(mask),X-1,Y);
	I_voisins_x4 = mapping_matrix(indices_voisins_x4);
	indices_voisins_y4=sub2ind(size(mask),X,Y+1);
	I_voisins_y4 = mapping_matrix(indices_voisins_y4);
	% Line of the central pixels
	values = epsilon(:,:,4,4);
	K = [K;values(indices_centre4)];
	I=[I;I_centre4(:)];
	J=[J;I_centre4(:)];
	values = epsilon(:,:,4,1);
	K = [K;values(indices_centre4)];
	I=[I;I_centre4(:)];
	J=[J;I_voisins_x4(:)];
	values = epsilon(:,:,4,2);
	K = [K;values(indices_centre4)];
	I=[I;I_centre4(:)];
	J=[J;I_voisins_y4(:)];
	% Line of the x-neighbor pixels
	values = epsilon(:,:,4,1);
	K = [K;values(indices_centre4)];
	I=[I;I_voisins_x4(:)];
	J=[J;I_centre4(:)];
	values = epsilon(:,:,4,5);
	K = [K;values(indices_centre4)];
	I=[I;I_voisins_x4(:)];
	J=[J;I_voisins_x4(:)];
	values = epsilon(:,:,4,3);
	K = [K;values(indices_centre4)];
	I=[I;I_voisins_x4(:)];
	J=[J;I_voisins_y4(:)];
	% Line of the y-neighbor pixels
	values = epsilon(:,:,4,2);
	K = [K;values(indices_centre4)];
	I=[I;I_voisins_y4(:)];
	J=[J;I_centre4(:)];
	values = epsilon(:,:,4,3);
	K = [K;values(indices_centre4)];
	I=[I;I_voisins_y4(:)];
	J=[J;I_voisins_x4(:)];
	values = epsilon(:,:,4,6);
	K = [K;values(indices_centre4)];
	I=[I;I_voisins_y4(:)];
	J=[J;I_voisins_y4(:)];
	
	clear values epsilon
	
	% Construction de A : 
	A = sparse(I,J,K);		
	clear I J K	
	
	disp('Solving...');
	
	% Some variables
	sum_Omega = sum(Omega,3);
	sum_Omega(sum_Omega==0)=Inf;
	%~ sum_Omega_rep = repmat(sum_Omega,[1 1 m]);
	sum_Omega_rep = sum_Omega(:,:,ones(m,1));
	deux_alpha = 2*alpha;
	
	% Initialisation
	z = zinit;z(indices_nmask)=NaN;
	d = zeros(nrows,ncols,m);
	e = zeros(nrows,ncols,m);
	zx = zeros(nrows,ncols,4);
	zy = zeros(nrows,ncols,4);	
	
	% Display
	if(display>0)
		if(display>2)
			h1 = figure();
			hdisp=imagesc(z);
			colormap jet
			axis equal
			axis off
			axis ij
			colorbar
			title('$z$')
			drawnow
			
			h1 = figure();
			hdisp2=imagesc(sum(e.^2,3));
			colormap jet
			axis equal
			axis off
			axis ij
			colorbar
			title('$||e||^2$','Interpreter','Latex')
			drawnow
			
			h3 = figure();
			hdisp3=imagesc(sum(d.^2,3));
			colormap jet
			axis equal
			axis off
			axis ij
			colorbar
			title('$||d||^2$','Interpreter','Latex')
			drawnow
			
			h4 = figure();
			hdisp4=imagesc((z-z));
			colormap jet
			axis equal
			axis off
			axis ij
			colorbar
			title('$z^{k+1}-z^k$','Interpreter','Latex')
			drawnow
        end
		if(display>3)
			rotate3d on
			h5 = figure();
			hdisp5=surf(fliplr(z));
			colormap gray
			shading flat
			axis equal
			axis off
			view(90,20)
			colorbar
			zoom(1.3)
			drawnow
        end
	end	
        
    residual = zeros(nrows,ncols,m);	
    	zxp = z([2:end end],:)-z;zxp(isnan(zxp))=0;
		zxm = z-z([1 1:end-1],:);zxm(isnan(zxm))=0;
		zx(:,:,1) = Omega(:,:,1).*zxp;
		zx(:,:,2) = Omega(:,:,2).*zxp;
		zx(:,:,3) = Omega(:,:,3).*zxm;
		zx(:,:,4) = Omega(:,:,4).*zxm;
		zyp = z(:,[2:end end])-z;zyp(isnan(zyp))=0;
		zym = z-z(:,[1 1:end-1]);zym(isnan(zym))=0;
		zy(:,:,1) = Omega(:,:,1).*zyp;
		zy(:,:,2) = Omega(:,:,2).*zym;
		zy(:,:,3) = Omega(:,:,3).*zym;
		zy(:,:,4) = Omega(:,:,4).*zyp;
				
		zx_rep = repmat(zx,[1 1 1 m]);
		zy_rep = repmat(zy,[1 1 1 m]);
		%residual(:)=0;        
       
        
		for k = 1:4
			%residual = residual + squeeze(Omega_rep(:,:,k,:)).*(squeeze(zx_rep(:,:,k,:)).*b(:,:,:,1)+squeeze(zy_rep(:,:,k,:)).*b(:,:,:,2)-s);	
            zx_rep_k=permute(zx_rep(:,:,k,:),[1,2,4,3]);
            b1=b(:,:,:,1);
            zy_rep_k=permute(zy_rep(:,:,k,:),[1,2,4,3]);
            b2=b(:,:,:,2);
            Omega_rep_k=permute(Omega_rep(:,:,k,:),[1,2,4,3]);            
            
            residual = residual + Omega_rep_k.*(zx_rep_k.*b1+zy_rep_k.*b2-s);
		end		
		residual = residual./sum_Omega_rep;		
		
		%%% d update
		if(display>1)
			toc
			disp('Bregman variables updates')
			tic
		end		
		res_plus_e = (residual+e);
		abs_res_plus_e = abs(res_plus_e);		
		d = res_plus_e.*max(abs_res_plus_e-deux_alpha,0)./(abs_res_plus_e);
		d(abs_res_plus_e==0)=0;		
		%%% e update
		e = res_plus_e-d;
				
    
		B(:) = -2*lambda*z0(indices_mask)/alpha;	
		s_rep = repmat(s+d-e,[1 1 1 4]);
		epsilon7 = -2*squeeze(sum(delta(:,:,:,:).*s_rep,3));
		epsilon8 = 2*squeeze(sum(gamma(:,:,:,:,1).*s_rep,3));
		epsilon9 = 2*squeeze(sum(gamma(:,:,:,:,2).*s_rep,3));	
		% Terms in Omega(:,:,1)
		values7 = epsilon7(:,:,1);
		values8 = epsilon8(:,:,1);
		values9 = epsilon9(:,:,1);
		B(I_centre1) = B(I_centre1)+values7(indices_centre1);
		B(I_voisins_x1) = B(I_voisins_x1)+values8(indices_centre1);
		B(I_voisins_y1) = B(I_voisins_y1)+values9(indices_centre1);
		% Terms in Omega(:,:,2)
		values7 = epsilon7(:,:,2);
		values8 = epsilon8(:,:,2);
		values9 = epsilon9(:,:,2);
		B(I_centre2) = B(I_centre2)+values7(indices_centre2);
		B(I_voisins_x2) = B(I_voisins_x2)+values8(indices_centre2);
		B(I_voisins_y2) = B(I_voisins_y2)+values9(indices_centre2);
		% Terms in Omega(:,:,3)
		values7 = epsilon7(:,:,3);
		values8 = epsilon8(:,:,3);
		values9 = epsilon9(:,:,3);
		B(I_centre3) = B(I_centre3)+values7(indices_centre3);
		B(I_voisins_x3) = B(I_voisins_x3)+values8(indices_centre3);
		B(I_voisins_y3) = B(I_voisins_y3)+values9(indices_centre3);
		% Terms in Omega(:,:,4)
		values7 = epsilon7(:,:,4);
		values8 = epsilon8(:,:,4);
		values9 = epsilon9(:,:,4);
		B(I_centre4) = B(I_centre4)+values7(indices_centre4);
		B(I_voisins_x4) = B(I_voisins_x4)+values8(indices_centre4);
		B(I_voisins_y4) = B(I_voisins_y4)+values9(indices_centre4);
		B = -B;
	
    x = z(indices_mask);
    %FOT
    if isa(x,'single')
        x=double(x);
    end
    r=B-A*x;
    trilA = tril(A);
    %~    

	%%% Split-Bregman
	for it_breg = 1:maxit
		z_before = z;	
		%%% z update
		% Compute systeme
		if(display>1)
			disp('System update...')
			tic
		end
		B(:) = -2*lambda*z0(indices_mask)/alpha;	
		s_rep = repmat(s+d-e,[1 1 1 4]);       
        
		epsilon7 = -2*permute(sum(delta.*s_rep,3),[1,2,4,3]);
		epsilon8 = 2*permute(sum(gamma(:,:,:,:,1).*s_rep,3),[1,2,4,3]);
		epsilon9 = 2*permute(sum(gamma(:,:,:,:,2).*s_rep,3),[1,2,4,3]);	
		% Terms in Omega(:,:,1)
		values7 = epsilon7(:,:,1);
		values8 = epsilon8(:,:,1);
		values9 = epsilon9(:,:,1);
		B(I_centre1) = B(I_centre1)+values7(indices_centre1);
		B(I_voisins_x1) = B(I_voisins_x1)+values8(indices_centre1);
		B(I_voisins_y1) = B(I_voisins_y1)+values9(indices_centre1);
		% Terms in Omega(:,:,2)
		values7 = epsilon7(:,:,2);
		values8 = epsilon8(:,:,2);
		values9 = epsilon9(:,:,2);
		B(I_centre2) = B(I_centre2)+values7(indices_centre2);
		B(I_voisins_x2) = B(I_voisins_x2)+values8(indices_centre2);
		B(I_voisins_y2) = B(I_voisins_y2)+values9(indices_centre2);
		% Terms in Omega(:,:,3)
		values7 = epsilon7(:,:,3);
		values8 = epsilon8(:,:,3);
		values9 = epsilon9(:,:,3);
		B(I_centre3) = B(I_centre3)+values7(indices_centre3);
		B(I_voisins_x3) = B(I_voisins_x3)+values8(indices_centre3);
		B(I_voisins_y3) = B(I_voisins_y3)+values9(indices_centre3);
		% Terms in Omega(:,:,4)
		values7 = epsilon7(:,:,4);
		values8 = epsilon8(:,:,4);
		values9 = epsilon9(:,:,4);
		B(I_centre4) = B(I_centre4)+values7(indices_centre4);
		B(I_voisins_x4) = B(I_voisins_x4)+values8(indices_centre4);
		B(I_voisins_y4) = B(I_voisins_y4)+values9(indices_centre4);
		B = -B;
		if(display>1)
			toc
			disp('Gauss-Seidel iterations...')
			tic
        end
        %% FOT matrix is symmetric positive definite, CG should work
%         for ii=1:10
%             rann=randn(size(B));
%             r= (rann')*A*rann;
%             fprintf(1, 'test positive def %f \n', r);
%         end
%         
  %% OR direct solve if we dont run out of memory (this was commented out, not sure why)      
		 z(indices_mask)=A\B;
% 		for iGS=1:nGS
% 			%%% Gauss-Seidel update of u					
% 			% EffectiveGauss-Seidel update
% 			x=x+trilA\r; 
% 			r=B-A*x;     
% 			z(indices_mask)=x;
% 		end				
		%%% Residual
		if(display>1)
			toc
			disp('Residual evaluation...')
			tic
		end
		zxp = z([2:end end],:)-z;zxp(isnan(zxp))=0;
		zxm = z-z([1 1:end-1],:);zxm(isnan(zxm))=0;    
        
		zx(:,:,1) = Omega(:,:,1).*zxp;
		zx(:,:,2) = Omega(:,:,2).*zxp;
		zx(:,:,3) = Omega(:,:,3).*zxm;
		zx(:,:,4) = Omega(:,:,4).*zxm;
		zyp = z(:,[2:end end])-z;zyp(isnan(zyp))=0;
		zym = z-z(:,[1 1:end-1]);zym(isnan(zym))=0;
		zy(:,:,1) = Omega(:,:,1).*zyp;
		zy(:,:,2) = Omega(:,:,2).*zym;
		zy(:,:,3) = Omega(:,:,3).*zym;
		zy(:,:,4) = Omega(:,:,4).*zyp;
				
		zx_rep = repmat(zx,[1 1 1 m]);
		zy_rep = repmat(zy,[1 1 1 m]);
		residual(:)=0;
		for k = 1:4
			%residual = residual + squeeze(Omega_rep(:,:,k,:)).*(squeeze(zx_rep(:,:,k,:)).*b(:,:,:,1)+squeeze(zy_rep(:,:,k,:)).*b(:,:,:,2)-s);		
            zx_rep_k=permute(zx_rep(:,:,k,:),[1,2,4,3]);
            b1=b(:,:,:,1);
            zy_rep_k=permute(zy_rep(:,:,k,:),[1,2,4,3]);
            b2=b(:,:,:,2);
            Omega_rep_k=permute(Omega_rep(:,:,k,:),[1,2,4,3]);            
            
            residual = residual + Omega_rep_k.*(zx_rep_k.*b1+zy_rep_k.*b2-s);
		end		
		residual = residual./sum_Omega_rep;		
		
		%%% d update
		if(display>1)
			toc
			disp('Bregman variables updates')
			tic
		end		
		res_plus_e = (residual+e);
		abs_res_plus_e = abs(res_plus_e);		
		d = res_plus_e.*max(abs_res_plus_e-deux_alpha,0)./(abs_res_plus_e);
		d(abs_res_plus_e==0)=0;		
		%%% e update
		e = res_plus_e-d;
				
		
		%%% Graphics update		
		if(display>2)	
			tic
			disp('Graphics update');
			set(hdisp,'CData',z)
			drawnow
			
			set(hdisp2,'CData',sum(e.^2,3))
			drawnow
			
			set(hdisp3,'CData',sum(d.^2,3))
			drawnow
			
			set(hdisp4,'CData',(z-z_before))
			drawnow
			toc
		end
		if(display>3)		
			tic
			disp('3D Graphics update');
			set(hdisp5,'ZData',z)
			set(hdisp5,'CData',1./(zx(:,:,1).^2+zy(:,:,1).^2+1))
			[az,el]=view;
			view(az+2,el);
			drawnow
			toc
		end
					
		%%% CV test
		if(display>0)
			tic
		end
		if(display>1)			
			disp('CV test')
		end		
		residual = norm(z_before(indices_mask)-z(indices_mask))/norm(z(indices_mask));
		if(display>0)
			fprintf(1,'Iteration %d : residual = %.08f\n',it_breg,residual);			
		end	
		if(display>1)
			toc
		end	
		if(residual<tol && it_breg>minit)
			break;
		end

	end
	
end
