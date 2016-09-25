function [A,mapping_matrix,Omega] = make_Matrix_freeboundary_ICCV(b,lambda,mask,nrows,ncols)
	%%%
	% Build matrix A
	% In : B, lambda, nrows, ncols
	% Out : A, indices_centre
	%%%
	if (~exist('nrows','var')||isempty(nrows)), nrows=size(b,1);end
	if (~exist('ncols','var')||isempty(ncols)), ncols=size(b,2);end
	if (~exist('mask','var')||isempty(mask)), mask=ones(nrows,ncols);end

	% Create masks
% 	Omega = zeros(size(mask,1),size(mask,2));
	Omega_padded = padarray(mask,[1 1],0);
	
	Omega = Omega_padded(3:end,2:end-1).*Omega_padded(2:end-1,3:end).*mask;
	
	indices_mask = find(mask>0);
% 	indices_nmask = find(mask==0);
	nb_pixels_mask = length(indices_mask);
	mapping_matrix = zeros([nrows ncols]);
	mapping_matrix(indices_mask)=1:nb_pixels_mask; 

	
	% Initialize system
	I = [];
	J = [];
	K = [];
% 	B = zeros(nb_pixels_mask,1);
	
	% Regularization 	
	[X,Y]=find(mask>0);
	indices_centre=sub2ind(size(mask),X,Y);
	I_centre = mapping_matrix(indices_centre);
	A_centre = 2*lambda*ones(length(I_centre),1);	
	I=[I;I_centre(:)]; % Which lines are we talking about ? 
	J=[J;I_centre(:)]; % Which columns ?
	K=[K;A_centre(:)]; % Which values ? 	
	clear I_centre indices_centre
    
	% Terms in Omega(:,:,1)
	[X,Y]=find(Omega>0);
	indices_centre1=sub2ind(size(mask),X,Y);
	I_centre1 = mapping_matrix(indices_centre1);
	indices_voisins_x1=sub2ind(size(mask),X+1,Y);
	I_voisins_x1 = mapping_matrix(indices_voisins_x1);
	indices_voisins_y1=sub2ind(size(mask),X,Y+1);
	I_voisins_y1 = mapping_matrix(indices_voisins_y1);

	
	% Line of the central pixels
	values = 2*sum((b(:,:,:,1)+b(:,:,:,2)).^2,3);
	K = [K;values(indices_centre1)];
	I=[I;I_centre1(:)];
	J=[J;I_centre1(:)];

	sum_b1_b2 = sum(b(:,:,:,1).*b(:,:,:,2),3);
	sum_b1_2 = sum(b(:,:,:,1).^2,3);
	sum_b2_2 = sum(b(:,:,:,2).^2,3);
	
	values = -2*sum_b1_2-2*sum_b1_b2;
	K = [K;values(indices_centre1)];
	I=[I;I_centre1(:)];
	J=[J;I_voisins_x1(:)];
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_x1(:)];
	J=[J;I_centre1(:)];
	values = -2*sum_b1_b2-2*sum_b2_2;
	K = [K;values(indices_centre1)];
	I=[I;I_centre1(:)];
	J=[J;I_voisins_y1(:)];
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_y1(:)];
	J=[J;I_centre1(:)];
	% Line of the x-neighbor pixels
	values = 2*sum_b1_2;
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_x1(:)];
	J=[J;I_voisins_x1(:)];
	values = 2*sum_b1_b2;
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_x1(:)];
	J=[J;I_voisins_y1(:)];
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_y1(:)];
	J=[J;I_voisins_x1(:)];
	% Line of the y-neighbor pixels	
	values = 2*sum_b2_2;
	K = [K;values(indices_centre1)];
	I=[I;I_voisins_y1(:)];
	J=[J;I_voisins_y1(:)];	
	
	% Construction de A : 
    %% FOT: HACK CAUSE NO SPARCE SINGLE       
    if isa(K, 'single')
       K=double(K);        
    end
	A = sparse(I,J,K); 
    
end % EOF
