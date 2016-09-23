function B = make_SecondMember_freeboundary_ICCV(b,s,lambda,z0,mask,mapping_matrix,Omega)
	%%%
	% Build second member B
	% In : b, s, lambda
	% Out : B
	%%%
	[nrows,ncols,~,~] = size(b);
	if (~exist('mask','var')||isempty(mask)), mask=ones(nrows,ncols);end

% 	[X,Y]=find(mask>0);
	nb_pixels_mask = sum(mask(:));
	B = zeros(nb_pixels_mask,1);
	
	% Regularization 	
% 	[X,Y]=find(mask>0);
% 	indices_centre=sub2ind(size(mask),X,Y);
% 	I_centre = mapping_matrix(indices_centre);
% 	A_centre = 2*lambda*ones(length(I_centre),1);	
	
	B(:) = 2*lambda*z0(mask);	
	
	% Terms in Omega(:,:,1)
	[X,Y]=find(Omega>0);
	indices_centre1=sub2ind(size(mask),X,Y);
	I_centre1 = mapping_matrix(indices_centre1);
	indices_voisins_x1=sub2ind(size(mask),X+1,Y);
	I_voisins_x1 = mapping_matrix(indices_voisins_x1);
	indices_voisins_y1=sub2ind(size(mask),X,Y+1);
	I_voisins_y1 = mapping_matrix(indices_voisins_y1);
		
	% Terms in Omega(:,:,1)
	b1s = sum(b(:,:,:,1).*s,3);
	b2s = sum(b(:,:,:,2).*s,3);
	values = -2*(b1s+b2s);
	B(I_centre1) = B(I_centre1)+values(indices_centre1);
	values = 2*b1s;
	B(I_voisins_x1) = B(I_voisins_x1)+values(indices_centre1);
	values = 2*b2s;
	B(I_voisins_y1) = B(I_voisins_y1)+values(indices_centre1);
	
end % EOF
