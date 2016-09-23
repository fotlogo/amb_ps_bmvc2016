function [ A, H ] = calculate_A_H_fields(X,Y,Z,C,indices_mask,S,Sd,Phi,mu,epsil )
%calculate_A_H_fields Summary of this function goes here
% A  Light attenuation
% H field in fundamental equation
nb_images =size(Phi,1);
[nrows,ncols]=size(X);

Lbar = zeros(nrows,ncols,nb_images,3); % Unit direction
A = zeros(nrows,ncols,nb_images); % Attenuation 

% Initial view field
 V = zeros(nrows,ncols,3);

% Initial H field
 H = zeros(nrows,ncols,3,nb_images);
 
 L_k_1=NaN(nrows,ncols); % Non unit direction components
 L_k_2=NaN(nrows,ncols);
 L_k_3=NaN(nrows,ncols); 
 
 Lnorm_k=NaN(nrows,ncols);  % Distance to the source
 
 Lbar_k_i=NaN(nrows,ncols); 
 
 A_k=NaN(nrows,ncols); 
 
% View
distance = sqrt(X.^2+Y.^2+Z.^2);

V(:,:,1) = -X./distance;
V(:,:,2) = -Y./distance;
V(:,:,3) = -Z./distance;

for k = 1:nb_images   
    L_k_1(indices_mask)= S(1,k)-X(indices_mask);  
    L_k_2(indices_mask)= S(2,k)-Y(indices_mask);
    L_k_3(indices_mask)= S(3,k)-Z(indices_mask);

    Lnorm_k(indices_mask)=sqrt( L_k_1(indices_mask).^2 + L_k_2(indices_mask).^2 +L_k_3(indices_mask).^2 );    
%re-use Lbar_k_i for all 3 components
    Lbar_k_i(indices_mask)= L_k_1(indices_mask) ./ Lnorm_k(indices_mask);
    Lbar(:,:,k,1) =  Lbar_k_i;
     
    Lbar_k_i(indices_mask)= L_k_2(indices_mask) ./ Lnorm_k(indices_mask);
    Lbar(:,:,k,2) =  Lbar_k_i;
     
    Lbar_k_i(indices_mask)= L_k_3(indices_mask) ./ Lnorm_k(indices_mask);
    Lbar(:,:,k,3) =  Lbar_k_i;        
%assume sd normalised
    cosfi= Lbar(:,:,k,1).*Sd(1,k) + Lbar(:,:,k,2).*Sd(2,k) + Lbar(:,:,k,3).*Sd(3,k) ;    
    cosfi=max(-cosfi,0);
    % Attenuation 
    A_k(indices_mask)=Phi(k).*(cosfi(indices_mask).^mu(k))./(Lnorm_k(indices_mask).^2);      
	A(:,:,k) = A_k;      
    
    % H vector
    for i = 1:3
       Hki = H(:,:,i,k);
       Lki = Lbar(:,:,k,i);
       Vk = V(:,:,i);
       Hki(indices_mask) = Lki(indices_mask) + Vk(indices_mask).*min(1,abs(C(indices_mask)-1)/epsil);     
       H(:,:,i,k) = Hki;       
    end      
   H(:,:,:,k) = real(H(:,:,:,k)./repmat(sqrt(sum(H(:,:,:,k).^2,3)),[1 1 3 1]));
end
A=A/max(A(:)); %normalise for better numerical stability
end