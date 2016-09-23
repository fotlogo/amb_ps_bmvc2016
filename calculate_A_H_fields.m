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
    
%     figure;
%     title(sprintf('LNormK for img %d \n', k));
%     %surf(Lnorm_k);
%     mesh(log(Lnorm_k));
%     AA=Lnorm_k;
%     axis equal  
%     save('a.mat', 'AA');
%     pause
    
%re-use Lbar_k_i for all 3 components
    Lbar_k_i(indices_mask)= L_k_1(indices_mask) ./ Lnorm_k(indices_mask);
    Lbar(:,:,k,1) =  Lbar_k_i;
     
    Lbar_k_i(indices_mask)= L_k_2(indices_mask) ./ Lnorm_k(indices_mask);
    Lbar(:,:,k,2) =  Lbar_k_i;
     
    Lbar_k_i(indices_mask)= L_k_3(indices_mask) ./ Lnorm_k(indices_mask);
    Lbar(:,:,k,3) =  Lbar_k_i;        
%assume sd normalised
    cosfi= Lbar(:,:,k,1).*Sd(k,1) + Lbar(:,:,k,2).*Sd(k,2) + Lbar(:,:,k,3).*Sd(k,3) ;    
%     nanmean(cosfi(:))    
    cosfi=max(-cosfi,0);
    
    A_k(indices_mask)=Phi(k).*(cosfi(indices_mask).^mu(k))./(Lnorm_k(indices_mask).^2);     
	% Attenuation 
% 	A_k(indices_mask)=Phi(k).*((Z(indices_mask)-S(3,k)).^mu(k))./(Lnorm_k(indices_mask).^(mu(k)+2));     
    
	A(:,:,k) = A_k;      
    
    % H vector
    for i = 1:3
       Hki = H(:,:,i,k);
       Lki = Lbar(:,:,k,i);
       Vk = V(:,:,i);
       Hki(indices_mask) = Lki(indices_mask) + Vk(indices_mask).*min(1,abs(C(indices_mask)-1)/epsil);
       %Hki(indices_mask) = Lki(indices_mask); %LAMBERTIAN HALF VECTOR
       H(:,:,i,k) = Hki;       
    end  
   % Hki=mean(Hki(:))
    
   H(:,:,:,k) = real(H(:,:,:,k)./repmat(sqrt(sum(H(:,:,:,k).^2,3)),[1 1 3 1]));
end
A=A/max(A(:));
end