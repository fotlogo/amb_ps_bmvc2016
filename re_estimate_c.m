function [C_refined]=re_estimate_c(I,C,N,H,A,mask,shadow_threshold)  
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
	disp('C estimation (ambient version)...');
    [nrows,ncols,nb_images] = size(I);
    indices_mask = (mask>0);
%% calculate n.H for each image per pixel and rank images accordingly. This prety much determines photometric parallax
    NDotH=zeros(nrows,ncols,nb_images);   
    
    Cmap=zeros(nrows,ncols);    
	for yy = 1:nb_images          
        Hi_dot_N=zeros(nrows,ncols);
        for k=1:size(N,3) %i.e. 3           
            Nk=N(:,:,k);
            Hki=H(:,:,k,yy);
            Hi_dot_N(indices_mask)= Hi_dot_N(indices_mask)+ Hki(indices_mask).*Nk(indices_mask);
        end  
        
        Hi_dot_N(I(:,:,yy)<shadow_threshold)=-1;
        
        NDotH(:,:,yy)=Hi_dot_N;
    end
    %% sort according to photo parallax
     [~,ind] = sort(NDotH,3,'descend');    
    
    %% loop over pixles : BAD PERFORMANCE
    for xx=1: ncols   
        for yy=1: nrows   
            if( mask(yy,xx) ==0)
                continue;
            end  
            if( isnan(N(yy,xx,1)))
                continue;
            end  
            N_H = squeeze(NDotH(yy,xx,ind(yy,xx,:)));  %    squeeze(NDotH(ind(ii,jj,:)))                
            last_indx=nb_images; 
            %% Discard light sources where they are in shadow
            while( N_H(last_indx)< 3*shadow_threshold)
                last_indx=last_indx-1;      
                
                if(last_indx<2)
                    break;
                end
            end   
            
            if(last_indx<2)
                continue;
            end            
           
            nunmber_of_quadruples=floor( (last_indx+1)/4 );
            
            crev=zeros(nunmber_of_quadruples,1);
            
            for ii=1 : nunmber_of_quadruples
%                 ii
                jj= 2*ii;
                q = last_indx+1-2*ii;
                r = last_indx+1-ii;                
                
                indx_i=ind(yy,xx,ii);
                indx_j=ind(yy,xx,jj);
                indx_q=ind(yy,xx,q);  
                indx_r=ind(yy,xx,r);  
                
                dI_ij= I(yy,xx,indx_i)-I(yy,xx,indx_j);
                dI_qr= I(yy,xx,indx_q)-I(yy,xx,indx_r);
                
                %write down the equation as sum( Ki*Mi^x)=0 where x=1/c
                %K1=dI_qr*ai Mi=NDotH_i etc                
                K=zeros(4,1);
                M=zeros(4,1);
                
                K(1)= dI_qr.*A(yy,xx,indx_i);
                K(2)= -dI_qr.*A(yy,xx,indx_j);
                K(3)= -dI_ij.*A(yy,xx,indx_q);
                K(4)= dI_ij.*A(yy,xx,indx_r);
                
                M(1)=NDotH(yy,xx,indx_i);
                M(2)=NDotH(yy,xx,indx_j);
                M(3)=NDotH(yy,xx,indx_q);
                M(4)=NDotH(yy,xx,indx_r);
                
                crev(ii)=newton_solve(K,M,1./C(yy,xx));    
            end
           
            Cmap(yy,xx)=1/nanmean(crev);
        
        end
    end       
    %% JUST DONT UPDATE BAD VALUES
    indx=Cmap > 1.0;
    Cmap (indx)= C(indx);    
    indx=Cmap < 0.1;
    Cmap (indx)= C(indx);   
    indx=isnan(Cmap);
    Cmap (indx)= C(indx);      
    %% BLUR
	Cmap=imresize(Cmap,[nrows,ncols]/4); 
    Cmap=imresize(Cmap,[nrows,ncols]); 
   
    %% Keep a bit of the old as a prior
       C_refined=0.05*C+0.95*Cmap;      
       C_mean_new = nanmean(C_refined(indices_mask));
       
       fprintf(1, 'C mean new is %.4f \n',C_mean_new);       
%        C_refined=C_mean_new*ones(size(Cmap)); %could be helpfull
       
end
%solve for sum( Ki*Mi^x)=0 
%enforce that x>1
function xk=newton_solve(K,M,x0)

lamda=0; % regulariser not really needed

xk=x0;
for k=1:10   
    fx=K(1).*(M(1).^xk)+K(2).*(M(2).^xk)+K(3).*(M(3).^xk)+K(4).*(M(4).^xk) ;
    fx_reg=fx.^2+lamda*((xk-1).^2);
    fdotx=K(1).*log(M(1)).*(M(1).^xk)+K(2).*log(M(2)).*(M(2).^xk)+K(3).*log(M(3)).*(M(3).^xk)+K(4).*log(M(4)).*(M(4).^xk);
    fdotx_reg=2*fx.*fdotx+2*(xk-1);
%     xk=xk-fx./fdotx;
    xk=xk-fx_reg./fdotx_reg;
    if(abs(fx)<1e-16)
        break;
    end
    if(xk<1)
        break;
    end
end
    if(xk<1)
        xk=NaN;
    end
end
