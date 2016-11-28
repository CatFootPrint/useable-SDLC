function output=SDLC(signal,param,signal_type)
%Initialization
dictionary =param.dictionary;
% step_1=param.step_1;
step_2=param.step_2;
lambda=param.lambda;
M=param.M;
% L=param.L;
dictionary = dictionary*diag(1./sqrt(sum(dictionary.*dictionary)));
coeff=param.coeff;
%calculation
% Psi=2*step_1*coeff-2/M*(dictionary'*dictionary-(1-lambda)*diag(ones(1,size(dictionary,2))))*coeff+2/M*dictionary'*signal;
% gap=Inf;
% threshold=0.01;
% error_ratio_sdlc=Inf;
coeff_ori=coeff;



m=M;
n=size(signal,2);
% CVX baseline
cvx_begin quiet
   variable coeff_vari(m,n)
   minimize(square_pos(norm(signal-dictionary*coeff_vari,'fro'))+lambda*square_pos(norm(coeff_vari,'fro'))+(1-lambda)*sum(sum(abs(coeff_vari))))
cvx_end
 
 % Custom method
%  coeff = prox_l1(dictionary*coeff, lambda)+prox_matrix(dictionary*coeff,1-lambda);
 
 % Comparison
 disp(['Error=',num2str(sum(sum((coeff_ori - coeff).^2)))]);









% while((gap>=threshold)||(error_ratio_sdlc>=1))
%     indicator_sparse=full(abs(coeff)>=1e-3);
%     step_1=norm(dictionary*dictionary')*10;
% Psi=coeff+1/M/step_1*indicator_sparse.*(dictionary'*(signal-dictionary*coeff))-(1-lambda)/M/step_1*coeff;    
%     coeffbuffer=coeff;
%     step_2=norm(full(coeff)*transpose(full(coeff)))*10;
%     dictionarybuffer=dictionary;
% coeff=Psi-indicator_sparse.*(1/4/M/step_1*sign(Psi-lambda/M/2/step_1*(ones(size(Psi))))-1/2*sign(Psi+lambda/M/2/step_1*(ones(size(Psi)))));
for k=1:M
    indicator_use_atom=full(abs(coeff(k,:))>=1e-2);
%     indicator_use_atom=full(abs(coeff(k,:))>=0);
    dictionary(:,k)=dictionary(:,k)+1/M/step_2*((signal-dictionary*coeff).*repmat(indicator_use_atom,size(signal,1),1))*transpose(coeff(k,:));
end
% gap=sum(sum((abs(dictionary-dictionarybuffer)).^2))/size(dictionary,1)/size(dictionary,2)/mean(mean(dictionarybuffer))+sum(sum(abs(coeff-coeffbuffer)).^2)/size(coeffbuffer,1)/size(coeffbuffer,2)/mean(mean(coeffbuffer));
[~,label_sdlc]=max(abs(coeff));
label=repmat(1:size(signal_type,2),1,size(label_sdlc,2)/size(signal_type,2));
error_ratio_sdlc=sum(sum(label_sdlc~=label))/size(label_sdlc,2);
% end
output.dictionary=dictionary;
output.coeff=coeff;
output.error_ratio=error_ratio_sdlc;
end








