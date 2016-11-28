function output=SDLC_proxi_advanced(signal,param,signal_type)
%Initialization
dictionary =param.dictionary;
lambda1=param.lambda1;
lambda2=param.lambda2;
M=param.M;
L=param.L;
iteration_times=param.numIteration;
dictionary = dictionary*diag(1./sqrt(sum(dictionary.*dictionary)));
coeff=param.coeff;
opts.tol = 1e-4; opts.maxit = iteration_times;
    opts.D0 = dictionary; opts.Y0 = coeff;
    opts.yType = 0;
    [dictionary,coeff,output] = dl_proxi_advanced(signal,M,lambda1,lambda2,opts);
[~,label_sdlc]=max(abs(coeff));
label=repmat(1:size(signal_type,2),1,size(label_sdlc,2)/size(signal_type,2));
error_ratio_sdlc=sum(sum(label_sdlc~=label))/size(label_sdlc,2);
output.dictionary=dictionary;
output.coeff=coeff;
output.error_ratio=error_ratio_sdlc;
end








