function output=frequencyoffset(input,phase)
length_of_data=size(input,1);
frequency_offset=transpose(exp(1j*2*pi*phase*(0:(length_of_data-1))));
frequency_offset=repmat(frequency_offset,1,size(input,2));
output=input.*frequency_offset;
% data_train_high_cum=data_train(33:end,:);
% [label_train_high_cum,centers_train_high_cum]=find_centers(data_train_high_cum);
% [label_train,centers_train]=find_centers(data_train);
% [data_test,data_test_high_cum]= samples_testing(SNR);
% data_test=data_test.*frequency_offset;
end