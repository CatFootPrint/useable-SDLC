function cumulants = High_Cumulants(signal)
%Calculating the high order cumulant
M20=mean(signal.^2);
M21=mean(signal.*conj(signal));
M40=mean(signal.^4);
M41=mean(signal.^3.*conj(signal));
M42=mean((signal.^2).*(conj(signal).^2));
M60=mean(signal.^6);
M63=mean(signal.^3.*(conj(signal).^3));
M80=mean(signal.^8);
C21=M21;
C40=M40-3*M20.^2;
C41=M41-3*M21.*M20;
C42=M42-abs(M20).^2-2*M21.^2;
C60=M60-15*M40.*M20+30*M20.^3;
C63=M63-9*C42.*C21-6*C21.^3;
C80=M80-28*M60.*M20-35*M40.^2+420*M40.*M20.^2-630*M20.^4;
f0=abs(C40./C42);
f1=abs(C42./C21.^2)  ;        % �����һ���ۻ���?% ��������������f2 f3
f2=(abs(C63)).^2./(abs(C42)).^3;   %f2��f3
f3=abs(C40)./abs(C21).^2;
f4=abs(C41)./abs(C21).^2;
f5=(abs(C60)).^2./(abs(C42)).^3;
f6=abs(C80)./(abs(C42)).^2;
cumulants=[f0;f1;f2;f3;f4;f5;f6];
end