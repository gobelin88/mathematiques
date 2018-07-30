pi=3.14159265359;

N = 1111;  
M = (N-1)/2;
K = 0.1;

fe=177556.0;
df=300;
fp = (150+df)/fe;
fs = (1000+df)/fe;

q  = [fp+K*(1-fs), fp*sinc(pi*fp*[1:2*M])-K*fs*sinc(pi*fs*[1:2*M])];

Q1 = toeplitz(q([0:M]+1));

C2=q([0:M]+1);
R2=q([M:2*M]+1);
M2 = size(C2,"*");
N2 = size(R2,"*");
COV2 = [matrix(C2,1,-1),matrix(R2(2:$),1,-1)];
Q2 = hank(M2,N2,COV2);

Q  = (Q1 + Q2)/2;

b  = fp*sinc(pi*fp*[0:M]');

a  = Q\b;

h  = [a(M+1:-1:2)/2; a(1); a(2:M+1)/2];

//plot(h);

[g,fr]=frmag(h,100000)
frt=fr(1:length(fr)/10);
gt=g(1:length(g)/10);
//plot2d(frt*fe,gt);
plot2d(frt*fe,20*log10(gt));

