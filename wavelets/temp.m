S=rand(1000000,1);
tic
N = length(S);
s1 = S(1:2:N-1) + sqrt(3)*S(2:2:N);
d1 = S(2:2:N) - sqrt(3)/4*s1 - (sqrt(3)-2)/4*[s1(N/2); s1(1:N/2-1)];
s2 = s1 - [d1(2:N/2); d1(1)];
s = (sqrt(3)-1)/sqrt(2) * s2;
d = -(sqrt(3)+1)/sqrt(2) * d1;
pic1 = [s;d];
toc
S=S';
tic
pic2 = wave_1d_multi(S,'db04',1);
toc