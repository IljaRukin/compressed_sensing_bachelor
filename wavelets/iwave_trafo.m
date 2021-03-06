function [result] = iwave_trafo(vec,trafo)

vec_length = numel(vec);

if trafo=='haar'
    %HAAR
    scaling = [1 1]'/sqrt(2);
    wavelet = [1 -1]'/sqrt(2);
elseif trafo=='db04'
    %DB4
    scaling = [3-sqrt(3), 3+sqrt(3), 1+sqrt(3), 1-sqrt(3)]'/(4*sqrt(2));
    wavelet = [1-sqrt(3), -(1+sqrt(3)), 3+sqrt(3), -(3-sqrt(3))]'/(4*sqrt(2));
end

vec_sorted = zeros(1,vec_length);
for k=1:vec_length/2
    vec_sorted(2*k-1) = vec(k);
    vec_sorted(2*k) = vec(k+vec_length/2);
end

trafo_length = numel(scaling);
overshoot = trafo_length-2;

if overshoot~=0
    %%%cyclic extension
    vec_extended = [vec_sorted(vec_length-overshoot+1:vec_length),vec_sorted]; %cyclic
    %%%mirror at endpoint
    %vec_extended = vec_sorted; %cyclic
    %vec_length = vec_length -2;
else
    vec_extended = vec_sorted;
end

result=zeros(vec_length,1);

for n=1:vec_length/2
    result(2*n-1) = vec_extended(2*n-1:2*n+trafo_length-2) * scaling;
    result(2*n) = vec_extended(2*n-1:2*n+trafo_length-2) * wavelet;
end

end