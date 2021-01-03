function [result] = wave_trafo(vec,trafo)

vec_length = numel(vec);

if trafo=='haar'
    %HAAR
    scaling = [1 1]'/sqrt(2);
    wavelet = [1 -1]'/sqrt(2);
elseif trafo=='db04'
    %DB4
    scaling = [1+sqrt(3), 3+sqrt(3), 3-sqrt(3), 1-sqrt(3)]'/(4*sqrt(2));
    wavelet = [1-sqrt(3), -(3-sqrt(3)), 3+sqrt(3), -(1+sqrt(3))]'/(4*sqrt(2));
end

trafo_length = numel(scaling);
overshoot = trafo_length-2;

if overshoot~=0
    %%%cyclic extension
    vec_extended = [vec,vec(1:overshoot)];
    %%%mirror at endpoint
    %vec_extended = [vec(overshoot:-1:1),vec,vec(vec_length:-1:vec_length-overshoot+1)];
    %vec_length = vec_length + 2;
else
    vec_extended = vec;
end

result=zeros(1,vec_length);

for n=1:vec_length/2
    result(n) = vec_extended(2*n-1:2*n+trafo_length-2) * scaling;
    result(vec_length/2+n) =  vec_extended(2*n-1:2*n+trafo_length-2) * wavelet;
end

end