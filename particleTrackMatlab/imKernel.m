function kernel = imKernel(w)
% Objective: filter/clean out noise from image using the same kernel that
% E. Weeks used for single pixel noise

% (1/B)*exp(-1*(i^2 + j^2 + k^2)/(4*lambda^2)) %eqn 8, M.C. Jenkins and
% Egelhaaf

lambda = 1; %noise correlation of 1 pixel

a = -w:w;
Bpre = zeros([length(a),1],'single');
weight = zeros([length(a) length(a) length(a)],'single');

for i = 1:length(a)
    Bpre(i) = exp((-1*(a(i)^2))/(4*(lambda^2)));
    for j = 1:length(a)
        for k = 1:length(a)
            weight(i,j,k) = exp((-1*((a(i)^2)+(a(j)^2)+(a(k)^2)))/(4*lambda^2));
        end
    end
end
B = sum(Bpre(:));
kernel = B*weight;
end