function [x0_new,X0_new,f_new,c_new] = Generate_undersampled_cube(x0,f,Ny,Nx,c,unds)

%Subsample the maps and reduce to desired size
x0_new = zeros(Ny,Nx,c/unds);
X0_new = zeros(Ny*Nx,c/unds);
f_new = zeros(c/unds,1);

counter = 1;
for i = round(linspace(1,c,c/unds))
    x0_new(:,:,counter)=imresize(x0(:,:,i), [Ny,Nx],'nearest');
    a = x0_new(:,:,counter);
    X0_new(:,counter) = a(:);
    f_new(counter) = f(i); 
    counter = counter + 1;
end
c_new = counter - 1;
