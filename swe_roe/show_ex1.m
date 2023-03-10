clear;close all;clc;

addpath('/Users/xudo627/developments/petsc/share/petsc/matlab/');

figure;
for i = 1 : 1000
    disp(i);
    filename = ['outputs/ex1_Nx_10_Ny_5_dt_0.010000_' num2str(i) '.dat'];
    h = PetscBinaryRead(filename);
    h = reshape(h,[3 50]);
    h = h(1,:);
    h = reshape(h,[10,5]);
    h = flipud(h');
    imagesc(h); caxis([0 10]); colorbar; colormap(jet);
    pause(0.1);
end