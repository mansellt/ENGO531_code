function calc_plot_ellipsoid(C,XYZ,S,n)
% to calculate and plot a confidence ellipsoid given the 3x3 point
% covariance matrix, C, the point coordinates, XYZ, the scale factor, S,
% and the number of faces for the ellipsoid, n. n=20 should suffice
% Derek Lichti 2023

[V,D]=eig(C);
abc=S*sqrt(diag(D));

% compute ellipsoid in standard orientation at the origin, rotate later
[X,Y,Z]=ellipsoid(0,0,0,abc(1),abc(2),abc(3),n);

xyz=[reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)];

[m,~]=size(xyz);
uvw=zeros(m,3);

% rotate the ellipsoid
for i=1:m
    uvw(i,:)=(V*xyz(i,:)')';
end
U=reshape(uvw(:,1),[n+1,n+1])+XYZ(1);
V=reshape(uvw(:,2),[n+1,n+1])+XYZ(2);
W=reshape(uvw(:,3),[n+1,n+1])+XYZ(3);

% plot
surf(U,V,W,ones(n+1,n+1)*max(abc)/S,'EdgeColor','none');
colormap turbo
