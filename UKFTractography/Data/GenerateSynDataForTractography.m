function GenerateSynDataForTractography()
paths;
u=icosahedron(2);
[S is_cross] = gen_2cross_w(u, 90, 1); %single tensor
%make slices
[nx ny dir]=size(S);
dir=dir/2;
S=S(:,:,1:dir);
u=u(1:dir,:);

S0=ones(nx,ny);
S=cat(3,S0,S);
S=reshape(S,[nx ny 1 dir+1]);
S=cat(3,S,S,S,S,S);
nrrdStrct.data=S;
nrrdStrct.gradients=[0 0 0;u];

mat2DWInhdr('single_tensor',nrrdStrct,'pnl');

clear S;
%single_tensor_free_water
[S is_cross] = gen_2cross_w(u, 90, 1,'fw'); %single tensor+fw
%make slices
[nx ny dir]=size(S);
S=cat(3,S0,S);
S=reshape(S,[nx ny 1 dir+1]);
S=cat(3,S,S,S,S,S);
nrrdStrct.data=S;
mat2DWInhdr('single_tensor_fw',nrrdStrct,'pnl');


clear S;
[S is_cross] = gen_2cross_w(u, 90, ukfHalf); %2-tensor
%make slices
[nx ny dir]=size(S);
S=cat(3,S0,S);
S=reshape(S,[nx ny 1 dir+1]);
S=cat(3,S,S,S,S,S);
nrrdStrct.data=S;
mat2DWInhdr('two_tensor',nrrdStrct,'pnl');


clear S;
[S is_cross] = gen_2cross_w(u, 90, ukfHalf,'fw'); %2-tensor+fw
%make slices
[nx ny dir]=size(S);
S=cat(3,S0,S);
S=reshape(S,[nx ny 1 dir+1]);
S=cat(3,S,S,S,S,S);
nrrdStrct.data=S;
mat2DWInhdr('two_tensor_fw',nrrdStrct,'pnl');


end
