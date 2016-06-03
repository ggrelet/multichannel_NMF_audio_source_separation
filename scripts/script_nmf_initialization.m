clear X2 X3 W H
%%
V = [1 2 3 ; 4 5 6];

%%
[W, H] = nmf_initialization(V, 2, 10^-7);
V2 = W*H