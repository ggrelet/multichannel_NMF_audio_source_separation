clear V W H
%%
V = [1 2 3 ; 4 5 6];

%%
[W, H] = nmf_initialization(M, 2, 10^-7);
V2 = W*H