Y1 = raylrnd(5,1e4,1); 
histogram(Y1,'FaceColor','g')
hold on
%% 
Y2 = raylrnd(4,1e4,1); 
hold on
histogram(Y2,'FaceColor','r')
%%
Y3 = raylrnd(3,1e4,1); 
hold on
histogram(Y3,'FaceColor','b')
%%
Y4 = raylrnd(2,1e4,1); 
hold on
histogram(Y4,'FaceColor','y')
%%
Y5 = raylrnd(1,1e4,1); 
hold on
histogram(Y5,'FaceColor','c')
legend('\sigma=1','\sigma=2','\sigma=3','\sigma=4','\sigma=5')