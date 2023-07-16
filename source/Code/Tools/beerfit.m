function k=beerfit(x,y)
% BEERFIT(x,y) : ajustement selon une loi de beer de x a y:
%  y=a+(b-a).*exp(-k*x))
% En retour, on obtient les parametres [a, b, k] par ajustement des moindres carres
% ainsi que les statistiques (rms et R2, n)
ini=[0.3 0.5 1];
exitflag=0;
while exitflag~=1
    [k,feval,exitflag]=fminsearch('modexp',ini,[],x,y);
    ini=rand(1,3)*2-1;
end
ys=k(1,1)+(k(1,2)-k(1,1)).*exp(-k(1,3).*x); 
k(1,4)=std(y-ys);
k(1,5)=1-norm(y-ys)/norm(y);
k(1,6)=length(y);
t=[min(x):(max(x)-min(x))/100:max(x)];
plot(x,y,'o',t,k(1,1)+(k(1,2)-k(1,1)).*(exp(-k(1,3)*t)),'-')

