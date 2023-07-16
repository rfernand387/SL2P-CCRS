function [u] =  uar(dbase)
uarr = [0.05,0.1,0.1,0.20,0.2,0.2];
uara = [ 0.025,.05,.05,.5,20,.02];
u.Albedo = (sum(((abs((dbase.Albedo.Estime - dbase.Albedo.Valid)./dbase.Albedo.Valid)<uarr(1))+(abs((dbase.Albedo.Estime - dbase.Albedo.Valid))<uara(1)))>0))/length(dbase.Albedo.Valid);
u.FAPAR = (sum(((abs((dbase.FAPAR.Estime - dbase.FAPAR.Valid)./dbase.FAPAR.Valid)<uarr(2))+(abs((dbase.FAPAR.Estime - dbase.FAPAR.Valid))<uara(2)))>0))/length(dbase.FAPAR.Valid);
u.FCOVER = (sum(((abs((dbase.FCOVER.Estime - dbase.FCOVER.Valid)./dbase.FCOVER.Valid)<uarr(3))+(abs((dbase.FCOVER.Estime - dbase.FCOVER.Valid))<uara(3)))>0))/length(dbase.FAPAR.Valid);
u.LAI = (sum(((abs((dbase.LAI.Estime - dbase.LAI.Valid)./dbase.LAI.Valid)<uarr(4))+(abs((dbase.LAI.Estime - dbase.LAI.Valid))<uara(4)))>0))/length(dbase.LAI.Valid);
u.LAI_Cab = (sum(((abs((dbase.LAI_Cab.Estime - dbase.LAI_Cab.Valid)./dbase.LAI_Cab.Valid)<uarr(5))+(abs((dbase.LAI_Cab.Estime - dbase.LAI_Cab.Valid))<uara(5)))>0))/length(dbase.LAI_Cab.Valid);
u.LAI_Cw = (sum(((abs((dbase.LAI_Cw.Estime - dbase.LAI_Cw.Valid)./dbase.LAI_Cw.Valid)<uarr(6))+(abs((dbase.LAI_Cw.Estime - dbase.LAI_Cw.Valid))<uara(6)))>0))/length(dbase.LAI_Cw.Valid);
return

