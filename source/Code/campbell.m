% ' Computes the leaf angle distribution function value (freq)
% '
% ' Ellipsoidal distribution function characterised by the average leaf
% ' inclination angle in degree (ala)
% ' Campbell 1986
% ' @param ala average leaf angle
% ' @return foliar_distrib list. lidf and litab
% ' @export
function foliar_distrib = campbell(ala)

  tx1 = [10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,82.0,84.0,86.0,88.0,90.0];
  tx2 = [0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,82.0,84.0,86.0,88.0];
  litab = (tx2+tx1)/2;
  n = length(litab);
  tl1 = tx1*(pi/180);
  tl2 = tx2*(pi/180);
  excent = exp(-1.6184e-5*ala^3+2.1145e-3*ala^2-1.2390e-1*ala+3.2491);
  sum0 = 0;

  freq=zeros(n,1);
  for (i = 1:n)
    x1 = excent/(sqrt(1.+excent^2.*tan(tl1(i))^2));
    x2 = excent/(sqrt(1.+excent^2.*tan(tl2(i))^2));
    if (excent==1)
      freq(i) = abs(cos(tl1(i))-cos(tl2(i)));
    else 
      alpha = excent/sqrt(abs(1-excent^2));
      alpha2 = alpha^2;
      x12 = x1^2;
      x22 = x2^2;
      alpx1   = 0*alpha2;
      alpx2   = 0*alpha2;
      almx1   = 0*alpha2;
      almx2   = 0*alpha2;
      if (excent>1)
        alpx1 = sqrt(alpha2(excent>1)+x12(excent>1));
        alpx2(excent>1) = sqrt(alpha2(excent>1)+x22(excent>1));
        dum = x1*alpx1+alpha2*log(x1+alpx1);
        freq(i) = abs(dum-(x2*alpx2+alpha2*log(x2+alpx2)));
      else 
        almx1 = sqrt(alpha2-x12);
        almx2 = sqrt(alpha2-x22);
        dum = x1*almx1+alpha2*asin(x1/alpha);
        freq(i) = abs(dum-(x2*almx2+alpha2*asin(x2/alpha)));
      end
    end
  end
  sum0 = sum(freq);
  freq0 = freq/sum0;
  foliar_distrib = struct("lidf", freq0,"litab",litab);
  return
