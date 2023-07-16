% ============================================================================= =
% prosail
% Lib_PROSAIL.R
% ============================================================================= =
% PROGRAMMERS:
% Jean-Baptiste FERET <jb.feret@teledetection.fr>
% Florian de BOISSIEU <fdeboiss@gmail.com>
% Copyright 2019/11 Jean-Baptiste FERET
% ============================================================================= =
% This Library includes functions dedicated to PROSAIL simulation
% SAIL versions available are 4SAIL and 4SAIL2
% ============================================================================= =




% ' Performs PROSAIL simulation based on a set of combinations of input parameters
% ' @param Spec_Sensor list. Includes optical constants required for PROSPECT
% ' refractive index, specific absorption coefficients and corresponding spectral bands
% ' @param Input_PROSPECT  list. PROSPECT input variables
% ' @param N numeric. Leaf structure parameter
% ' @param CHL numeric. Chlorophyll content (microg.cm-2)
% ' @param CAR numeric. Carotenoid content (microg.cm-2)
% ' @param ANT numeric. Anthocyain content (microg.cm-2)
% ' @param BROWN numeric. Brown pigment content (Arbitrary units)
% ' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
% ' @param LMA numeric. Leaf Mass per Area (g.cm-2)
% ' @param PROT numeric. protein content  (g.cm-2)
% ' @param CBC numeric. NonProtCarbon-based constituent content (g.cm-2)
% ' @param alpha numeric. Solid angle for incident light at surface of leaf (simulation of roughness)
% ' @param TypeLidf numeric. Type of leaf inclination distribution function
% ' @param LIDFa numeric.
% ' if TypeLidf ==1, controls the average leaf slope
% ' if TypeLidf ==2, controls the average leaf angle
% ' @param LIDFb numeric.
% ' if TypeLidf ==1, controls the distribution's bimodality
% ' if TypeLidf ==2, unused
% ' @param lai numeric. Leaf Area Index
% ' @param q numeric. Hot Spot parameter
% ' @param tts numeric. Sun zeith angle
% ' @param tto numeric. Observer zeith angle
% ' @param psi numeric. Azimuth Sun / Observer
% ' @param rsoil numeric. Soil reflectance
% ' @param fraction_brown numeric. Fraction of brown leaf area
% ' @param diss numeric. Layer dissociation factor
% ' @param Cv numeric. vertical crown cover percentage
% ' = % ground area covered with crowns as seen from nadir direction
% ' @param Zeta numeric. Tree shape factor
% ' = ratio of crown diameter to crown height
% ' @param SAILversion character. choose between 4SAIL and 4SAIL2
% ' @param BrownVegetation list. Defines optical properties for brown vegetation, if not NULL
% ' - WVL, Reflectance, Transmittance
% ' - Set to NULL if use PROSPECT to generate it
% '
% ' @return list. rdot,rsot,rddt,rsdt
% ' rdot: hemispherical-directional reflectance factor in viewing direction
% ' rsot: bi-directional reflectance factor
% ' rsdt: directional-hemispherical reflectance factor for solar incident flux
% ' rddt: bi-hemispherical reflectance factor
% ' @import prospect
% ' @export
function Ref  = PRO4SAIL2(Spec_Sensor,Input_PROSPECT,N,CHL, ...
                     CAR,ANT,BROWN,EWT, ...
                     LMA,PROT,CBC,alpha, ...
                     TypeLidf,LIDFa,LIDFb,lai, ...
                     q,tts,tto,psi,rsoil, ...
                     fraction_brown, diss, Cv,Zeta, ...
                     SAILversion,BrownVegetation)

  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  % 	LEAF OPTICAL PROPERTIES	% % 
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  if (isempty(Input_PROSPECT))
    Input_PROSPECT = struct('CHL', CHL, 'CAR', CAR, 'ANT', ANT, 'BROWN', BROWN, 'EWT', EWT, ...
                                'LMA', LMA, 'PROT', PROT, 'CBC', CBC, 'N', N, 'alpha', alpha);
  end

%   GreenVegetation = PROSPECT(Spec_Sensor, ...
%                              Input_PROSPECT.N(1), ...
%                              Input_PROSPECT.CHL(1), ...
%                              Input_PROSPECT.CAR(1),  ...
%                              Input_PROSPECT.ANT(1), ...
%                              Input_PROSPECT.BROWN(1), ...
%                              Input_PROSPECT.EWT(1), ...
%                              Input_PROSPECT.LMA(1), ...
%                              Input_PROSPECT.PROT(1), ...
%                              Input_PROSPECT.CBC(1), ...
%                              Input_PROSPECT.alpha(1));
  GreenVegetation = prospect_DB(Input_PROSPECT.N(1), ...
                             Input_PROSPECT.CHL(1), ...
                             Input_PROSPECT.CAR(1),  ...
                             Input_PROSPECT.ANT(1), ...
                             Input_PROSPECT.BROWN(1), ...
                             Input_PROSPECT.EWT(1), ...
                             Input_PROSPECT.LMA(1));
  if (strcmp(SAILversion, '4SAIL2'))
    % 4SAIL2 requires one of the following combination of input parameters
    % Case % 1: valid optical properties for brown vegetation
    if (~isempty(BrownVegetation))
      % need to define Reflectance and Transmittance for BrownVegetation
      if (isempty(BrownVegetation.Reflectance)| isempty(BrownVegetation.Transmittance))
        print('Please define BrownVegetation as a list including "Reflectance" and "Transmittance"')
        return
      end
      % check if spectral domain for optical properties of brown vegetation match
      % with spectral domain for optical properties of green vegetation
      if (sum(Spec_Sensor.lambda-BrownVegetation.wvl)>0)
        print('Please define same spectral domain for BrownVegetation and SpecPROSPECT')
        return
      end
      % if all PROSPECT input parameters have the same length
      if (length(unique(structfun(@numel,Input_PROSPECT)))==1)
        % if NOT ALL OF PROSPECT input parameters are unique 
        if (~unique(structfun(@numel,Input_PROSPECT))==1)
          print('BrownVegetation defined along with multiple leaf chemical properties')
          print('Only first set of leaf chemical properties will be used to simulate green vegetation')
        end
      end
    % if no leaf optical properties brown vegetation defined
    elseif (isempty(BrownVegetation))
      % if all PROSPECT input parameters have the same length
      if (length(unique(structfun(@numel,Input_PROSPECT)))==1)
        % if all PROSPECT input parameters are unique (no possibility to simulate 2 types of leaf optics)
        if (unique(structfun(@numel,Input_PROSPECT))==1)
          % if fraction_brown set to 0, then assign green vegetation optics to brown vegetation optics
          if (fraction_brown==0)
            BrownVegetation = GreenVegetation
          % else run 4SAIL
          else 
            print('4SAIL2 needs two sets of optical properties for green and brown vegetation')
            print('Currently one set is defined. will run 4SAIL instead of 4SAIL2')
            SAILversion = '4SAIL'
          end
        % if all PROSPECT parameters have at least 2 elements
        elseif (unique(structfun(@numel,Input_PROSPECT))>=2)
          % compute leaf optical properties
%           BrownVegetation = PROSPECT(Spec_Sensor, ...
%                                         Input_PROSPECT.N(2), ...
%                                         Input_PROSPECT.CHL(2), ...
%                                         Input_PROSPECT.CAR(2), ...
%                                         Input_PROSPECT.ANT(2), ...
%                                         Input_PROSPECT.BROWN(2), ...
%                                         Input_PROSPECT.EWT(2), ...
%                                         Input_PROSPECT.LMA(2), ...
%                                         Input_PROSPECT.PROT(2), ...
%                                         Input_PROSPECT.CBC(2), ...
%                                         Input_PROSPECT.alpha(2))
         BrownVegetation = prospect_DB(Input_PROSPECT.N(2), ...
                                        Input_PROSPECT.CHL(2), ...
                                        Input_PROSPECT.CAR(2), ...
                                        Input_PROSPECT.ANT(2), ...
                                        Input_PROSPECT.BROWN(2), ...
                                        Input_PROSPECT.EWT(2), ...
                                        Input_PROSPECT.LMA(2));
          if (unique(structfun(@numel,Input_PROSPECT))>2)
            print('4SAIL2 needs two sets of optical properties for green and brown vegetation')
            print('Currently more than 2 sets are defined. will only use the first 2')
          end
        end
      end
    end
  end
  if (strcmp(SAILversion , '4SAIL'))
    if (length(unique(structfun(@numel,Input_PROSPECT)))==1)
      if (unique(structfun(@numel,Input_PROSPECT))>1)
        print('4SAIL needs only one set of optical properties')
        print('Currently more than one set of leaf chemical constituents is defined.')
        print('Will run 4SAIL with the first set of leaf chemical constituents')
      end
    end
  end

  if (strcmp(SAILversion , '4SAIL'))
    % run 4SAIL
    Ref = fourSAIL(GreenVegetation,TypeLidf, LIDFa, LIDFb, lai, q, tts, tto, psi, rsoil)
  elseif (strcmp(SAILversion , '4SAIL2'))
    % run 4SAIL2
    Ref = fourSAIL2(GreenVegetation, BrownVegetation,  ...
                     TypeLidf, LIDFa, LIDFb, lai, q, tts, tto, psi, rsoil,  ...
                     fraction_brown, diss, Cv, Zeta)
  end
  return

% ' Computes bidirectional reflectance factor based on outputs from PROSAIL and sun position
% '
% ' The direct and diffuse light are taken into account as proposed by:
% ' Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo
% ' measurements into total shortwave broadband albedo using a canopy
% ' radiative transfer model, Agronomie
% ' Es = direct
% ' Ed = diffuse
% '
% ' @param rdot numeric. Hemispherical-directional reflectance factor in viewing direction
% ' @param rsot numeric. Bi-directional reflectance factor
% ' @param tts numeric. Solar zenith angle
% ' @param SpecATM_Sensor list. direct and diffuse radiation for clear conditions
% ' @return BRF numeric. Bidirectional reflectance factor
% ' @export
function BRF = Compute_BRF(rdot,rsot,tts,SpecATM_Sensor)

  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  % % 	direct / diffuse light	% % 
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  Es = SpecATM_Sensor.Direct_Light;
  Ed = SpecATM_Sensor.Diffuse_Light;
  rd = pi/180;
  skyl = 0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd); % diffuse radiation (Francois et al.0, 2002)
  PARdiro = (1-skyl)*Es;
  PARdifo = skyl*Ed;
  BRF = (rdot*PARdifo+rsot*PARdiro)/(PARdiro+PARdifo);
  return 

% ' Performs PROSAIL simulation based on a set of combinations of input parameters
% ' @param LeafOptics list. Includes leaf optical properties (reflectance and transmittance)
% ' and corresponding spectral bands
% ' @param TypeLidf numeric. Type of leaf inclination distribution function
% ' @param LIDFa numeric.
% ' if TypeLidf ==1, controls the average leaf slope
% ' if TypeLidf ==2, controls the average leaf angle
% ' @param LIDFb numeric.
% ' if TypeLidf ==1, controls the distribution's bimodality
% ' if TypeLidf ==2, unused
% ' @param lai numeric. Leaf Area Index
% ' @param q numeric. Hot Spot parameter
% ' @param tts numeric. Sun zeith angle
% ' @param tto numeric. Observer zeith angle
% ' @param psi numeric. Azimuth Sun / Observer
% ' @param rsoil numeric. Soil reflectance
% '
% ' @return list. rdot,rsot,rddt,rsdt
% ' rdot: hemispherical-directional reflectance factor in viewing direction
% ' rsot: bi-directional reflectance factor
% ' rsdt: directional-hemispherical reflectance factor for solar incident flux
% ' rddt: bi-hemispherical reflectance factor
% ' @export

function my_list  = fourSAIL(LeafOptics, TypeLidf, LIDFa, LIDFb, lai, q, tts, tto, psi, rsoil)

  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  % 	LEAF OPTICAL PROPERTIES	% % 
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  rho = LeafOptics.Reflectance;
  tau = LeafOptics.Transmittance;

  % 	Geometric quantities
  rd = pi/180;
  cts = cos(rd*tts);
  cto = cos(rd*tto);
  ctscto = cts*cto;
  tants = tan(rd*tts);
  tanto = tan(rd*tto);
  cospsi = cos(rd*psi);
  dso = sqrt(tants*tants+tanto*tanto-2.*tants*tanto*cospsi);

  % 	Generate leaf angle distribution from average leaf angle (ellipsoidal) or (a,b) parameters
  if (TypeLidf==1)
    foliar_distrib = dladgen(LIDFa,LIDFb);
    lidf = foliar_distrib.lidf;
    litab = foliar_distrib.litab;
  elseif (TypeLidf==2)
    foliar_distrib = campbell(LIDFa);
    lidf = foliar_distrib.lidf;
    litab = foliar_distrib.litab;
  end

  % angular distance, compensation of shadow length
  % 	Calculate geometric factors associated with extinction and scattering
  % 	Initialise sums
  ks = 0;
  ko = 0;
  bf = 0;
  sob = 0;
  sof = 0;

  % 	Weighted sums over LIDF
  na = length(litab);
  for (i=1:na)
    ttl = litab(i);	    % leaf inclination discrete values
    ctl = cos(rd*ttl);
    % 	SAIL volume scattering phase function gives interception and portions to be
    % 	multiplied by rho and tau
    resVolscatt = volscatt(tts,tto,psi,ttl);
    chi_s = resVolscatt.chi_s;
    chi_o = resVolscatt.chi_o;
    frho = resVolscatt.frho;
    ftau = resVolscatt.ftau;

    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    % *                   SUITS SYSTEM COEFFICIENTS
    % *
    % *	ks  : Extinction coefficient for direct solar flux
    % *	ko  : Extinction coefficient for direct observed flux
    % *	att : Attenuation coefficient for diffuse flux
    % *	sigb : Backscattering coefficient of the diffuse downward flux
    % *	sigf : Forwardscattering coefficient of the diffuse upward flux
    % *	sf  : Scattering coefficient of the direct solar flux for downward diffuse flux
    % *	sb  : Scattering coefficient of the direct solar flux for upward diffuse flux
    % *	vf   : Scattering coefficient of upward diffuse flux in the observed direction
    % *	vb   : Scattering coefficient of downward diffuse flux in the observed direction
    % *	w   : Bidirectional scattering coefficient
    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*

    % 	Extinction coefficients
    ksli = chi_s/cts;
    koli = chi_o/cto;

    % 	Area scattering coefficient fractions
    sobli = frho*pi/ctscto;
    sofli = ftau*pi/ctscto;
    bfli = ctl*ctl;
    ks = ks+ksli*lidf(i);
    ko = ko+koli*lidf(i);
    bf = bf+bfli*lidf(i);
    sob = sob+sobli*lidf(i);
    sof = sof+sofli*lidf(i);
  end

  % 	Geometric factors to be used later with rho and tau
  sdb = 0.5*(ks+bf);
  sdf = 0.5*(ks-bf);
  dob = 0.5*(ko+bf);
  dof = 0.5*(ko-bf);
  ddb = 0.5*(1.+bf);
  ddf = 0.5*(1.-bf);

  % 	Here rho and tau come in
  sigb = ddb*rho+ddf*tau;
  sigf = ddf*rho+ddb*tau;
  att = 1-sigf;
  m2 = (att+sigb)*(att-sigb);
  m2(m2<=0)=0;
  m = sqrt(m2);

  sb = sdb*rho+sdf*tau;
  sf = sdf*rho+sdb*tau;
  vb = dob*rho+dof*tau;
  vf = dof*rho+dob*tau;
  w = sob*rho+sof*tau;

  % 	Here the LAI comes in
  %   Outputs for the case LAI = 0
  if (lai<0)
    tss = 1;
    too = 1;
    tsstoo = 1;
    rdd = 0;
    tdd = 1;
    rsd = 0;
    tsd = 0;
    rdo = 0;
    tdo = 0;
    rso = 0;
    rsos = 0;
    rsod = 0;

    rddt = rsoil;
    rsdt = rsoil;
    rdot = rsoil;
    rsodt = 0*rsoil;
    rsost = rsoil;
    rsot = rsoil;
  else 
    % 	Other cases (LAI > 0)
    e1 = exp(-m*lai);
    e2 = e1*e1;
    rinf = (att-m)/sigb;
    rinf2 = rinf*rinf;
    re = rinf*e1;
    denom = 1.-rinf2*e2;

    J1ks = Jfunc1(ks,m,lai);
    J2ks = Jfunc2(ks,m,lai);
    J1ko = Jfunc1(ko,m,lai);
    J2ko = Jfunc2(ko,m,lai);

    Ps = (sf+sb*rinf)*J1ks;
    Qs = (sf*rinf+sb)*J2ks;
    Pv = (vf+vb*rinf)*J1ko;
    Qv = (vf*rinf+vb)*J2ko;

    rdd = rinf*(1.-e2)/denom;
    tdd = (1.-rinf2)*e1/denom;
    tsd = (Ps-re*Qs)/denom;
    rsd = (Qs-re*Ps)/denom;
    tdo = (Pv-re*Qv)/denom;
    rdo = (Qv-re*Pv)/denom;

    tss = exp(-ks*lai);
    too = exp(-ko*lai);
    z = Jfunc3(ks,ko,lai);
    g1 = (z-J1ks*too)/(ko+m);
    g2 = (z-J1ko*tss)/(ks+m);

    Tv1 = (vf*rinf+vb)*g1;
    Tv2 = (vf+vb*rinf)*g2;
    T1 = Tv1*(sf+sb*rinf);
    T2 = Tv2*(sf*rinf+sb);
    T3 = (rdo*Qs+tdo*Ps)*rinf;

    % 	Multiple scattering contribution to bidirectional canopy reflectance
    rsod = (T1+T2-T3)/(1.-rinf2);

    % 	Treatment of the hotspot-effect
    alf = 1e6;
    % 	Apply correction 2/(K+k) suggested by F.-M. Breon
    if (q>0)
      alf = (dso/q)*2./(ks+ko);
    end
    if (alf>200)
      % inserted H. Bach 1/3/04
      alf = 200;
    end
    if (alf==0)
      % 	The pure hotspot - no shadow
      tsstoo = tss;
      sumint = (1-tss)/(ks*lai);
    else
      % 	Outside the hotspot
      fhot = lai*sqrt(ko*ks);
      % 	Integrate by exponential Simpson method in 20 steps
      % 	the steps are arranged according to equal partitioning
      % 	of the slope of the joint probability function
      x1 = 0;
      y1 = 0;
      f1 = 1;
      fint = (1.-exp(-alf))*0.05;
      sumint = 0;
      for (i =1:20)
        if (i<20)
          x2 = -log(1.-i*fint)/alf;
        else 
          x2 = 1;
        end
        y2 = -(ko+ks)*lai*x2+fhot*(1.-exp(-alf*x2))/alf;
        f2 = exp(y2);
        sumint = sumint+(f2-f1)*(x2-x1)/(y2-y1);
        x1 = x2;
        y1 = y2;
        f1 = f2;
      end
      tsstoo=f1;
    end
    % 	Bidirectional reflectance
    % 	Single scattering contribution
    rsos = w*lai*sumint;
    % 	Total canopy contribution
    rso = rsos+rsod;
    % 	Interaction with the soil
    dn = 1.-rsoil*rdd;
    % rddt: bi-hemispherical reflectance factor
    rddt = rdd+tdd*rsoil*tdd/dn;
    % rsdt: directional-hemispherical reflectance factor for solar incident flux
    rsdt = rsd+(tsd+tss)*rsoil*tdd/dn;
    % rdot: hemispherical-directional reflectance factor in viewing direction
    rdot = rdo+tdd*rsoil*(tdo+too)/dn;
    % rsot: bi-directional reflectance factor
    rsodt = rsod+((tss+tsd)*tdo+(tsd+tss*rsoil*rdd)*too)*rsoil/dn;
    rsost = rsos+tsstoo*rsoil;
    rsot = rsost+rsodt;
  end
  my_list = struct("rdot", rdot,"rsot", rsot,"rddt", rddt,"rsdt", rsdt);
  return

% ' Performs PRO4SAIL2 simulation based on a set of combinations of input parameters
% ' @param leafgreen list. includes reflectance and transmittance for vegetation % 1 (e.g. green vegetation)
% ' @param leafbrown list. includes reflectance and transmittance for vegetation % 2 (e.g. brown vegetation)
% ' @param TypeLidf numeric. Type of leaf inclination distribution function
% ' @param LIDFa numeric.
% ' if TypeLidf ==1, controls the average leaf slope
% ' if TypeLidf ==2, controls the average leaf angle
% ' @param LIDFb numeric.
% ' if TypeLidf ==1, unused
% ' if TypeLidf ==2, controls the distribution's bimodality
% ' @param lai numeric. Leaf Area Index
% ' @param hot numeric. Hot Spot parameter = ratio of the correlation length of leaf projections in the horizontal plane and the canopy height (doi:10.1016/j.rse.2006.12.013)
% ' @param tts numeric. Sun zeith angle
% ' @param tto numeric. Observer zeith angle
% ' @param psi numeric. Azimuth Sun / Observer
% ' @param rsoil numeric. Soil reflectance
% ' @param fraction_brown numeric. Fraction of brown leaf area
% ' @param diss numeric. Layer dissociation factor
% ' @param Cv numeric. vertical crown cover percentage
% ' = % ground area covered with crowns as seen from nadir direction
% ' @param Zeta numeric. Tree shape factor
% ' = ratio of crown diameter to crown height
% '
% ' @return list. rdot,rsot,rddt,rsdt
% ' rdot: hemispherical-directional reflectance factor in viewing direction
% ' rsot: bi-directional reflectance factor
% ' rsdt: directional-hemispherical reflectance factor for solar incident flux
% ' rddt: bi-hemispherical reflectance factor
% ' alfast: canopy absorptance for direct solar incident flux
% ' alfadt: canopy absorptance for hemispherical diffuse incident flux
% ' @export

function my_list = fourSAIL2(leafgreen, leafbrown, ...
                       TypeLidf,LIDFa,LIDFb, ...
                       lai, hot,tts,tto,psi,rsoil, ...
                       fraction_brown, diss, Cv,Zeta)

  % 	This version does not include non-Lambertian soil properties.
  % 	original codes do, and only need to add the following variables as input
  rddsoil = rsoil;
  rdosoil = rsoil;
  rsdsoil = rsoil;
  rsosoil = rsoil;

  % 	Geometric quantities
  rd = pi/180;

  % 	Generate leaf angle distribution from average leaf angle (ellipsoidal) or (a,b) parameters
  if (TypeLidf==1)
    foliar_distrib = dladgen(LIDFa,LIDFb);
    lidf = foliar_distrib.lidf;
    litab = foliar_distrib.litab;
  elseif (TypeLidf==2)
    foliar_distrib = campbell(LIDFa);
    lidf = foliar_distrib.lidf;
    litab = foliar_distrib.litab;
  end

  if (lai<0) 
    print('Please define positive LAI value')
    rddt = rsoil;
    rsdt = rsoil;
    rdot = rsoil;
    rsost = rsoil;
    rsot = rsoil;
    alfast =  0*rsoil;
    alfadt = 0*rsoil;
  elseif (lai==0)
    tss = 1.0;
    too = 1.0;
    tsstoo = 1.0;
    tdd = 1.0;
    rdd =  0.0;
    rsd = 0.0;
    tsd = 0.0;
    rdo =  0.0;
    tdo = 0.0;
    rso = 0.0;
    rsos = 0.0;
    rsod = 0.0;
    rsodt = 0.0;
    rddt = rsoil;
    rsdt = rsoil;
    rdot = rsoil;
    rsost = rsoil;
    rsot = rsoil;
    alfast =  0*rsoil;
    alfadt = 0*rsoil;
  elseif (lai>0)
    cts = cos(rd*tts);
    cto = cos(rd*tto);
    ctscto = cts*cto;
    tants = tan(rd*tts);
    tanto = tan(rd*tto);
    cospsi = cos(rd*psi);
    dso = sqrt(tants*tants+tanto*tanto-2.0*tants*tanto*cospsi);

    % Clumping effects
    Cs = 1.0;
    Co = 1.0;
    if (Cv<=1.0)
      Cs = 1.0-(1.0-Cv)^(1.0/cts);
      Co = 1.0-(1.0-Cv)^(1.0/cto);
    end
    Overlap = 0.0;
    if (Zeta>0.0)
      Overlap = min(Cs*(1.0-Co),Co*(1.0-Cs))*exp(-dso/Zeta);
    end
    Fcd = Cs*Co+Overlap;
    Fcs = (1.0-Cs)*Co-Overlap;
    Fod = Cs*(1.0-Co)-Overlap;
    Fos = (1.0-Cs)*(1.0-Co)+Overlap;
    Fcdc = 1.0-(1.0-Fcd)^(0.5/cts+0.5/cto);

    % 	Part depending on diss, fraction_brown, and leaf optical properties
    % 	First save the input fraction_brown as the old fraction_brown, as the following change is only artificial
    % Better define an fraction_brown that is actually used: fb, so that the input is not modified!

    fb = fraction_brown;
    % if only green leaves
    if (fraction_brown==0.0)
      fb = 0.5; % changed to below
      fb = 0.0001;
      leafbrown.Reflectance = leafgreen.Reflectance;
      leafbrown.Transmittance = leafgreen.Transmittance;
    end
    if (fraction_brown==1.0)
      fb = 0.5;
      leafgreen.Reflectance = leafbrown.Reflectance;
      leafgreen.Transmittance = leafbrown.Transmittance;
    end
    s = (1.0-diss)*fb*(1.0-fb);
    % rho1 & tau1 : green foliage
    % rho2 & tau2 : brown foliage (bottom layer)
    rho1 = ((1-fb-s)*leafgreen.Reflectance+s*leafbrown.Reflectance)/(1-fb);
    tau1 = ((1-fb-s)*leafgreen.Transmittance+s*leafbrown.Transmittance)/(1-fb);
    rho2 = (s*leafgreen.Reflectance+(fb-s)*leafbrown.Reflectance)/fb;
    tau2 = (s*leafgreen.Transmittance+(fb-s)*leafbrown.Transmittance)/fb;

    % angular distance, compensation of shadow length
    % 	Calculate geometric factors associated with extinction and scattering
    % 	Initialise sums
    ks = 0;
    ko = 0;
    bf = 0;
    sob = 0;
    sof = 0;

    % Weighted sums over LIDF

    for (i = 1:length(litab))
      ttl = litab(i);
      ctl = cos(rd*ttl);
      % SAIL volscatt function gives interception coefficients
      % and two portions of the volume scattering phase function to be
      % multiplied by rho and tau, respectively
      resVolscatt = volscatt(tts,tto,psi,ttl);
      chi_s = resVolscatt.chi_s;
      chi_o = resVolscatt.chi_o;
      frho = resVolscatt.frho;
      ftau = resVolscatt.ftau;
      % Extinction coefficients
      ksli = chi_s/cts;
      koli = chi_o/cto;
      % Area scattering coefficient fractions
      sobli = frho*pi/ctscto;
      sofli = ftau*pi/ctscto;
      bfli = ctl*ctl;
      ks = ks+ksli*lidf(i);
      ko = ko+koli*lidf(i);
      bf = bf+bfli*lidf(i);
      sob = sob+sobli*lidf(i);
      sof = sof+sofli*lidf(i);
    end
    % Geometric factors to be used later in combination with rho and tau
    sdb = 0.5*(ks+bf);
    sdf = 0.5*(ks-bf);
    dob = 0.5*(ko+bf);
    dof = 0.5*(ko-bf);
    ddb = 0.5*(1.+bf);
    ddf = 0.5*(1.-bf);

    % LAIs in two layers
    lai1 = (1-fb)*lai;
    lai2 = fb*lai;

    tss = exp(-ks*lai);
    ck = exp(-ks*lai1);
    alf = 1e6;
    if (hot>0.0)
      alf = (dso/hot)*2.0/(ks+ko);
    end
    if (alf>200.0)
      alf= 200.0;     % inserted H. Bach 1/3/04
    end
    if (alf==0.0)
      % The pure hotspot
      tsstoo = tss;
      s1 = (1-ck)/(ks*lai);
      s2 = (ck-tss)/(ks*lai);
    else 
      % Outside the hotspot
      fhot = lai*sqrt(ko*ks);
      % Integrate 2 layers by exponential simpson method in 20 steps
      % the steps are arranged according to equal partitioning
      % of the derivative of the joint probability function
      x1 = 0.0;
      y1 = 0.0;
      f1 = 1.0;
      ca = exp(alf*(fb-1.0));
      fint = (1.0-ca)*.05;
      s1 = 0.0;
      for istep = 1:20
        if (istep<20)
          x2 = -log(1.-istep*fint)/alf;
        else 
          x2 = 1.-fb;
        end
        y2 = -(ko+ks)*lai*x2+fhot*(1.0-exp(-alf*x2))/alf;
        f2 = exp(y2);
        s1 = s1+(f2-f1)*(x2-x1)/(y2-y1);
        x1 = x2;
        y1 = y2;
        f1 = f2;
      end
      fint = (ca-exp(-alf))*.05;
      s2 = 0.0;
      for (istep = 1:20)
        if (istep<20)
          x2 = -log(ca-istep*fint)/alf;
        else 
          x2 = 1.0;
        end
        y2 = -(ko+ks)*lai*x2+fhot*(1.0-exp(-alf*x2))/alf;
        f2 = exp(y2);
        s2 = s2+(f2-f1)*(x2-x1)/(y2-y1);
        x1 = x2;
        y1 = y2;
        f1 = f2;
      end
      tsstoo = f1
    end

    % Calculate reflectances and transmittances
    % Bottom layer
    tss = exp(-ks*lai2);
    too = exp(-ko*lai2);
    sb = sdb*rho2+sdf*tau2;
    sf = sdf*rho2+sdb*tau2;

    vb = dob*rho2+dof*tau2;
    vf = dof*rho2+dob*tau2;

    w2 = sob*rho2+sof*tau2;

    sigb = ddb*rho2+ddf*tau2;
    sigf = ddf*rho2+ddb*tau2;
    att = 1.0-sigf;
    m2 = (att+sigb).*(att-sigb);
    m2(m2<0) = 0;
    m = sqrt(m2);
    Which_NCS = find(m>0.01);
    Which_CS = find(m<=0.01);

    tdd = 0*m;
    rdd = 0*m;
    tsd = 0*m;
    rsd = 0*m;
    tdo = 0*m;
    rdo = 0*m;
    rsod = 0*m;
    if (~isempty(Which_NCS))
      resNCS = NonConservativeScattering(m(Which_NCS),lai2,att(Which_NCS),sigb(Which_NCS), ...
                                          ks,ko,sf(Which_NCS),sb(Which_NCS),vf(Which_NCS),vb(Which_NCS),tss,too);
      tdd(Which_NCS) = resNCS.tdd;
      rdd(Which_NCS) = resNCS.rdd;
      tsd(Which_NCS) = resNCS.tsd;
      rsd(Which_NCS) = resNCS.rsd;
      tdo(Which_NCS) = resNCS.tdo;
      rdo(Which_NCS) = resNCS.rdo;
      rsod(Which_NCS) = resNCS.rsod;
    end
    if (~isempty(Which_CS))
      resCS = ConservativeScattering(m(Which_CS),lai2,att(Which_CS),sigb(Which_CS), ...
                                      ks,ko,sf(Which_CS),sb(Which_CS),vf(Which_CS),vb(Which_CS),tss,too);
      tdd(Which_CS) = resCS.tdd;
      rdd(Which_CS) = resCS.rdd;
      tsd(Which_CS) = resCS.tsd;
      rsd(Which_CS) = resCS.rsd;
      tdo(Which_CS) = resCS.tdo;
      rdo(Which_CS) = resCS.rdo;
      rsod(Which_CS) = resCS.rsod;
    end

    % Set background properties equal to those of the bottom layer on a black soil
    rddb = rdd;
    rsdb = rsd;
    rdob = rdo;
    rsodb = rsod;
    tddb = tdd;
    tsdb = tsd;
    tdob = tdo;
    toob = too;
    tssb = tss;
    % Top layer
    tss = exp(-ks*lai1);
    too = exp(-ko*lai1);

    sb = sdb*rho1+sdf*tau1;
    sf = sdf*rho1+sdb*tau1;

    vb = dob*rho1+dof*tau1;
    vf = dof*rho1+dob*tau1;

    w1 = sob*rho1+sof*tau1;

    sigb = ddb*rho1+ddf*tau1;
    sigf = ddf*rho1+ddb*tau1;
    att = 1.0-sigf;

    m2 = (att+sigb).*(att-sigb);
    m2(m2<0) = 0;
    m = sqrt(m2);
    Which_NCS = find(m>0.01);
    Which_CS = find(m<=0.01);

    tdd = 0*m;
    rdd = 0*m;
    tsd = 0*m;
    rsd = 0*m;
    tdo = 0*m;
    rdo = 0*m;
    rsod = 0*m;
    if (~isempty(Which_NCS))
      resNCS = NonConservativeScattering(m(Which_NCS),lai1,att(Which_NCS),sigb(Which_NCS),  ...
                                          ks,ko,sf(Which_NCS),sb(Which_NCS),vf(Which_NCS),vb(Which_NCS),tss,too);
      tdd(Which_NCS) = resNCS.tdd;
      rdd(Which_NCS) = resNCS.rdd;
      tsd(Which_NCS) = resNCS.tsd;
      rsd(Which_NCS) = resNCS.rsd;
      tdo(Which_NCS) = resNCS.tdo;
      rdo(Which_NCS) = resNCS.rdo;
      rsod(Which_NCS) = resNCS.rsod;
    end
    if (~isempty(Which_CS))
      resCS = ConservativeScattering(m(Which_CS),lai1,att(Which_CS),sigb(Which_CS),  ...
                                      ks,ko,sf(Which_CS),sb(Which_CS),vf(Which_CS),vb(Which_CS),tss,too);
      tdd(Which_CS) = resCS.tdd;
      rdd(Which_CS) = resCS.rdd;
      tsd(Which_CS) = resCS.tsd;
      rsd(Which_CS) = resCS.rsd;
      tdo(Which_CS) = resCS.tdo;
      rdo(Which_CS) = resCS.rdo;
      rsod(Which_CS) = resCS.rsod;
    end

    % Combine with bottom layer reflectances and transmittances (adding method)
    rn = 1.0-rdd.*rddb;
    tup = (tss.*rsdb+tsd.*rddb)./rn;
    tdn = (tsd+tss.*rsdb.*rdd)./rn;
    rsdt = rsd+tup.*tdd;
    rdot = rdo+tdd.*(rddb.*tdo+rdob.*too)./rn;
    rsodt = rsod+(tss.*rsodb+tdn.*rdob).*too+tup.*tdo;

    rsost = (w1.*s1+w2.*s2)*lai;

    rsot = rsost+rsodt;

    % Diffuse reflectances at the top and the bottom are now different
    rddt_t = rdd+tdd.*rddb.*tdd./rn;
    rddt_b = rddb+tddb.*rdd.*tddb./rn;

    % Transmittances of the combined canopy layers
    tsst = tss.*tssb;
    toot = too.*toob;
    tsdt = tss.*tsdb+tdn.*tddb;
    tdot = tdob.*too+tddb.*(tdo+rdd.*rdob.*too)./rn;
    tddt = tdd.*tddb./rn;

    % Apply clumping effects to vegetation layer
    rddcb = Cv*rddt_b;
    rddct = Cv*rddt_t;
    tddc = 1-Cv+Cv*tddt;
    rsdc = Cs*rsdt;
    tsdc = Cs*tsdt;
    rdoc = Co*rdot;
    tdoc = Co*tdot;
    tssc = 1-Cs+Cs*tsst;
    tooc = 1-Co+Co*toot;

    % New weight function Fcdc for crown contribution (W. Verhoef, 22-05-08)
    rsoc = Fcdc*rsot;
    tssooc = Fcd*tsstoo+Fcs*toot+Fod*tsst+Fos;
    % Canopy absorptance for black background (W. Verhoef, 02-03-04)
    alfas = 1.-tssc-tsdc-rsdc;
    alfad = 1.-tddc-rddct;
    % Add the soil background
    rn = 1-rddcb.*rddsoil;
    tup = (tssc.*rsdsoil+tsdc.*rddsoil)./rn;
    tdn = (tsdc+tssc.*rsdsoil.*rddcb)./rn;

    rddt = rddct+tddc.*rddsoil.*tddc./rn;
    rsdt = rsdc+tup.*tddc;
    rdot = rdoc+tddc.*(rddsoil.*tdoc+rdosoil.*tooc)./rn;
    rsot = rsoc+tssooc.*rsosoil+tdn.*rdosoil.*tooc+tup.*tdoc;

    % Effect of soil background on canopy absorptances (W. Verhoef, 02-03-04)
    alfast = alfas+tup.*alfad;
    alfadt = alfad.*(1.+tddc.*rddsoil./rn);
  end
  my_list = struct("lambda",leafgreen.lambda,"rdot", rdot,"rsot" , rsot,"rddt" , rddt,"rsdt" , rsdt,  ...
                  "alfast" , alfast, "alfadt" , alfadt);
  return



% ' computes non conservative scattering conditions
% ' @param m numeric.
% ' @param lai numeric. Leaf Area Index
% ' @param att numeric.
% ' @param sigb numeric.
% ' @param ks numeric.
% ' @param ko numeric.
% ' @param sf numeric.
% ' @param sb numeric.
% ' @param vf numeric.
% ' @param vb numeric.
% ' @param tss numeric.
% ' @param too numeric.
% '
% ' @return list. tdd, rdd, tsd, rsd, tdo, rdo, rsod
% '
% ' @export
function my_list = NonConservativeScattering(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too)

  e1 = exp(-m*lai);
  e2 = e1.*e1;
  rinf = (att-m)./sigb;
  rinf2 = rinf.*rinf;
  re = rinf.*e1;
  denom = 1.0-rinf2.*e2;

  J1ks = Jfunc1(ks,m,lai);
  J2ks = Jfunc2(ks,m,lai);
  J1ko = Jfunc1(ko,m,lai);
  J2ko = Jfunc2(ko,m,lai);

  Ps = (sf+sb.*rinf).*J1ks;
  Qs = (sf.*rinf+sb).*J2ks;
  Pv = (vf+vb.*rinf).*J1ko;
  Qv = (vf.*rinf+vb).*J2ko;

  tdd = (1.0-rinf2).*e1./denom;
  rdd = rinf.*(1.0-e2)./denom;
  tsd = (Ps-re.*Qs)./denom;
  rsd = (Qs-re.*Ps)./denom;
  tdo = (Pv-re.*Qv)./denom;
  rdo = (Qv-re.*Pv)./denom;

  z = Jfunc2(ks,ko,lai);
  g1 = (z-J1ks.*too)./(ko+m);
  g2 = (z-J1ko.*tss)./(ks+m);

  Tv1 = (vf.*rinf+vb).*g1;
  Tv2 = (vf+vb.*rinf).*g2;

  T1 = Tv1.*(sf+sb.*rinf);
  T2 = Tv2.*(sf.*rinf+sb);
  T3 = (rdo.*Qs+tdo.*Ps).*rinf;

  % Multiple scattering contribution to bidirectional canopy reflectance
  rsod = (T1+T2-T3)./(1.0-rinf2);
  my_list = struct("tdd", tdd, "rdd", rdd, "tsd", tsd, ...
                  "rsd", rsd, "tdo", tdo, "rdo",  rdo, "rsod", rsod);
  return

% ' computes conservative scattering conditions
% ' @param m numeric.
% ' @param lai numeric. Leaf Area Index
% ' @param att numeric.
% ' @param sigb numeric.
% ' @param ks numeric.
% ' @param ko numeric.
% ' @param sf numeric.
% ' @param sb numeric.
% ' @param vf numeric.
% ' @param vb numeric.
% ' @param tss numeric.
% ' @param too numeric.
% '
% ' @return list. tdd, rdd, tsd, rsd, tdo, rdo, rsod
% '
% ' @export
function my_list = ConservativeScattering(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too)

  % Near or complete conservative scattering
  J4 = Jfunc4(m,lai);
  amsig = att-sigb;
  apsig = att+sigb;
  rtp = (1-amsig*J4)/(1+amsig*J4);
  rtm = (-1+apsig*J4)/(1+apsig*J4);
  rdd = 0.5*(rtp+rtm);
  tdd = 0.5*(rtp-rtm);

  dns = ks*ks-m*m;
  dno = ko*ko-m*m;
  cks = (sb*(ks-att)-sf*sigb)/dns;
  cko = (vb*(ko-att)-vf*sigb)/dno;
  dks = (-sf*(ks+att)-sb*sigb)/dns;
  dko = (-vf*(ko+att)-vb*sigb)/dno;
  ho = (sf*cko+sb*dko)/(ko+ks);

  rsd = cks*(1-tss*tdd)-dks*rdd;
  rdo = cko*(1-too*tdd)-dko*rdd;
  tsd = dks*(tss-tdd)-cks*tss*rdd;
  tdo = dko*(too-tdd)-cko*too*rdd;
  % Multiple scattering contribution to bidirectional canopy reflectance
  rsod = ho*(1-tss*too)-cko*tsd*too-dko*rsd;

  my_list = list("tdd", tdd, "rdd", rdd, "tsd", tsd, ...
                  "rsd", rsd, "tdo", tdo, "rdo",  rdo, "rsod", rsod);
  return







% ' Computes the leaf angle distribution function value (freq)
% '
% ' Using the original bimodal distribution function initially proposed in SAIL
% '  References
% '  ----------
% '  (Verhoef1998) Verhoef, Wout. Theory of radiative transfer models applied
% '  in optical remote sensing of vegetation canopies.
% '  Nationaal Lucht en Ruimtevaartlaboratorium, 1998.
% '  http://library.wur.nl/WebQuery/clc/945481.
% ' @param a controls the average leaf slope
% ' @param b controls the distribution's bimodality
% ' LIDF type 		  a 		b
% ' Planophile 	  1		  0
% ' Erectophile    -1	 	0
% ' Plagiophile 	  0		  -1
% ' Extremophile 	0		  1
% ' Spherical 	    -0.35 -0.15
% ' Uniform        0     0
% ' requirement: |LIDFa| + |LIDFb| < 1
% '
% ' @return foliar_distrib list. lidf and litab
% ' @export
function foliar_distrib = dladgen(a,b)
  litab=c(5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,81.0,83.0,85.0,87.0,89.0);
  freq=c();
  for (i1 = 1:8)
    t = i1*10;
    freq(i1) = dcum(a,b,t);
  end
  for (i2 = 9:12)
    t = 80.+(i2-8)*2.;
    freq(i2) = dcum(a,b,t);
  end
  freq(13) = 1;
  i=13;
  freq(i) = freq(i)-freq(i-1);
  foliar_distrib = struct("lidf", freq,"litab", litab);
  return

% ' dcum function
% ' @param a numeric. controls the average leaf slope
% ' @param b numeric. controls the distribution's bimodality
% ' @param t numeric. angle
% ' @return f
% ' @export
function f = dcum(a,b,t)
  rd = pi/180;
  if (a>=1)
    f = 1-cos(rd*t);
  else 
    eps = 1e-8;
    delx = 1;
    x = 2*rd*t;
    p = x;
    while (delx >= eps)
      y = a*sin(x)+.5*b*sin(2.*x);
      dx = .5*(y-x+p);
      x = x+dx;
      delx = abs(dx);
    end
    f = (2.*y+p)/pi;
  end
  return

% ' J1 function with avoidance of singularity problem
% '
% ' @param k numeric. Extinction coefficient for direct (solar or observer) flux
% ' @param l numeric.
% ' @param t numeric. Leaf Area Index
% ' @return Jout numeric.
% ' @export
function Jout = Jfunc1(k,l,t)
  % J1 function with avoidance of singularity problem
  del = (k-l)*t;
  Jout = 0*l;
  list1 = find(abs(del)>1e-3);
  list2 = find(abs(del)<=1e-3);

  Jout(list1) = (exp(-l(list1)*t)-exp(-k*t))./(k-l(list1));
  Jout(list2) = 0.5*t*(exp(-k*t)+exp(-l(list2)*t)).*(1-del(list2).*del(list2)/12);
  return


% ' J2 function with avoidance of singularity problem
% '
% ' @param k numeric. Extinction coefficient for direct (solar or observer) flux
% ' @param l numeric.
% ' @param t numeric. Leaf Area Index
% ' @return Jout numeric.
% ' @export
function Jout = Jfunc2(k,l,t)
  % 	J2 function
  Jout = (1.-exp(-(k+l)*t))./(k+l);
  return


% ' J3 function with avoidance of singularity problem
% '
% ' @param k numeric. Extinction coefficient for direct (solar or observer) flux
% ' @param l numeric.
% ' @param t numeric. Leaf Area Index
% ' @return Jout numeric.
% ' @export
function Jout = Jfunc3(k,l,t)
  out = (1.-exp(-(k+l)*t))./(k+l);
  return



% ' J4 function for treating (near) conservative scattering
% '
% ' @param m numeric. Extinction coefficient for direct (solar or observer) flux
% ' @param t numeric. Leaf Area Index
% ' @return Jout numeric.
% ' @export
function Jout = Jfunc4(m,t)

  del = m*t;
  out = 0*del;
  out(del>1e-3) = (1-exp(-del))/(m*(1+exp(-del)));
  out(del<=1e-3) = 0.5*t*(1.-del*del/12.);
  return



% ' Compute volume scattering functions and interception coefficients
% ' for given solar zenith, viewing zenith, azimuth and leaf inclination angle.
% '
% ' @param tts numeric. solar zenith
% ' @param tto numeric. viewing zenith
% ' @param psi numeric. azimuth
% ' @param ttl numeric. leaf inclination angle
% ' @return res list. includes chi_s, chi_o, frho, ftau
% ' @export
function res = volscatt(tts,tto,psi,ttl)
  % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  % *	chi_s	= interception functions
  % *	chi_o	= interception functions
  % *	frho	= function to be multiplied by leaf reflectance rho
  % *	ftau	= functions to be multiplied by leaf transmittance tau
  % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  % 	Wout Verhoef, april 2001, for CROMA

  rd = pi/180;
  costs = cos(rd*tts);
  costo = cos(rd*tto);
  sints = sin(rd*tts);
  sinto = sin(rd*tto);
  cospsi = cos(rd*psi);
  psir = rd*psi;
  costl = cos(rd*ttl);
  sintl = sin(rd*ttl);
  cs = costl*costs;
  co = costl*costo;
  ss = sintl*sints;
  so = sintl*sinto;

  % c ..............................................................................
  % c     betas -bts- and betao -bto- computation
  % c     Transition angles (beta) for solar (betas) and view (betao) directions
  % c     if thetav+thetal>pi/2, bottom side of the leaves is observed for leaf azimut
  % c     interval betao+phi<leaf azimut<2pi-betao+phi.
  % c     if thetav+thetal<pi/2, top side of the leaves is always observed, betao=pi
  % c     same consideration for solar direction to compute betas
  % c ..............................................................................

  cosbts = 5;
  if (abs(ss)>1e-6)
    cosbts = -cs/ss;
  end
  cosbto = 5;
  if (abs(so)>1e-6)
    cosbto = -co/so;
  end

  if (abs(cosbts)<1)
    bts = acos(cosbts);
    ds = ss;
  else
    bts = pi;
    ds = cs;
  end
  chi_s = 2./pi*((bts-pi*.5)*cs+sin(bts)*ss);
  if (abs(cosbto)<1)
    bto = acos(cosbto);
    doo = so;
  elseif(tto<90) 
    bto = pi;
    doo = co;
  else 
    bto = 0;
    doo = -co;
  end
  chi_o = 2./pi*((bto-pi*.5)*co+sin(bto)*so);

  % c ..............................................................................
  % c   Computation of auxiliary azimut angles bt1, bt2, bt3 used
  % c   for the computation of the bidirectional scattering coefficient w
  % c .............................................................................

  btran1 = abs(bts-bto);
  btran2 = pi-abs(bts+bto-pi);

  if (psir<=btran1)
    bt1 = psir;
    bt2 = btran1;
    bt3 = btran2;
  else 
    bt1 = btran1;
    if (psir<=btran2) 
      bt2 = psir;
      bt3 = btran2;
    else 
      bt2 = btran2;
      bt3 = psir;
    end
  end
  t1 = 2.*cs*co+ss*so*cospsi;
  t2 = 0;
  if (bt2>0) 
    t2 = sin(bt2)*(2.*ds*doo+ss*so*cos(bt1)*cos(bt3));
  end

  denom = 2.*pi*pi;
  frho = ((pi-bt2)*t1+t2)/denom;
  ftau = (-bt2*t1+t2)/denom;

  if (frho<0)
    frho = 0;
  end
  if (ftau<0)
    ftau = 0;
  end
  res = struct("chi_s", chi_s,"chi_o", chi_o,"frho", frho,"ftau", ftau);
  return
