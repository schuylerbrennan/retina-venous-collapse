function [timedata,spacedata] = ...
  DAsolver(Pc,Dc,Ac,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0in,Lengthconstant,...
  multiplier,capdensity,respswitches,krogh,formatvars,condrespsignalc,...
  tauc, indicator, rref,IOPin,Pvar)

%%% Appears to make a large grid for time data and space data that have
%%% various attributes we want to measure. After setting those up, we solve
%%% the ODEs using DAall, and fill out the data after those are solved.

%initial guess/conditions

  u = [Dc.DLA Ac.ALA Dc.DSA Ac.ASA];
  %  These tolerances match up with above components in u
  tolvec = [1e-7 1e-4  1e-7    1e-4];

  %dynamic method of finding activation and diameter; loop through 
  %small time interval first (for non-oscillatory solutions)
  vessel_fields = {'LA','SA'};
  vessel_fields_cap = {vessel_fields{:},'C'};
  vessel_fields_veins = {'LA','SA','C','SV','LV'};

  %  Name of variables we want to store.  These are the variables that we
  %  are integrating, plus time...none of the others.
  timenames.var = {'D','A'};
  %  Names of flow variables we want to store.
  timenames.flow = {'Q','P','tau', 'R'};
  %  Names of response variables we want to store.
  timenames.resp = {'shear','myo','meta','total'};
  %  Names of metabolic response related variables.
  timenames.meta = {'sat','conc','scr','CO2','cons','PO2'};

  %  A lot of these fields will end up being empty...they don't take up
  %  much memory so that is ok
  for vfn = vessel_fields_veins
    for sfn = {'var'}
      for ssfn = timenames.(sfn{1})
        timedata.(sfn{1}).(ssfn{1}).(vfn{1}) = [];
      end
    end
  end
  timedata.t = [];
  ut = [];

  %  Store these at beginning and end of simulation
  spacenames.grid = {'x'};
  spacenames.meta = {'O2Hbsat','ATPconc','x','CO2','qO2cons','PO2'};
  for vfn = vessel_fields_veins
    for sfn = fieldnames(spacenames)'
      for ssfn = spacenames.(sfn{1})
        spacedata.(sfn{1}).(ssfn{1}).(vfn{1}) = [];
      end
    end
  end

  %  Also keep track of response data

  struccount = 0;

  for NN = 1:Prmc.NNmax
    tvec = Prmc.tmax*[NN-1,NN]./Prmc.NNmax;

    %%% Repeatedly integrate ODEs (both Diameter and Activation) until we
    %%% reach an equilibrium
    [t,u] = ode15s(@DAall,tvec,u(end,:),...
      odeset('InitialStep',0.1,'NonNegative',ones(size(u(end,:))),...
      'MaxStep',1000,'RelTol',min(tolvec),'AbsTol',min(tolvec)),...
      Pc,Dc,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0in,Lengthconstant,multiplier,...
      respswitches,krogh,formatvars,condrespsignalc,tauc, indicator, ...
      rref, IOPin,Pvar);

    ind = 1;
    for vfn = vessel_fields
      for sfn = timenames.var
        timedata.var.(sfn{1}).(vfn{1}) = ...
          [timedata.var.(sfn{1}).(vfn{1});u(:,ind)];
        ind = ind+1;
      end
    end
    timedata.t = [timedata.t;t];
    ut = [ut;u];

    %  A stopping criterion.  I think there are nicer ones out there (e.g.
    %  stop when norm(u'./u) < tol but this is probably fine for now.
    avguall = mean(u,1);
    %  Matlab now (2016 on I believe) supports basic singleton expansion
    %  (previously bsxfun) for all objects with singleton dimensions...here
    %  u-avguall = u-ones(size(t),1)*avguall.
    stduall = mean(abs(u-avguall),1);

    %testing
    stduallquant = sum(stduall < tolvec);
    fprintf('NN: %d; stduall stuff: %d\n',NN,stduallquant);

    if (stduallquant == numel(avguall)) && (~Prmc.fulltime)
      break
    elseif NN == Prmc.NNmax
      warning(['M0in = %g may be oscillating, using averaging',...
        ' techniques.\n'],M0in);
      %  Assuming relatively periodic oscillations, find location of peaks
      %  (we assume peaks occur only once every period).  Not really sure
      %  which variable is best to use for this (i.e. which variable has
      %  largest and most regular oscillations)...we try diameter of mids for
      %  now.
      [~,locs] = findpeaks(timedata.var.D.MID);
      %  Use average over last several peaks (minimum of 10 and (number of
      %  peaks)/10, if available)
      locstart = max(1,numel(locs)-min(10,ceil(numel(locs)/10)));
      if locstart == numel(locs)
        fprintf(['Solution doesn''t seem to be oscillating!  ',...
          'No averaging done!']);
      else
        inds = locs(locstart):locs(end);
        for vfc = 1:numel(vessel_fields)
          sfn = vessel_fields{vfc};
          %  The main variables
          timedata.var.D.(sfn)(end) = mean(timedata.var.D.(sfn)(inds));
          timedata.var.A.(sfn)(end) = mean(timedata.var.A.(sfn)(inds));
        end
      end
    else
      avgupast = avguall;
    end 
  end

  for tc = 1:numel(timedata.t)    
    [tmp,timedatastruc(tc),spatialdatastruc(tc)] = ...
      DAall(timedata.t(tc),ut(tc,:),Pc,Dc,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0in,...
      Lengthconstant,multiplier,respswitches,krogh,...
      formatvars,condrespsignalc,tauc,indicator, rref,IOPin,Pvar);
  end




  for vfnc = 1:numel(vessel_fields_veins)
    vfn = vessel_fields_veins{vfnc};
    for sfn = spacenames.meta
      spacedata.meta(1).(sfn{1}).(vfn) = ...
        spatialdatastruc(1).grid.(sfn{1})(:,vfnc);
      spacedata.meta(2).(sfn{1}).(vfn) = ...
        spatialdatastruc(end).grid.(sfn{1})(:,vfnc);
    end
    spacedata.grid.x.(vfn) = spatialdatastruc(1).grid.x(:,vfnc);
  end

  %  Store time data in structures
  %  Store flow data (Q-flow,P-pressure,tau-shear)
  for sc = 1:numel(timedatastruc)
    for fn = setdiff(fieldnames(timenames),'var')'
      for fn2 = timenames.(fn{1})
        for fn3 = fieldnames(timedatastruc(sc).(fn{1}).(fn2{1}))'
          timedata.(fn{1}).(fn2{1}).(fn3{1})(sc) = ...
            timedatastruc(sc).(fn{1}).(fn2{1}).(fn3{1});
        end
      end
    end
  end


end