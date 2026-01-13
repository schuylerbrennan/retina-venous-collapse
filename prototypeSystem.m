function [timedata,spacedata,metadata] = prototypeSystem(...
  myoFEEDin,myoMIDin,myoCOLin,...
  shearFEEDin,shearMIDin,shearCOLin,...
  Cmeta,Length_constant,multiplier,capdensity)

  %Author: Myson Burch
  %Last Modified: 14 July 2016
  %Prototype model based on Arciero 2008 Blood Flow Regulation publications
  %Edit/cleanup: Jared Barber 3/17/2017ish
  %
  % Code description
  %  Currently this model is of the flow and diameter dynamics of the
  %  vasculature in a rat hind limb starting from the central iliac artery,
  %  branching into an external and internal iliac artery and then continues
  %  down through a COLateral artery, FEED arterioles, MID arterioles,
  %  capillaries and then back through veins.
  %
  % Inputs
  %  myoFEEDin-on or off:  myogenic response in FEED vessels
  %  myoMIDin-on or off:  myogenic response in MID vessels
  %  myoCOLin-on or off:  myogenic response in COLateral vessels
  %  SFEEDin-on or off:  shear response in FEED vessels
  %  SMIDin-on or off:  shear response in MID vessels
  %  SCOLin-on or off:  shear response in COLateral vessels
  %  Cmeta-default {30}:  magnitude of the metabolic response
  %  Length_constant-default {1cm}:  length constant associated with
  %    exponential decay of the conducted response upstream
  %  multiplier-default {1}:  as conducted response signal moves upstream it
  %    gets multiplied from 1 node into multiplier nodes along the way
  %    somehow?  Not yet used
  %  capdensity-default {500}:

%   close all
  %  Default values for inputs (if not fed in, these will be used,
  %  particularly helpful for debugging purposes)
  if nargin == 0
    myoFEEDin = true; myoMIDin = true; myoCOLin = true;
    shearFEEDin = true; shearMIDin = true; shearCOLin = true;
    Cmeta = 30; Length_constant = 3; multiplier = 1; capdensity = 500;
  end
  if multiplier > 1
    error('multiplier > 1 does not seem to be implemented yet');
  end
  %  Throw myogenic and shear response switches into a nice structure that we
  %  can pass back and forth.  Note that the control state uses all three
  %  responses while non-control states either allow the response levels
  %  to change (true) or always use the control state response levels
  %  (false; response levels constant in this case).
  respswitches.myo.FEED = myoFEEDin;
  respswitches.myo.MID = myoMIDin;
  respswitches.myo.COL = myoCOLin;
  respswitches.shear.FEED = shearFEEDin;
  respswitches.shear.MID = shearMIDin;
  respswitches.shear.COL = shearCOLin;
  %  Just in case we want to turn off the metabolic response later on.
  respswitches.meta.FEED = true; %metaFEEDin;
  respswitches.meta.MID = true; %metaMIDin;
  respswitches.meta.COL = true; %metaCOLin;

  %  Default values for nice pictures...shouldn't affect calculations
  formatvars = initializeformatvars;

  %  Define system control state values using experimental-based info.  Note
  %  that some of the values are known given the other values (e.g. given
  %  pressures (Pc) and lengths (Lc) of vessels, we can easily obtain the
  %  shear stress on the walls (tauc)).  Nonetheless, we store those values
  %  to make for easier computations later.
  [Dc,tauc,Pc,nc,Lc,muc,M0c,Prmc] = SetCState;

  %  Use previous values along with some constraints (e.g. activation
  %  levels are 0.5 in the MID and FEED arteries in the control state and
  %  0.99 in the COLateral arteries...changeable ContAct values inside of
  %  SetCState) to obtain other control state values (e.g. diameters of
  %  MID, FEED, and COLateral arteries)
  [Dc,tauc,Pc,Ac,nc,Lc,muc,M0c,Prmc,Rc,condrespsignalc,krogh] = ...
    GetCState(Cmeta,Length_constant,multiplier,...
    Pc,nc,Dc,Lc,tauc,muc,M0c,Prmc,formatvars,respswitches);

  %  Below we loop through a different set of occlusions (partial occlusion
  %  can now occur).  For each "occlusion level" we loop through multiple
  %  M0 values.  Might or might not want to change this later so that one
  %  can have different values for M0 depending on the occlusion level but
  %  right now that doesn't happen.  Also note that in the past we also
  %  looped through different Pa values.  Additional parameters can be
  %  added using the general paradigm below.

  %  Vector of tissue demand levels
  % M0v = [1 4 M0c 12 20];
  M0v = [1 M0c 20];
%   M0v = M0c;
%   M0v = [1e-3,1:0.1:20];
  %M0v= [M0c 12 20];
  %M0v = [M0c  M0c 4 4 12 12];
  formatvars.M0v = M0v;
  %  occlusion_level-1 total occlusion, 0 no occlusion
  occlusion_levels = [0,1];  
  
  [metadata.ocM,metadata.M0M] = ndgrid(occlusion_levels,M0v);
  nocc = numel(occlusion_levels);
  nM0 = numel(M0v);
%   sols(nocc,nM0).misc.occlev = occlusion_levels(end);
  %  Runs all the simulations
  %  Note for faster runtimes, switch "for" to "parfor"
  for simc = 1:numel(metadata.ocM)

    %  Note: finalsystemstate is final system state at end time (one time)
    %  across all x's (formatvars.nodeslength) while sols contains solution
    %  information across all times in the run at a select few points in
    %  the spatial grid (midpoint, maybe midpoint and a couple of
    %  endpoints, etc.)
    %  This is the function where the ODE solvers are used to find D and A
    [timedata{simc},spacedata{simc}] =...
      DAsolver(Pc,Dc,Ac,nc,Lc,muc,Prmc,Rc,metadata.M0M(simc),metadata.ocM(simc),...
      Length_constant,multiplier,capdensity,respswitches,krogh,...
      formatvars,condrespsignalc,tauc);
      
  end
  
  if formatvars.save_curr_set_of_runs
    save curr_run timedata spacedata metadata
  end
  
  %  Plots info from those simulations
  postprocessor(timedata,spacedata,metadata)
    
  return
  
  %  Some handy post-processing
  
      %         keyboard
      %
      %         load DAvals.mat
      %         load DAallvals2.mat
      %  ??? Perfusion in calf I believe
      % perfusion = QFEED*n.nFEED/krogh.totalvolume*6000;
      flow_tissue(iamhere,simc) = QSS.QTotalCalf/krogh.totalvolume*6000;
      consja(iamhere,simc) = constotncu;

      %for autoreg graphs:
      CRsignalMID(iamhere,simc) = goodmeta.MID;
      CRsignalFEED(iamhere,simc) = goodmeta.FEED;
      CRsignalCOL(iamhere,simc) = goodmeta.COL;

      StonevalueMID(iamhere,simc) = SSDTtPAQ.MID(iamhere,2)*Prmc.MID.Cmyo-...
        SSDTtPAQ.MID(iamhere,3)*Prmc.MID.Cshear-goodmeta.MID+...
        Prmc.MID.Cpptone;
      StonevalueFEED(iamhere,simc) = SSDTtPAQ.FEED(iamhere,2)*Prmc.FEED.Cmyo-...
        SSDTtPAQ.FEED(iamhere,3)*Prmc.FEED.Cshear-goodmeta.FEED+...
        Prmc.FEED.Cpptone;
      StonevalueCOL(iamhere,simc) = SSDTtPAQ.COL(iamhere,2)*Prmc.COL.Cmyo-...
        SSDTtPAQ.COL(iamhere,3)*Prmc.COL.Cshear-goodmeta.COL+...
        Prmc.COL.Cpptone;

      %satstarta(ii) = saturationtotal(1);
      satstartCOL(simc) = finalsystemstate(simc).endpts.O2Hbsat(1);
      satstartFEED(simc) = finalsystemstate(simc).endpts.O2Hbsat(2);
      satstartMID(simc) = finalsystemstate(simc).endpts.O2Hbsat(3);
      satstartcap(simc) = finalsystemstate(simc).endpts.O2Hbsat(4);
      satstartvMID(simc) = finalsystemstate(simc).endpts.O2Hbsat(5);

      %  Tracks tension, shear stress, metabolic signal, Stone, diameter, activation, flow, and pressures with consumption

      TensionMID(simc) = SSDTtPAQ.MID(2);
      TensionFEED(simc) = SSDTtPAQ.FEED(2);
      TensionCOL(simc) = SSDTtPAQ.COL(2);
      myoMID(simc) = TensionMID(simc)*Prmc.MID.Cmyo;
      myoFEED(simc) = TensionFEED(simc)*Prmc.FEED.Cmyo;
      myoCOL(simc) = TensionCOL(simc)*Prmc.COL.Cmyo;
      shearMID(simc) = SSDTtPAQ.MID(3);
      shearFEED(simc) = SSDTtPAQ.FEED(3);
      shearCOL(simc) = SSDTtPAQ.COL(3);
      tautermMID(simc) = shearMID(simc)*Prmc.MID.Cshear;
      tautermFEED(simc) = shearFEED(simc)*Prmc.FEED.Cshear;
      tautermCOL(simc) = shearCOL(simc)*Prmc.COL.Cshear;
      CondMID(simc) = goodmeta.MID;
      CondFEED(simc) = goodmeta.FEED;
      CondCOL(simc) = goodmeta.COL;
      SignalCOL(simc) = int.COL;
      SignalFEED(simc) = int.FEED;
      SignalMID(simc) = int.MID;
      StoneMID(simc) = StonevalueMID(iamhere,simc);
      StoneFEED(simc) = StonevalueFEED(iamhere,simc);
      StoneCOL(simc) = StonevalueCOL(iamhere,simc);
      DiamMID(simc) = SSDTtPAQ.MID(1);
      DiamFEED(simc) = SSDTtPAQ.FEED(1);
      DiamCOL(simc) = SSDTtPAQ.COL(1);
      ActMID(simc) = SSDTtPAQ.MID(5);
      ActFEED(simc) = SSDTtPAQ.FEED(5);
      ActCOL(simc) = SSDTtPAQ.COL(5);
      FlowMID(simc) = SSDTtPAQ.MID(6);
      FlowFEED(simc) = SSDTtPAQ.FEED(6);
      FlowCOL(simc) = SSDTtPAQ.COL(6);
      Flowcalf(simc) = SSDTtPAQ.totalcalf(1);
      PressMID(simc) = SSDTtPAQ.MID(4);
      PressFEED(simc) = SSDTtPAQ.FEED(4);
      PressCOL(simc) = SSDTtPAQ.COL(4);
      Totaltime = SSDTtPAQ.totaltime;

      %STOPPED HERE MYSON 1/17/2017 5pm
      %comparing with rest
      %        changeDFEED(ii) = (DiamFEED(ii) - 0.006659)/0.006659;
      %        changeDMID(ii) = (DiamMID(ii) - 0.001151)/0.001151;
      %comparing with control state
      %        changeDFEED(ii) = (DiamFEED(ii) - 0.008525)/0.008525;
      %        changeDMID(ii) = (DiamMID(ii) - 0.001875)/0.001875;
      %comparing with M_0 = 2
      %         changeDFEED(ii) = (DiamFEED(ii) - 0.007042)/0.007042;
      %         changeDMID(ii) = (DiamMID(ii) - 0.001362)/0.001362;

      %gamA(ii) = gammaA;
      gamCOL(simc) = GSS.GCOL;
      gamFEED(simc) = GSS.GFEED;
      gamMID(simc) = GSS.GMID;
      gamCAP(simc) = GSS.GCAP;
      gamvMID(simc) = GSS.GvMID;
      gamvFEED(simc) = GSS.GvFEED;
      %gamV(ii) = gammaV;



%     end

    save my_comp times sols;

    %         Again = Again - 1;
    %         Start = Start2;
    %         End = End2;
    %         Inc = Inc2;
    % end

    % end


    %consja
    %flow_tissue

    % MYSON: Keyboard here to see what plots actually matter
    % keyboard

    tmpx = consja(:).*6000;  tmpxlabel = 'Consumption';
    tmpx = M0M;  tmpxlabel = 'Demand';

    if all(ocM)
      my_suff = '--';
      occludedy = Flowcalf;
    else
      my_suff = '-';
      nonoccludedy = Flowcalf;
    end
  
  %  Helpful things to store before trying to plot stuff
  
  my_plotter(sols);

end