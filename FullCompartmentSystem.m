function [timedata,spacedata,metadata] = FullCompartmentSystem(...
    myoLAin,myoSAin, shearLAin,shearSAin,...
    Cmeta,metaLAin, metaSAin, GCO2_LA, GCO2_SA,Length_constant,multiplier,capdensity)

%%% This is the main function that runs the whole simulation

% Inputs
%  myoLAin-on or off:  myogenic response in LA vessels
%  myoSAin-on or off:  myogenic response in SA vessels
%  shearLAin-on or off:  shear response in LA vessels
%  shearSAin-on or off:  shear response in SA vessels
%  Cmeta-default {200}:  magnitude of the metabolic response
%  Length_constant-default {1cm}:  length constant associated with
%    exponential decay of the conducted response upstream
%  multiplier-default {1}:  as conducted response signal moves upstream it
%    gets multiplied from 1 node into multiplier nodes along the way
%    somehow?  Not yet used
%  capdensity-default {500}:

% close all
%  Default values for inputs (if not fed in, these will be used,
%  particularly helpful for debugging purposes)

%%% This section is to set up what responses we want to run. They are the
%%% boolean variables
if nargin == 0
    myoLAin = true; myoSAin = true;
    shearLAin = true; shearSAin = true;
    Cmeta = 200; metaLAin = true; metaSAin = true;
    GCO2_LA = true; GCO2_SA = true;
    M0vary = false; Pavary = true;
    Length_constant = 1; multiplier = 1; capdensity = 500;
    %%% Drop the capillary density by x%
    % reduction = 0;
    % capdensity = (1-reduction)*capdensity;
    % needtosavecontrol = true; %make true when control state needs to be updated/run, otherwise can be false

end
if multiplier > 1
    error('multiplier > 1 does not seem to be implemented yet');
end
%  Throw myogenic and shear response switches into a nice structure that we
%  can pass back and forth.  Note that the control state uses all
%  responses while non-control states either allow the response levels
%  to change (true) or always use the control state response levels
%  (false; response levels constant in this case).

%%% Store these boolean variable responses into a structure so we can use
%%% them later
respswitches.myo.LA = myoLAin;
respswitches.myo.SA = myoSAin;
respswitches.shear.LA = shearLAin;
respswitches.shear.SA = shearSAin;
respswitches.CO2.LA = GCO2_LA;
respswitches.CO2.SA = GCO2_SA;
%  Just in case we want to turn off the metabolic response later on.
respswitches.meta.LA = metaLAin; %metaFEEDin;
respswitches.meta.SA = metaSAin; %metaMIDin;
% Whether we want to test for elevated IOP or not

%  Default values for nice pictures...shouldn't affect calculations
formatvars = initializeformatvars;


%defining vector of IOP values for varying IOP
IOPv = [15*1333 25*1333];
% IOPv = [15*1333 25*1333];
nIOP = numel(IOPv);


%  Define system control state values using experimental-based info.
%%% Get the control state values and put them in their respective structures
[tauc,Pc,muc,M0c,Prmc] = SetCState(capdensity);
%  Use previous values along with some constraints (e.g. activation
%  levels are 0.5 in the LA and SA arterioles in the control state...15
%  changeable ContAct values inside of
%  SetCState) to obtain other control state values

%%% This uses the control state values from SetCState, along with other
%%% equations to calculate more Cstate values
%%% Lots of stuff using the equations defined in publications 2 and 3
[Dc,tauc,Pc,Ac,nc,Lc,muc,M0c,Prmc,Rc,Cc,Qc,...
    condrespsignalc,krogh, rref, etac] = GetCState(Cmeta,Length_constant,multiplier,...
    Pc,tauc,muc,M0c,Prmc,formatvars,respswitches,Pc.IOP);
%  Below we loop through a different set of occlusions (partial occlusion
%  can now occur).  For each "occlusion level" we loop through multiple
%  M0 values.  Might or might not want to change this later so that one
%  can have different values for M0 depending on the occlusion level but
%  right now that doesn't happen.  Also note that in the past we also
%  looped through different Pa values.  Additional parameters can be
%  added using the general paradigm below.
%
%
% %  Vector of tissue demand levels
if M0vary
    % M0v = [0.5 M0c 2 3 4 5 6];
     % M0v = [M0c 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3];
    % M0v = [0.3 1 3];
    M0v = [.3 .6 M0c  1.3  1.7  2  2.3  2.5 2.7  3];
    % M0v = [1];
    formatvars.M0v = M0v;

    nM0 = numel(M0v);

    %a string of color codes for M0 plot curves
    colors = 'bgcmgcbgr';

    for j = 1:nIOP


        for simc = 1:nM0  %make this PARFOR!!!!!!!!!!!!!!

            [timedata{simc},spacedata{simc}] =...
                DAsolver(Pc,Dc,Ac,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0v(simc),...
                Length_constant,multiplier,capdensity,respswitches,krogh,...
                formatvars,condrespsignalc,tauc,1,rref,IOPv(j),Pc);

            %                 [timedata{simc},spacedata{simc}] =...
            %                 DAsolver(Pvar,Dnew,Anew,nc,Lc,muc,Prmcnew,Rnew,Cnew,Qnew,M0v(simc),...
            %                 Length_constant,multiplier,capdensity,respswitches,kroghnew,...
            %                 formatvars,condrespsignalnew,taunew, 1, rref);
            %
            %             if needtosavecontrol && j==1
            %                 [timedatac{simc},spacedatac{simc}] =...
            %                 DAsolver(Pc,Dc,Ac,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0v(simc),...
            %                 Length_constant,multiplier,capdensity,respswitches,krogh,...
            %                 formatvars,condrespsignalc,tauc,0, rref);
            %
            %                 if IOPv(j)==15*1333
            %                     colors = 'b';
            %                 end
            %
            %             end
            %         end
            %
            %         if needtosavecontrol
            %             save('controldataM0', 'timedatac', 'spacedatac');
            %         else
            %             load controldataM0.mat;
            %         end
            %
            % postprocessorM0elevated(timedatac,spacedatac, timedata, spacedata, M0v, muc, colors(j))


            %
            %
            %
        end

        postprocessorM0(timedata, spacedata, M0v, muc, colors(j), Pc.Pa, IOPv(j))
        %
    end
end

%
%
%

if Pavary
    %%% Vector of many different Pa values to work with
        % Pav = [25*1333 26*1333 27*1333 28*1333 29*1333 30*1333 31*1333 32*1333 33*1333 34*1333 35*1333 36*1333 37*1333 38*1333 39*1333 40*1333 41*1333 42*1333 43*1333 44*1333 45*1333 46*1333 47*1333 48*1333 49*1333 50*1333 51*1333 52*1333 53*1333 54*1333 55*1333];
        % Pav = [35*1333 36*1333 37*1333 38*1333 39*1333 40*1333 41*1333 42*1333 43*1333 44*1333 45*1333 46*1333 47*1333 48*1333 49*1333 50*1333 51*1333 52*1333 53*1333 54*1333 55*1333];
        % Pav = [32*1333 33*1333 34*1333 35*1333 36*1333 37*1333 38*1333 39*1333 40*1333 41*1333 42*1333 43*1333 44*1333 45*1333 46*1333 47*1333 48*1333 49*1333 50*1333 51*1333 52*1333 53*1333 54*1333 55*1333];


     Pav = [35*1333 38*1333 40*1333  42*1333 45*1333 48*1333 50*1333 55*1333];

     % Pav = [38*1333 40*1333  42*1333 45*1333 48*1333];
     % Pav = [30*1333 31*1333 32*1333 33*1333 34*1333 35*1333  38*1333 40*1333  42*1333 45*1333]; 
     % Pav = [50*1333 51*1333 52*1333 53*1333 54*1333];

     % Pav = [35*1333 40*1333 50*1333];

    % % Pav = [ 20*1333 23*1333 26*1333 29*1333 32*1333 35*1333 ];
        % Pav = [25*1333 27*1333  30*1333 32*1333  35*1333  38*1333 40*1333  42*1333 45*1333 48*1333 50*1333 55*1333];
    % Pav = [33*1333 34*1333 35*1333 36*1333];

    formatvars.Pav = Pav;

    nPa = numel(Pav);
    %%% Repeat this for each element of the Pa values vector
    for j = 1:nIOP

        % [Pvar, taunew, Dnew, Tnew, Cnew, Rnew, etanew, Qnew, kroghnew, condrespsignalnew, Anew, Prmcnew] = Get_Elev_IOP_State(Cmeta,Length_constant,multiplier,Pc,Qc,nc,Lc,Dc,tauc,muc,M0c,...
        % Prmc, formatvars, respswitches, IOPv(j), rref);

        %string of colors for plots with varying Pa;
        colors = 'bgmrg';

        for simc = 1:nPa
            %%% Set the control state Pa value to the current one in the
            %%% iteration
            % Pc.Pa = Pav(simc);
            Pvar.Pa = Pav(simc);

            %  Note: finalsystemstate is final system state at end time (one time)
            %  across all x's (formatvars.nodeslength) while sols contains solution
            %  information across all times in the run at a select few points in
            %  the spatial grid (midpoint, maybe midpoint and a couple of
            %  endpoints, etc.)
            %  This is the function where the ODE solvers are used to find D and A

            % [timedata{simc},spacedata{simc}] =...
            %     DAsolver(Pnew,Dnew,Anew,nc,Lc,muc,Prmcnew,Rnew,Cnew,Qnew,M0c,...
            %     Length_constant,multiplier,capdensity,respswitches,kroghnew,...
            %     formatvars,condrespsignalnew,taunew,1, rref);


            [timedata{simc},spacedata{simc}] =...
                DAsolver(Pc,Dc,Ac,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0c,...
                Length_constant,multiplier,capdensity,respswitches,krogh,...
                formatvars,condrespsignalc,tauc,1,rref,IOPv(j),Pvar);




            % if needtosavecontrol && j==1
            %     [timedatac{simc},spacedatac{simc}] =...
            %         DAsolver(Pc,Dc,Ac,nc,Lc,muc,Prmc,Rc,Cc,Qc,M0c,...
            %         Length_constant,multiplier,capdensity,respswitches,krogh,...
            %         formatvars,condrespsignalc,tauc,0, rref,IOPv(j),Pvar);
            %
            %
            %
            %     if IOPv(j)==15*1333
            %         colors = 'b';
            %     end
            % end



        end



        % if needtosavecontrol
        %     save('controldataPa', 'timedatac', 'spacedatac');
        % else
        %     load controldataPa.mat;
        % end

        % postprocessorPaelevated(timedata,spacedata, timedatac, spacedatac, Pav,Prmc.LA.Cpptone,Prmc.SA.Cpptone, muc, colors(j))
        postprocessorPa(timedata,spacedata, Pav,Prmc.LA.Cpptone,Prmc.SA.Cpptone, muc, colors(j), Prmc, IOPv(j))


    end


    %%% Plot perfusion
    %     figure(24);
    % %  flow_tissue(iamhere,simc) = QSS.QTotalCalf/krogh.totalvolume*6000;
    % hold on;
    %    %%% Calculate perfusion for each Pa value
    %    for j = 1:length(Pav)
    % %        perfusion(j) = timedata{1, j}.meta.perfusion.LA(end);
    %        perfusion(j) = timedata{1, j}.flow.Q.QLA(end)*nc.nLA/krogh.totalvolume*6000;
    %    end
    %    plot(Pav./1333, perfusion, 'o-');
    %    xlabel('Pa')
    %    ylabel('Perfusion')
    %    title('Perfusion vs Pa')

    %%% Save the space data so we can use it elsewhere
    %%% This only works when running with only Pa = 40
    %writestruct(spacedata{1,1}, "1layerdata.xml")


    %
    % %%% Plot a bunch of graphs




end


return
