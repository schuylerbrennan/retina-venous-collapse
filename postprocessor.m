function postprocessor(varargin)

  %  I envision this program will take in either
  %    a string-name of the file where time and space data are stored
  %    three structures-timedata, spacedata, and metadata
  %  In addition, I envision that other program options (like which plots to
  %  include and which ones not to include) might also be handed in through
  %  the "varargin" capability

  %  First argument is a string and the name of the mat file where time and
  %  space data are stored
  if ischar(varargin{1})
   %  Should load in three structures, timedata, spacedata,
   %  and metadata (see prototypeSystem.m)
   load(varargin{1});
   varargin = {varargin{2:end}};
  else
   timedata = varargin{1};
   spacedata = varargin{2};
   metadata = varargin{3};
   varargin = {varargin{4:end}};
  end

  %  Define various default program options
  options.do_plots = true;

  %  Handle other trailing varargin arguments (if there are any)
  for vac = 1:2:numel(varargin)
   options.(varargin{vac}) = varargin{vac+1};
  end

  %  Get some useful info for later
  vessel_fields_veins = fieldnames(spacedata{1}.meta(1).O2Hbsat);
  %  Collateral removed!
  vessel_fields = {'MID','FEED'};
  cols = 'rbkgmc';
  fs = 8;

  %  Do some plots...still to flush out more
  %%  Saturation levels across network at beginning and end of simulations
  close(figure(16)); figure(16);
  %  Four quantities: Oxygen saturation in the RBCs, ATP concentration in
  %  the vessels, the signal contribution (how much a given location
  %  downstream contributes to Scr in an upstream location, the "integrand"
  %  evaluated at a given downstream location for a given upstream
  %  location), and the signal accumulation (Scr after integration at a
  %  given upstream location)
  %  First do a precalculation to find siggene and sigaccu before plotting
  %  Collect x values (assume same for each simulation)
  xs = [];
  M0s = unique(metadata);
  mycols = 'rgbcmy';
  for simc = 1:numel(metadata)
    ATPconc = cell(2,1);
    if simc == 1
      fn = fieldnames(spacedata{simc}.grid.x);
      for fnc = 1:numel(fn)
        x = spacedata{simc}.grid.x.(fn{fnc});
        xs = [xs;x];
        xmid.(fn{fnc}) = mean(x);
      end
    end
    [xs,ia,ic] = unique(xs);
    distmatrix = xs-xs';
    for timepoints = 1:2
      for fnc = 1:numel(fn)
        ATPconc{timepoints} = [ATPconc{timepoints};...
          spacedata{simc}.meta(timepoints).ATPconc.(fn{fnc})];
      end
      ATPconc{timepoints} = ATPconc{timepoints}(ia);
    end
    for timepoints = 1:2
%       integrandmatrix{timepoints} = ATPconc{timepoints}.*...
%         (distmatrix >= 0).*exp(-distmatrix./...
%         metadata.inputs.Lengthconstant);
      %  Signal contribution relative to 3 locations, collateral, feed, and
      %  mid
%       for vn = {'COL','FEED','MID'}
%         spacedata{simc}.meta(timepoints).sigcont.(vn{1}) = ...
%           interp1(xs,integrandmatrix{timepoints}',xmid.(vn{1}));
%       end
%       tmpsigaccu = sum(diff(xs).*(integrandmatrix{timepoints}(1:end-1,:)...
%         +integrandmatrix{timepoints}(2:end,:))./2,1);
      for fnc = 1:numel(fn)
        spacedata{simc}.meta(timepoints).sigaccu.(fn{fnc}) = interp1(...
          xs,tmpsigaccu,spacedata{simc}.grid.x.(fn{fnc}),'linear',...
          'extrap');
      end
    end
  end
  myfields = {'O2Hbsat','ATPconc','sigcont','sigaccu'};
  for simc = 1:numel(metadata.ocM)
    mycol = mycols(find(metadata.M0M(simc) == M0s));
   %  Non-occluded is on top and occluded on bottom
   for myfc = 1:numel(myfields)
     myfield = myfields{myfc};
     subplot(4,2,1+metadata.ocM(simc)+2*(myfc-1));
     for vfvc = 1:numel(vessel_fields_veins)
       name = vessel_fields_veins{vfvc};
       if isequal(myfields{myfc},'sigcont')
         if vfvc < 4
           fiebef = spacedata{simc}.meta(1).sigcont.(name);
           fieaft = spacedata{simc}.meta(1).sigcont.(name);
           inds = fiebef > 0;
           x = xs(inds);
           fiebef = fiebef(inds);
           fieaft = fieaft(inds);
         else
           break
         end
       else
         x = spacedata{simc}.grid.x.(name);
         fiebef = spacedata{simc}.meta(1).(myfield).(name);
         fieaft = spacedata{simc}.meta(2).(myfield).(name);
       end
       line_hands = plot(x,fiebef,[mycol,':'],...
         x,fieaft,mycol,...
         x([1,numel(x)]),fiebef([1,numel(x)]),[mycol,'x'],...
         x([1,numel(x)]),fieaft([1,numel(x)]),[mycol,'x'],...
         'LineWidth',2);
       %      plot(x,sataft,cols(vfvc),x,sataft,cols(vfvc),'LineWidth',2);
       if (simc == 1) && (vfvc == 1)
         legend('Before','After','Location','southwest');
         for j = 3:4
           line_hands(j).Annotation.LegendInformation.IconDisplayStyle...
             = 'off';
         end
       else
         for j = 1:4
           line_hands(j).Annotation.LegendInformation.IconDisplayStyle...
             = 'off';
         end
       end
       hold on
     end
     text(x(1),fiebef(1),['M0 = ',num2str(metadata.M0M(simc))],...
       'FontSize',fs);
     text(x(1),fieaft(1),['M0 = ',num2str(metadata.M0M(simc))],...
       'FontSize',fs);
     if simc == find(~metadata.ocM(:),1,'last')
       title([myfield,' non-occluded']);
     elseif simc == find(metadata.ocM(:),1,'last')
       title([myfield,' occluded']);
     end
   end
  end

  %%  Total flow to the calf/thigh as function of time
  %  close(figure(17)); 
  figure(17);
  hold on
  for simc = 1:numel(metadata.ocM)
   %  Non-occluded is on top and occluded on bottom
   mycell = {'COL','CIA','Calf','EF','IIA','THIGH'};
   t = timedata{simc}.t;
   for mycellc = 1:numel(mycell)
     subplot(6,2,2*(mycellc-1)+2-metadata.ocM(simc));
     hold on
     Qtmp = timedata{simc}.flow.Q.(['QTotal',mycell{mycellc}]);
     plot(t,Qtmp,cols(mycellc));
     text(t(end),Qtmp(end),['M0 = ',num2str(metadata.M0M(simc))],...
       'FontSize',fs);
     if metadata.ocM(simc)
       title([mycell{mycellc},'-total flow occluded']);
     else
       title([mycell{mycellc},'-total flow non-occluded']);
     end
     hold on
   end
  end
 
  %%  Plot the diameters
  for vfc = 1:numel(vessel_fields)  
    figure(17+vfc);
    for simc = 1:numel(metadata.ocM)
      t = timedata{simc}.t;
      fn = vessel_fields{vfc};
      subplot(4,2,0+1+metadata.ocM(simc));
      plot(t,timedata{simc}.var.D.(fn).*10.^4,'LineWidth',2);
      hold on
      if metadata.ocM(simc)
        title('Occluded');
      else
        title('Non-Occluded');
      end
      ylabel([fn,' diameter (\mu{m})']);
      subplot(4,2,2+1+metadata.ocM(simc));
      plot(t,timedata{simc}.var.A.(fn),'LineWidth',2);
      hold on
      ylabel([fn,' activation (frac)'])
      subplot(4,2,4+1+metadata.ocM(simc));
      plot(t,timedata{simc}.var.DR.(fn),'LineWidth',2);
      hold on
      ylabel([fn,' ref diameter (\mu{m})'])
      subplot(4,2,6+1+metadata.ocM(simc));
      plot(t,timedata{simc}.var.n.(fn),'LineWidth',2);
      hold on
      ylabel(['# ',fn])
      xlabel('time (sec)')
    end
    legend(num2str(metadata.M0M(1,:)'));
  end
  keyboard
  
end