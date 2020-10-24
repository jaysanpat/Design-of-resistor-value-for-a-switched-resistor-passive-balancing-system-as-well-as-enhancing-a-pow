
% First, make sure that the ESC toolbox functions are in the path
addpath readonly
load ./readonly/E2model; % The ESC model of the cell used for this project

load ./readonly/packData
hist(100*packData.storez); % Plot histogram of cell SOC in this pack (in percent)
title('Histogram of battery-pack individual SOC values before balancing');
xlabel('State-of-charge (%)'); ylabel('Count (in 100-cell battery pack)');
xlim([85 95])

% function [storez_output,pk,maxpk,tbal] = cellBalance(packData,model,Rbal)
%
% - packData - Contains parameter data for all cells before balancing
% - model    - The ESC model for the cell used in the battery pack
% - Rbal     - The balancing resistance you have chosen for this design
%
% - storez_output - The final SOC values for each cell after balancing
% - pk            - The total power dissipated by all balacing resistors
%                   in the pack as a function of time (in W)
% - maxpk         - The maximum power dissipated by any individual balancing 
%                   resistor in the pack as a function of time (in W)
% - tbal          - The time taken to balance the pack to within 2mV (in s)
% - storeV        - All cell voltages versus time... note that cells have
%                   different self-discharge rates, different leakage current,
%                   etc, so will decrease in voltage even when not being 
%                   balanced actively.
function [storez_output,pk,maxpk,tbal,saveV] = cellBalance(packData,model,Rbal)
  % ------------------------------------------------------------------------
  % Initialize some simulation configuration parameters ...
  % ------------------------------------------------------------------------
  maxtime = 36*3600; % Maximum simulation run time in simulated seconds
  pk = zeros(maxtime,1);
  maxpk = zeros(maxtime,1);
  saveV = zeros(length(packData.storez),maxtime);
  
  % ------------------------------------------------------------------------
  % Initialize states for ESC cell model. These values came from instructor-
  % executed code that simulated more than 200h of drive cycles and charging
  % (without balancing)... these are the final states of those cells after
  % that simulation, and so are the initial states of these cells before
  % balancing.
  % ------------------------------------------------------------------------
  z = packData.storez;
  irc = packData.storeirc;
  
  % ------------------------------------------------------------------------
  % Default initialization for cells within the pack
  % ------------------------------------------------------------------------
  T = packData.T;       % Temperature for each cell, assumed to be constant
  Tsd = packData.Tsd;   % Self-discharge "temperature" for each cell
  leak = packData.leak; % Leakage current for each cell
  q = packData.q;       % Total capacity for each cell (all may be different)
  rc = packData.rc;     % R-C time constant for each cell
  r = packData.r;       % Diffusion-resistor values for each cell
  r0 = packData.r0;     % Series resistance values for each cell
  
  % ------------------------------------------------------------------------
  for k = 1: maxtime
    % Calculate cell voltages
    v = OCVfromSOCtemp(z,T,model); % get OCV for each cell
    v = v - r.*irc; % add in capacitor voltages (ignore hysteresis to simplify sim)
    saveV(:,k) = v(:); % save all cells' terminal voltage for later analysis
    % Cell Simulation
    ik = zeros(size(v)); % no string current
    % Simulate self discharge via variable resistor in parallel
    rsd = ((-20+0.4*Tsd).*z + (35-0.5*Tsd))*1e3; ik = ik + v./rsd;
    % Simulate leakage current
    ik = ik + leak;
    % Check to see which cells are 2mV or more above minimum cell voltage
    checkBalance = (v - min(v)) - 2e-3 >= 0;
    if sum(checkBalance) == 0, % balancing is complete, so return
      saveV = saveV(:,1:k-1);
      pk = pk(1:k-1);
      maxpk = maxpk(1:k-1);
      tbal = k-1;
      storez_output = z; % Output only final SOC for each cell after 4 hours of simulation time
      return
    end
    % cells 2mV or more from the minimum will have array value 1, otherwise 0
    % Add balancing resistors and calculate resulting cell current
    v_balance = v.*checkBalance;
    % Set non-balance cell voltage to 0 for calculation (to ensure balance current = 0 for no-balance cells)
    i_balance = (v_balance./Rbal); % Current calculated for balance cell, with parallel resistor 
    ik = ik + i_balance; % Add balance current to externally applied cell current
    
    % Compute power
    pk(k) = sum(i_balance.^2*Rbal);    % total power dissipated by all cells being balanced in pack 
    maxpk(k) = max(i_balance.^2*Rbal); % maximum single-cell power dissipated
    % Calculate new SOC for each cell
    z = z - (1/3600)*ik./q; % Update each cell SOC
    % Update diffusion-resistor currents
    irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
    if mod(k,3600) == 0,
      fprintf('Completed %dh balancing (%d "unbalanced" cells remain).\n',round(k/3600),sum(checkBalance));
    end
  end % for k
  % If we get to this point, then must have balanced for more than 200h (not good)
  tbal = k;
  storez_output = z; % Output only final SOC for each cell 
end


% function rBalance = tuneBalancer
%
% rBalance - the value you choose for your design

function rBalance = tuneBalancer

  rBalance = 170; % [ohms] ... This is a sample value. You will need to change it.
end  

rBalance = tuneBalancer;
[storez_output,pk,maxpk,tbal,saveV] = cellBalance(packData,model,rBalance);

t = 1:length(pk);
subplot(2,2,1);
plot(t/3600,pk); xlabel('Time (h)'); ylabel('Total pack balancing power'); 
title('Total pack power versus time');

subplot(2,2,2)
plot(t/3600,maxpk); xlabel('Time (h)'); ylabel('Max. cell balancing power'); 
title('Max. cell balancing power versus time');

subplot(2,2,3)
hist(100*storez_output); % Plot histogram of post-balancing cell SOC in this pack (in percent)
title('Post-balancing histogram of SOC values'); xlabel('State-of-charge (%)');
ylabel('Count (in 100-cell battery pack)'); xlim([85 95]);

subplot(2,2,4)
plot(t/3600,saveV'); % Plot trace of all cell voltages over time
title('Cell voltages during balancing'); xlabel('Time (h)'); ylabel('Voltage (V)');

fprintf('During balancing, maximum pack balancing power was %fW (should be less than 10W)\n',max(pk));
fprintf('During balancing, maximum cell balancing power was %fW (should be less than 0.1W)\n',max(maxpk));
fprintf('Time to balance was %fh (should be less than 29.5h)\n',tbal/3600);

% Compute expected grade
grade = [];
if (max(pk) <= 10) && max(maxpk) <= 0.1 && tbal/3600 < 32.5,
  grade = max(0,min(10,find((29.5:-0.5:22)<tbal/3600,1,'first') - 1));
end
if isempty(grade), grade = 0; end
fprintf(['If you submit the assignment with this balancing resistance,\n' ...
          '  you should expect to receive a grade of %d/10\n'],grade);


%% PART 2

% First, make sure that the ESC toolbox functions are in the path
addpath readonly
load ./readonly/CellModel.mat % load ESC model


% function [pChg,pDis] = HPPCpower(z0,T,dT,eta,ns,np,model,limits)
%
% z0  -  the SOC to use when computing discharge and charge resistances RDis and RChg, 
%       and when computing available power based on SOC limits; the same value is 
%       used for every cell in the battery pack
% T   - the temperature to use when computing discharge and charge resistances RDis and RChg, 
%       and when computing available power; the same value is used for every cell in the battery pack
% dT  - the pulse duration to use when computing discharge and charge resistances RDis and RChg;
%       note that the pulse magnitude should use a 10C rate (just like example code from lesson 5.3.3)
% eta - the coulombic efficiency to use when computing available power based on SOC limits
% ns  - the number of cells in series in the battery pack
% np  - the number of cells in parallel in the battery pack
% model - an ESC model type
% limits - design limits on SOC, voltage, current, and power (see code for how these are stored)
%
% pChg - your computed value for charge power (W)
% pDis - your computed value for discharge power (W)
function [pChg,pDis] = HPPCpower(z0,T,dT,eta,ns,np,model,limits)


  zmin = limits.zMin; zmax = limits.zMax; % Retrieve SOC limits [unitless]
  vmin = limits.vMin; vmax = limits.vMax; % Retrieve voltage limits [V]
  imin = limits.iMin; imax = limits.iMax; % Retrieve current limits [A]
  pmin = limits.pMin; pmax = limits.pMax; % Retrieve design power limits [W]
  

  % NOTE: Resistance calculated using a 10C dis/charge pulse for dT samples
  % You will need to modify this code to work correctly with the specific input
  % parameter list to this function: z0, T, dT ...
  % Note that rChg and rDis are calculated correctly for the default input parameters
  % but are not calculated correctly for non-default inputs. You will need to change
  % this code to calculate rChg and rDis correctly for arbitrary function inputs
  
  Q = getParamESC('QParam',T,model); 
  iChgPulse = 10*Q*[zeros(dT/2,1); -ones(dT,1); zeros(dT/2,1)];  % [A] charge pulse
  iDisPulse = 10*Q*[zeros(dT/2,1);  ones(dT,1); zeros(dT/2,1)];  % [A] discharge pulse
  [vk,~,~,~,~] = simCell(iChgPulse,T,model,1,z0,0,0);
  rChg  = abs((max(vk)-vk(1))/min(iChgPulse));
  [vk,~,~,~,~] = simCell(iDisPulse,T,model,1,z0,0,0);
  rDis  = abs((min(vk)-vk(1))/max(iDisPulse));


  % equations need to be implemented.)
 % pDis = 1;  % You will need to change this to compute it correctly
 % pChg = -1; % You will need to change this to compute it correctly

soc = z0 ;

% HPPC Power Estimation: Truth
OCV      = OCVfromSOCtemp(soc,T,model);
iDisMaxV = (OCV-vmin)/rDis;
iDisMaxZ = (soc - zmin)*3600*Q/dT;
iDisMax  = max(0,min([iDisMaxV;iDisMaxZ;imax*ones(size(soc))]));
pDisMax  = min(vmin*iDisMax,pmax*ones(size(soc)));
iChgMinV = (OCV-vmax)/rChg;
iChgMinZ = (soc - zmax)*3600*Q/eta/dT;
iChgMin  = max([iChgMinV;iChgMinZ;imin*ones(size(soc))]);
pChgMin  = min(0,max(vmax*iChgMin,pmin*ones(size(soc))));
pDis = pDisMax;
pChg = pChgMin;


end


default.z0 = 0.5;
default.T = 25;
default.dT = 10;
default.eta = 1;
default.ns = 1;
default.np = 1;
limits.zMin = 0.1;
limits.zMax = 0.9;
limits.vMin = 2.8;
limits.vMax = 4.3;
limits.iMin = -200;
limits.iMax = 350;
limits.pMin = -1000;
limits.pMax =  1000;
default.limits = limits;
[pChg,pDis] = HPPCpower(default.z0,default.T,default.dT,default.eta,default.ns,default.np,model,default.limits)
% Note that the correct answer for the default set of limits is:
% pChg = -385.18
% pDis = 885.55
% It may also be helpful to know that the correct values for rChg and rDis for the default set of parameters is
% rChg = 3.6787 mOhm
% rDis = 3.7009 mOhm
