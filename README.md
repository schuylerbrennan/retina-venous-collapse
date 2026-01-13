# retina-venous-collapse
The code in this repository aims to simulate blood flow and oxygenation in a compartmental model of the retinal vasculature with venous collapse. This code was used to produce the results found in:

>Schuyler Brennan, Tajkera Khatun, Brendan Fry, Charlotte Weiss, Brent Siesky, Alice Verticchio, Alon Harris, & Julia Arciero (2026). _Modeling the impact of venous collapsibility on retinal oxygenation_. Manuscript under review. 

A more detailed description of the methods, including parameter values and sources, can be found there. 

# 🏁 How to run 

- To run the compartmental model code, call FullCompartmentSystem. The main components to be changed in this file are M0_vary and Pa_vary, which can be set to true or false. This will choose whether the simulation varies oxygen demand or incoming arterial pressure.

- M0v and Pav are two vectors that can be changed to contain the values of M0 (oxygen demand) and Pa (incoming arterial pressure) that you would like to include in each simulation.

- The vector IOPv should contain a list of the IOP values you would like to run. The IOP values should be in dyn/cm^2.

- The variable colors can be changed to contain a string of color codes indicating the color of the plot curves. When plotting, the order of the color codes will follow the order of the IOP values in IOPv.

# :bulb: Turning venous collapsibility on/off :bulb:

- To turn venous collapsibility on and off, you will need to go to obtainnewflowquantselev. To turn collapsibility on, uncomment the if statements in lines 147-170 and comment out lines 174-180. To turn collapsibility off, comment out lines 147-170 and uncomment lines 174-180.
- These lines need to be uncommented to turn venous collapsibility on:

```
  if Pvar.TMLV >= 0
       Rvar.RLV = (((Prmc.kr.LV*Lc.Lefflv)/((pi*rref.LV^2))^2)*(1+(Pvar.TMLV/(Prmc.kp.LV*Prmc.kl.LV)))^-4)/nvar.nLV;
       Rvar.alphaLV = ((Pvar.TMLV)/(Prmc.kp.LV*Prmc.kl.LV) + 1)^2;
       Rvar.ALV = Rvar.alphaLV*pi*rref.LV^2;
  else
       Rvar.RLV = (((Prmc.kr.LV*Lc.Lefflv)/((pi*rref.LV^2))^2)*(1-(Pvar.TMLV/(Prmc.kp.LV)))^(4/3))/nvar.nLV;
       Rvar.alphaLV = (1 - (Pvar.TMLV)/(Prmc.kp.LV))^(-2/3);
       Rvar.ALV = Rvar.alphaLV*pi*rref.LV^2;
  end
```

	
- These lines need to be uncommented to turn venous collapsibility off:


```
        Rvar.RSV = (((Prmc.kr.SV*Lc.Leffsv)/((pi*rref.SV^2))^2)*(1+(Pvar.TMSV/(Prmc.kp.SV*Prmc.kl.SV)))^-4)/nvar.nSV;
        Rvar.alphaSV = ((Pvar.TMSV)/(Prmc.kp.SV*Prmc.kl.SV) + 1)^2;
        Rvar.ASV = Rvar.alphaSV*pi*rref.SV^2;

        Rvar.RLV = (((Prmc.kr.LV*Lc.Lefflv)/((pi*rref.LV^2))^2)*(1+(Pvar.TMLV/(Prmc.kp.LV*Prmc.kl.LV)))^-4)/nvar.nLV;
        Rvar.alphaLV = ((Pvar.TMLV)/(Prmc.kp.LV*Prmc.kl.LV) + 1)^2;
        Rvar.ALV = Rvar.alphaLV*pi*rref.LV^2;
```

- In postprocessorPa and postprocessorM0, there are for loops that are commented out that will prevent the code from plotting trials with negative transmural pressure in the capillaries past a certain threshold. To use this feature, uncomment the for loop, change the input to the function to read Pa1 = varargin{3} or M01 = varargin{3}, and change to n = numel(Pa1) or n = numel(M01).


