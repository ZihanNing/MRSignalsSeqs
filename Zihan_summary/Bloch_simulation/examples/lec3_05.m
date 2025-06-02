%	Short-TR IR Signal Example
%
%	(Use abprop to do this more compactly.)
%	Note we start applying A,B from TE onward as this is steady-state.

% Simple Inversion-Recovery sequence where TR is not long enough for
% full recovery so steady-state must be calculated.

%% Original code
%-- Times are all in seconds here.
TR=1;
TI = 0.5;
TE = 0.05;
T1 = 0.5;
T2 = 0.1;

%-- Shorthand - often describe E1 and E2 like this:

E1a = exp(-(TI/T1));
E2a = exp(-(TI/T2));
E1b = exp(-(TE/T1));
E2b = exp(-(TE/T2));
E1c = exp(-((TR-TI-TE)/T1));
E2c = exp(-((TR-TI-TE)/T2));

%-- Define A,B as in lecture example

A1 = diag([E2a E2a E1a]) * xrot(180);
B1 = [0;0;1-E1a];               
A2 = diag([E2b E2b E1b]) * xrot(90);
B2 = [0;0;1-E1b];
A3 = diag([E2c E2c E1c]);
B3 = [0;0;1-E1c];

%-- Could use abprop.m here

A = A2*A1*A3;                        
B = B2+A2*(B1+A1*B3);

Mss = inv(eye(3)-A)*B  

%% ZN: expressed with functions 
%-- Times are all in seconds here.
TR=1;
TI = 0.5;
TE = 0.05;
T1 = 0.5;
T2 = 0.1;

% period 1 (3 in the img; TR-TE-TI)
[A1,B1] = relax(TR-TE-TI,T1,T2,0);
A2 = xrot(180);
% period 2 (TI)
[A3,B3] = relax(TI,T1,T2,0);
% period 3 (TE)
A4 = xrot(90);
[A5,B5] = relax(TE,T1,T2,0);

% final
[A,B,Mss] = abprop([A1,B1],A2,[A3 B3],A4,[A5,B5]);
Mss