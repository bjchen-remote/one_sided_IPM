function tableau = rk6Tableau()
%IPM.EVOLVE.RK6TABLEAU Dormand--Prince eight-stage RK6(5)8M tableau.
%   The sixth-order propagation row is used by the solver; the embedded
%   fifth-order row is retained for independent defect tests and future
%   step control.  Its stability polynomial is the degree-seven Taylor
%   polynomial of exp(z), giving nonzero imaginary-axis stability.  This is
%   preferable for advection to the seven-stage Butcher RK6 process, whose
%   stability region touches the imaginary axis only at the origin.
%
%   Source: P. J. Prince and J. R. Dormand, "High order embedded Runge-Kutta
%   formulae", J. Comput. Appl. Math. 7 (1981), 67--75.

persistent storedTableau
if isempty(storedTableau)
    A = zeros(8);
    A(2,1) = 1/10;
    A(3,1:2) = [-2/81,20/81];
    A(4,1:3) = [615/1372,-270/343,1053/1372];
    A(5,1:4) = [3243/5500,-54/55,50949/71500,4998/17875];
    A(6,1:5) = [-26492/37125,72/55,2808/23375, ...
        -24206/37125,338/459];
    A(7,1:6) = [5561/2376,-35/11,-24117/31603, ...
        899983/200772,-5225/1836,3925/4056];
    A(8,1:7) = [465467/266112,-2945/1232,-5610201/14158144, ...
        10513573/3212352,-424325/205632,376225/454272,0];
    b = [61/864,0,98415/321776,16807/146016,1375/7344, ...
        1375/5408,-37/1120,1/10];
    bEmbedded = [821/10800,0,19683/71825,175273/912600, ...
        395/3672,785/2704,3/50,0];
    storedTableau = struct('A',A,'b',b,'c',sum(A,2), ...
        'bEmbedded',bEmbedded,'embeddedOrder',5, ...
        'order',6,'stages',8,'sspCoefficient',0, ...
        'source','DormandPrince1981_RK6_5_8M');
end
tableau = storedTableau;
end
