function [alpha, x, f, grad, fail, nsteps] = lszoomAMPL(...
      lo, hi, flo, fhi, glo, ghi, f0, g0, x0, d, fgname, c1, c2, prtlevel)
% call: [alpha, x, f, grad, fail, nsteps] = lszoom(...
%     lo, hi, flo, fhi, glo, ghi, f0, g0, x0, d, pars, c1, c2, prtlevel)
% based on Nocedal Wright line search Zoom subroutine, see their book
% strong Wolfe line search with cubic interpolation, does accurate line
% search efficiently if 2nd Wolfe parameter c2 is small
% intended only to be called by linesch_sw.m
% input
%  lo and hi are two points bracketing the interval we want to search,
%  with corresponding function values flo <= fhi and derivatives glo, ghi
%   (these are all passed to avoid reevaluating the function unnecessarily)
%  f0 and g0 are the function and gradient values at the initial point alpha=0
%  x0 and d are the current x-vector and direction d respectively
%  pars is a parameter structure, including pars.fgname, the function name
%  c1 and c2 are Wolfe line search parameters
%  prtlevel is the print level
% output
%  alpha, x, f, grad: final steplength, x, function and gradient values
%  fail: 0 if Wolfe conditions satisfied, 1 if not
%  nsteps: number of steps taken by zoom

%  The following is an enhancement to Nocedal-Wright's routine.
%  When doing successive cubic interpolation, if we want to achieve
%  quadratic convergence, we have to keep track of 3 points:
%     lo: the best so far
%     lo2: the previous value of lo (initialized to hi)
%     hi: we guarantee a point satisfying the Wolfe conditions lies between
%         the points lo and hi
%  We must use the two best points for the interpolation.  This might not 
%  involve the point hi.  This is because the successive cubic interpolation 
%  may converge from one side, so that hi is not updated at all.  
%  But we still need hi, to make sure our cubic interpolant doesn't go 
%  across the bisection point between lo and hi.
%  If it does, we'll use the bisection point instead.
% Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
% with a subject header containing the string "nlcg".
% NLCG Version 1.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  NLCG 1.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fail = 0;
lo2 = hi;  
flo2 = fhi;
glo2 = ghi;
nsteps = 0;
while lo ~= hi & nsteps < 50  % machine precision limits make lo = hi possible
   nsteps = nsteps + 1;
   bisect = (lo + hi)/2;   % bisection 
   interp = cubic_interp(lo, lo2, flo, flo2, glo, glo2); % cubic interpolation
   % want new point to be between best point lo and the bisection point
   % it's possible that interp is infinity or even NaN: inside checks this too
   if inside(interp, lo, bisect)  
      atry = interp; % alpha_try ("try" is a reserved word)
   else
      atry = bisect; % fall back on bisection if interpolation was no good
   end;
   xtry = x0 + atry*d;          % try this x next
   % [ftry,gradtry] = feval(fgname, xtry);  % get function and gradient
   [ftry, c] = spamfunc(xtry, 0);
   [gradtry, A] = spamfunc(xtry, 1);
   gtry = gradtry'*d;          % directional derivative
   if ftry > f0 + c1*atry*g0 | ftry >= flo    % gone too far
      hi = atry; 
      lo2 = hi; flo2 = ftry; glo2 = gtry;  % since hi is updated, make this lo2
   else
      if abs(gtry) <= -c2*g0     % strong Wolfe conditions are satisfied
         alpha = atry; x = xtry; f = ftry; grad = gradtry;
         return
      end
      if gtry*(hi - lo) >= 0     % check which way to look next
         hi = lo; 
      end
      lo2 = lo; flo2 = flo; glo2 = glo;   % second best is previous value of lo
      lo = atry; flo = ftry; glo = gtry;  % best is new value of lo
   end % else part
end % loop
if prtlevel > 1
    fprintf('linesch_sw: failed to satisfy strong Wolfe conditions: loop in lszoom ran 50 times\n')
end
alpha = atry; x = xtry; f = ftry; grad = gradtry; fail = 1; 
% grad at lo wasn't saved, but atry, lo and hi should all be (nearly) same if this happens
% in the case that the line search fails because it could not find a lower
% point, it seems possible that ftry is much bigger
% than f0, yet we still return this value.  This seems like a bad idea,
% since f could be crazily varying due to ill conditioning or rounding.
% However, it has not been a problem so we didn't change this.  In
% linesch_ww, we return the lowest function value, which might be f0.