% SQP    Schittkowski's Sequential Quadratic Programming method
%        to find the constrained minimum of a function of several variables
%
% Copyright (c) 2020, Robert A. Canfield, Mark Spillman. All rights reserved.
%                     See accompanying LICENSE.txt file for conditions.
%
% Schittkowski (1985) "NLPQL: A FORTRAN Subroutine Solving Constrained
% Nonlinear Programming Problems," Annals Ops. Research, 5:485-500.
%
%          Optimization Toolbox Version 1-compatible input arguments:
%
%  usage: [x,out,v,H,status]=sqp(Fun,x0,opts,vlb,vub,Grd,P1,P2,...)
%
%  input: Fun    - string name of a function file which returns the
%                  value of the objective function and a vector of
%                  constraints (i.e. [f,g]=fun(x)).  f is minimized
%                  such that g<zeros(g). Set g=[] for unconstrained.
%
%         x0     - initial vector of design variables
%
%         opts   - (optional) vector of parameters, or optimset structure
%
%                opts(1)-Display parameter (Default: 0, none)
%                opts(2)-Termination tolerance for X. (Default: 1e-4).
%                opts(3)-Termination tolerance on F.  (Default: 1e-4).
%                opts(4)-Termination criterion on constraint violation. (Default: 1e-6)
%                opts(5)-scale design variables if <0  (or opts.Scale)
%                        scale functions if f,g>abs(opts(5))
%                opts(6)-change termination criteria
%                        (or opts.Termination)
%                opts(7)-maximum function evaluations in line search
%                        (or opts.MaxLineSearchFun)
%                opts(9)-Derivative Check (Default 0, off)
%                opts(14)-max number of function evaluations
%                opts(15)-max iterations
%                Type help foptions for more details.
%
%            Or, structure of fmincon options (optimset)
%                opts.Display         = opts(1)
%                opts.TolX            = opts(2)
%                opts.TolFun          = opts(3)
%                opts.TolCon          = opts(4)
%                opts.DerivativeCheck = opts(9)
%                opts.MaxFunEvals     = opts(14)
%                opts.MaxIter         = opts(15)
%                opts.DiffMinChange   = opts(16)
%                opts.DiffMaxChange   = opts(17)
%                opts.OutputFcn = function handle for user output/plot
%
%                In addition to optimset options, opts may contain:
%                opts.foptions    - vector (<=18 length) of old style foptions
%                opts.Scale       - opts(5) variable/function scaling
%                opts.Termination - opts(6) convergence criteria
%                opts.nec         - opts(13) number equality constraints
%                opts.LagrangeMultipliers - initial Lagrange multiplier estimate
%                opts.HessMatrix          - initial positive-definite Hessian estimate
%                opts.HessFun     - user-supplied Hessian function handle
%                     H=Hessian(x,LagrangeMultipliers)
%
%         vlb    - (optional) vector of lower bounds on the design
%                  variables
%         vub    - (optional) vector of upper bounds on the design
%                  variables
%         Grd    - (optional) string name of a function file which
%                  returns a vector of function gradients and a
%                  matrix of constraint gradients
%                  (i.e. [fp,gp]=grd(x)).
%         Pn     - (optional) variables directly passed to fun and grd
%                  optional inputs Pn can be skipped by inputing []
%
%  output: x      - vector of design variables at the optimal solution
%          out    - final program parameters, vector if opts was vector:
%                   out(8)  = value of the function at the solution
%                   out(10) = number of function evaluations
%                   out(11) = number of gradient evaluations
%                   out(15) = number of iterations
%                   fields same as fmincon output, if opts was structure
%          v      - vector of Lagrange multipliers at the solution
%          H      - Hessian at the solution
%          status - Termination status: 1=converged
%
%  fmincon-compatible problem structure input argument
%          optimtool GUI option "Export to Workspace" dialog box
%          sends problem information to the MATLAB workspace as a structure
%
%          usage: [x,opts,v,H,status]=sqp( problem )
%
%          input: problem - Data structure with fields:
%                 objective - Objective function
%                 x0        - Initial point for x
%                 Aineq     - Matrix for linear inequality constraints
%                 bineq     - Vector for linear inequality constraints
%                 Aeq       - Matrix for linear equality constraints
%                 beq       - Vector for linear equality constraints
%                 lb        - Vector of lower bounds
%                 ub        - Vector of upper bounds
%                 nonlcon   - Nonlinear constraint function
%                 options   - Options created with optimset
%
%  Written by:   Capt Mark Spillman and Maj Robert A. Canfield
%                Air Force Institute of Technology, Virginia Tech
%                Octave-compatiby functions by Blake van Winkle
%  e-mail:       bob.canfield@vt.edu
%
%  Created:      12/5/94
%  Modified:      6/4/20
%
% The function format is based on the MATLAB function constr.m written
% by Andy Grace of MathWorks, 7/90.  The algorithm is based the FORTRAN
% routine NLPQL written by Klaus Schittkowski, 6/91.
%
%---------------------------------------------------------------------
% Explain the different possible termination criteria
%
% Three different termination criterias can be selected with opts(6):
%
% 1.  If opts(6)=(-1), Schittkowski's criteria is used:
%        KTO=abs(s'*fp)+sum(abs(u.*gv))  <=  opts(3)
%                       SCV=sum(g(g>0))  <=  sqrt(opts(3))
%
% 2.  If opts(6)=1, Andy Grace's criteria is used:
%                     ms=.5*max(abs(s))  <   opts(2)
%                      AG=.5*abs(fp'*s)  <   opts(3)
%                                max(g)  <   opts(4)
%
% 3.  If opts(6)~=(-1) & opts(6)~=1, the default criteria is used:
%                          max(abs(dx))  <=  opts(2)
%                                   KTO  <=  opts(3)
%                                max(g)  <=  opts(4)
%
% 4.  If opts(6)==2, add Slowed convergence criterion to (3) above.
%                    NLG = norm(Lagrangian gradient) instead of KTO
%---------------------------------------------------------------------
% Explain trouble shooting information
%
% If opts(1)=2 the following information will also be displayed
% when applicable:
%
%     'dH' - Hessian has been perturbed for improved conditioning
%     'aS' - The augmented Lagrangian type Search direction was used
%     'mS' - The modified Search direction problem was used
%     'sx' - Design variables are being scaled
%     'sf' - Objective function is being scaled
%     'sg' - One or more constraint functions scaled
%     'L1' - L1 merit (exterior penalty)
%--------------------------------------------------------------------------
% Copyright (c) 2015, Robert A. Canfield, Mark Spillman. All rights reserved.
%                     Virginia Tech and Air Force Institute of Technology
%                     bob.canfield@vt.edu
%                    <http://www.aoe.vt.edu/people/faculty/canfield.html>
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal with the Software without restriction, including without
% limitation the rights to use, copy, modify, merge, publish, distribute,
% sublicense, and/or sell copies of the Software, and to permit persons
% to whom the Software is furnished to do so, subject to the following
% conditions:
%
% * Redistributions of source code must retain the above copyright notice,
%   this list of conditions and the following disclaimers.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimers in the
%   documentation and/or other materials provided with the distribution.
%
% * Neither the names of Robert A. Canfield, Virginia Tech, Mark Spillman,
%   Air Force Institute of Technology, nor the names of its contributors
%   may be used to endorse or promote products derived from this Software
%   without specific prior written permission.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
% THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
%--------------------------------------------------------------------------
type LICENSE.txt

