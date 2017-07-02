# ofli-OH-trapping
Calculate Orthogonal Fast Lyapunov Indicators for OH molecule traps

Following up on some work by Rosario Gonzalez-Ferez et al in a 2014 PRE entitled "Nonlinear dynamics of atoms in a crossed optical dipole trap", I wish to calculate chaos indicators for traps used for OH molecules here in Jun Ye's lab in JILA, University of Colorado / NIST.

One of the biggest challenges has been representing the traps in an analytical, exactly differentiable way. Last summer, 10 months ago, I worked on this problem extensively, with the approach of fitting the traps with an analytic expansion which could then be differentiated appropriately. I eventually found a suitable expansion, but the resulting expression still seemed beyond MATLAB's numerical capabilities.

Now I'm wondering if it wouldn't be possible to simply represent the traps as cubic or even higher order splines. They would be differentiable to the required order and fully general for all trapping configurations. We'll see.

For now, I want to try to get this working on GIThub, because I'm tired of dropping a project and then having to search home, lab, shares, clusters and email in order to pick up where I left off. I need to see if I can run a git on the cluster, since that will be the limiting factor in the usefulness of the technique.

# update 6/20/2017
It’s working on github. I’ve got a spline version that is close to working. A few months ago it was working but giving blurred results, and just yesterday I discovered it’s because I didn’t use a long enough propagation time. (10 instead of 1000 seconds). I’m not sure of this, because I’m still waiting for the longer sims to finish. (1 day so far, expected time… 1 week?).

In the meanwhile I got a spline of the pin trap and submitted that also. We’ll see how it looks. More of a guess on the relevant energy scale.

# update 7/1/2017
The spline version of the crossed dipole trap works great. The pin version wasn't finishing though. I implemented mirror symmetry a few days ago and the crossed dipole finished even faster but the pin was still stuck. Just this evening I realized that the pin was getting totally bogged down when molecules were falling out of the trap. I think the potential was quintically diverging so the velocity was going nuts and the time step got super slow. In comparison the x-dipole doesn't get slow when the particle escapes because the potential doesn't diverge.

I addressed this by implementing MATLAB's builtin event handling functionality for ODE solving. I've re-submitted the pin jobs, and I'm confident of better results. (But not confident that there won't be some small bug that requires me to resubmit another few times yet)

In the meanwhile, what else does this project need to become a paper?
-Obviously super pretty pictures from all slices of the pin trap at various energies.
-Can I take any hard data or ask any chaos related questions?
-Is there a trap with chaos at arbitrarily low energies? What would this do to evaporation or BEC?