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

# update 7/7/17
Some of the pin data finished, but the ofli metric was getting absurdly large (overflowing in some cases, like 10^600 something for matlab!!!). I’m not sure what this means. I’ve been spot checking a few trajectories, and messing with the precision of the ODE solver. I’ve definitely observed some dependence on the precision, one time specifically at a zero crossing of one coordinate, but it wasn’t clear that it was dramatically changing the classification of an orbit as stable, quasiperiodic, or chaotic.

Hypotheses:
1. Even at E=-0.8, essentially all orbits are chaotic, explaining the static-y panel and the large ofli2 values.
2. The oscillation frequency is larger than the crossed dipole trap, thus requiring less integration time to differentiate the orbits.
3. Low precision could be the culprit, preventing any truly quasi-periodic or periodic orbits from showing as such.

# update 8/2/17
Plan of action regarding the above:
1. Tell solver to give up after the OFLI goes above 10^10.
2. Get frame to save some intermediate time data for each frame. Will vastly increase data generated, but will give options for determining best run-time for distinguishing chaotic, quasiperiodic, and periodic.
3. Run a dense grid at E=-0.9.
4. (optional) get things working in another plane, maybe xz plane instead of xy.

# update 8/16/17
I've achieved all of the above except the different plane. Running a dense grid at E=-0.9 and E=-0.8

# update 10/7/17
So the big news… waited forever for the previous things to finish. Finally realized that part of why they were slow was that parfor preallocates and doesn’t redistribute if something finishes first. Addressed this by randomizing the index sent to parfor.

Also implemented intermediate saving of rows of data, which also enables the code to be run on multiple nodes, since they check for a saved version of the row before computing it. Used this to crank through data at E=-0.9 through 0.6.

Noticed that the molecules were all flying out of the trap by -0.6, so revisited the definition of the energy spline. Set E=0 to the trap depth. Before E=-0.7 was the trap depth.

Simultaneously improved the spline by re-exporting the fields from COMSOL with a higher mesh density in the relevant region.

Also added plane capability, although so far I’ve only tested it in the original plane.

CONFIRMED that doing 1e-6 precision on the ODE stepper gives no noticeable change relative to 1e-9 precision. Compare watchprogress to watchprogress6.

Changed the saving to happen for each pixel. This improves computation time especially near the end since a single core doesn’t end up stuck plowing through the most difficult row. Also reduces the overlap when running the same job on many cores. Downside is that the directory on the data drive gets loaded up with files, but who cares? Just don’t try to open it in a GUI.

# update 10/8/17
Analyzed point by point frame. Shows green on axes… normally the axes are points of stability because of the reduced dimensionality presumably. I’m nervous it changed because of the recent change to a single octant spline with numerous “signs” and “abs” to handle the change. To ensure this isn’t hampering the results, I need to setup a side-by-side with and without single octant splining.

# update 10/21/17
I didn't find any blips in the ODE solver associated with the absolute values. Assuming it was working, I let it finish. Took a lot of cores a lot of work. The results suggest that more run time is needed. Its time to investigate a way to reduce the computation time further. In the meantime, I've added a figure file, octant_results_100.fig, and a video watch0p5.avi.