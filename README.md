CHQuant is a protocol for analysing molecular dynamics (MD) trajectories using convex hull (CH) volumes to quantify conformational space coverage.

This package can:

- Read in an MD trajectory, find the convex hulls for each atom and save the convex hull vertices 
- Generate Pandas DataFrames containing convex hull data (convex hull volumes, surface areas, percentage contributions) for a molecular dynamics trajectory
- Check whether atoms in a functional group have overlapping convex hulls
- Calculate convex hull volumes for functional groups

Using additional scripts, the outputs of this package can be used to:
- Generate mist plots with atomic convex hulls overlaid
- Generate "convex hull" plots where atoms are coloured based on atomic convex hull volume.

The additional scripts require matplotlib, scikit-learn, plotly, and ICHOR.

The tutorials.ipynb file contains examples showing how to use the package.