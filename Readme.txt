UETC 2013
Joanes Lizarraga
Matlab scripts for UETC post-production analysis.


First of all, set the path to the directory where UETC data is stored, using a global variable.
e.g: gpath = ‘UETCdata/’ Do not forget the final ‘/’
All of them need plotc.m (in order to plot eigenvector dot products), multiPlot.m and multiPlotZoom.m.

/////////////////////////////////
1. Stats file reader: statsFile.m
/////////////////////////////////
Usage: tOffset = statsFile(pNum, id, run, tFit, dx, N, path)

Outputs:
• Plots: Energy, field modulus, scaling...
• Mean scaling slope, both from lagrangian and winding.
• Mean tOffset

	e.g: tOffsetmat = (-1,’mat_%2‘, 1:3, [500 1000], 0.5, 4094, gpath)

/////////////////////////////////
2. Plotting:
/////////////////////////////////
%%%%%%%%%%%%%%%%%%
2.1. ETCplot.m
%%%%%%%%%%%%%%%%%%
Usage: ETCplot(CName, id, run, tRef, tOffset, tLimit, path, parent)

Output: Individual ETC plot
	e.g: ETCplot(‘scalar11‘,’mat_2%’,1:3,199.930,toffsetmat,[150 1000], gpath)

Dependences:
	ETCload.m
	Store kt, r, Correlators and calculates the standard deviation of the correlators.

%%%%%%%%%%%%%%%%%%
2.2. UETCplot.m
%%%%%%%%%%%%%%%%%%
Usage: UETCplot(CName, id, run, tRef, path, parent)

Output: - Individual UETC plot
	e.g: UETCplot(‘scalar11‘,’mat_%2‘,1:3,199.930,toffsetmat,gpath)
	- Coherence functions, just adding ‘*CName’.
	e.g: UETCplot(‘*scalar11‘,’mat_%2‘,1:3,199.930,toffsetmat,gpath)

Dependences:
	• UETCload.m
	Used to store kt, r and original Correlators from the simulation outputs. If we
	consider more than one run, it averages the Correlators.
	Usage: [kt,r,C] = UETCload(filenamePre, CName, id, run, tRef, tOffset)
		- filenamePre = gpath
		- kt and r using UETCinfoRead.m
		- If tOffset is set, performs time offset correction
	• UETCtimeOffset.m
	Used to correct kt, r and Correlators, considering the time offset. This version has
	correct vector time offset.

/////////////////////////////////
3. Eigenvector decomposition: UETCeigen.m
/////////////////////////////////
Usage: UETCeigen(id, run, tRef, tOffset, inPath, outPath, era, number, bootstrap, Ni, before, rMax,
plateauMode, highKtMode,plotLevel)
- inPath = gpath
- outPath relative to inPath, usually outPath = ‘’

Summary:
1. Depending on the value of Ni set the dimensions of the interpolated matrix, Ni x Ni. Usually kt’s separated linearly.
2. From original UETC’s, where data is represented using k*sqrt(tt’) and r, rewrite UETC’s by kt and kt’. In order to accomplish this, interpolate UETC’s into a Ni x Ni matrix.
3. Average over number of runs.
4. Obtain the eigenvectors and eigenvalues.
5. Output interpolated eigenvectors of dimension (8*Ni for scalars and 4*Ni for vector and tensor)

Additional options:
1. Output non-interpolated eigenvectors, before=1.
2. Extrapolate data at high kt, using ETC’s.
3. Addition of extra data at low kt, plateau mode.
4. Plots
	3.1. Firts eigenvalues and eigenvectors
	3.2. Reconstruct ETC’s from eigenvectors, convergence
	3.3. UETC reconstruction


Internal functions:
1. eigen(): Sort eigenvalues and eigenvectors in descending eigenvalue.
2. UETCktr2xx(): Map kt and r into KT and R, constructed from x(Ni).
3. UETCmatrix(): Map UETC’s into the matrix defined by KT and R.
4. UETCextrap(): Extrapolation at high kt, usually not used.
5. output(): Create output files for interpolated eigenvectors


/////////////////////////////////
4. Eigenvector post analysis
/////////////////////////////////
%%%%%%%%%%%%%%%%%%
4.1 UETCoutputCMBeasy.m JO 2013
%%%%%%%%%%%%%%%%%%
Usage: UETCoutputCMBeasy(id, inPath, outPath, Ni, Evalmin, Evaldim, interp, plot, sign,
change)

Summary:
1. Read non-interpolated (Before) eigenvectors and eigenvalues.
2. If sign=1, check and correct radiation vs matter eigenvector relative sign.
3. If interp=1, interpolate eigenvectors using more kt points. Initially non-interpolated eigenvectors are represented by Ni dimensional vectors, after interpolation they are represented
by (8*Ni or 4*Ni) dimensional vectors. Output interpolated eigenvectors, number of eigenvector in
output files determined by Evaldim.
4. Plot rad vs mat eigenvectors dot product from Evalmin to Evaldim. Useful to determine whether eigenvector order is mixed or is correct.

%%%%%%%%%%%%%%%%%%
4.2 UETCeigenCheck.m
%%%%%%%%%%%%%%%%%%
Usage: UETCeigenCheck(path, ID, lim, Ni, BefAf)

Summary:
1. BefAf=0 Before interpolation. BefAF=1 After interpolation.
2. Check rad vs mat eigenvector corrrespondences:
2.1. Plotting eigenvectors and also checking correspondences with previous and posterior eigenvectors.
2.2. Plotting 2D eigenvector dot product surfaces.

