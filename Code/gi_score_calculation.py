import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib
import matplotlib.pyplot as plt
import sys
import os
plt.switch_backend('agg')
almost_black = '#262626'
cdict = {'red':((0.0,0.0,0.0),
				(0.5,0.0,0.0),
				(1.0,1.0,1.0)),
		'green':((0.0,0.0,0.0),
				(0.5,0.0,0.0),
				(1.0,1.0,1.0)),
		'blue': ((0.0,1.0,1.0),
				(0.5,0.0,0.0),
				(1.0,0.0,0.0))}
blue_yellow = matplotlib.colors.LinearSegmentedColormap('BlueYellow',cdict)
plt.register_cmap(cmap=blue_yellow)


def generate_fc_matrix(matrix_inf):
	#fc_matrix = pd.read_csv("Puro12_Lif12_counts_prop_logfc_matrix.txt",sep="\t")
	fc_matrix = pd.read_csv(matrix_inf,sep="\t")
	singlesTable = pd.DataFrame([s.split(';')[0] for s in fc_matrix.index], index=fc_matrix.index, columns=['gene'])
	single_fc = pd.concat((fc_matrix.loc[singlesTable['gene'] == "negative_control", :].apply(np.nanmean, axis=0),
						fc_matrix.loc[singlesTable['gene'] == "negative_control", :].apply(np.nanstd, axis=0),
						fc_matrix.loc[:,singlesTable['gene'] == "negative_control"].apply(np.nanmean, axis=1),
						fc_matrix.loc[:,singlesTable['gene'] == "negative_control"].apply(np.nanstd, axis=1)),
						axis=1,keys=['b.mean', 'b.std', 'a.mean', 'a.std'])
	return fc_matrix, singlesTable, single_fc

def abbaAverageDepletionScore(fc_matrix, singlesTable):
	fc_matrix_abba = (fc_matrix + fc_matrix.T) / 2

	single_fc_abba = pd.concat((fc_matrix_abba.loc[singlesTable['gene'] == 'negative_control',:].apply(np.nanmean, axis=0), 
								  fc_matrix_abba.loc[singlesTable['gene'] == 'negative_control',:].apply(np.nanstd, axis=0), 
								  fc_matrix_abba.loc[:, singlesTable['gene'] == 'negative_control'].apply(np.nanmean, axis=1),
								  fc_matrix_abba.loc[:, singlesTable['gene'] == 'negative_control'].apply(np.nanstd, axis=1)), 
								 axis=1, keys=['b.mean','b.std','a.mean','a.std'])

	return fc_matrix_abba, single_fc_abba


def getXYB(sgRNA, single_fc, fc_matrix, variablePosition, fixedPosition, returnXerr=False):
	if not returnXerr:
		return single_fc[variablePosition+'.mean'], \
			fc_matrix.loc[sgRNA,:] if fixedPosition == 'a' else fc_matrix.loc[:,sgRNA], \
			single_fc.loc[sgRNA, fixedPosition +'.mean']
	else:
		return single_fc[variablePosition+'.mean'], \
			fc_matrix.loc[sgRNA,:] if fixedPosition == 'a' else fc_matrix.loc[:,sgRNA], \
			single_fc.loc[sgRNA, fixedPosition +'.mean'], single_fc[variablePosition+'.std']
	


#calculate epistasis interactions, optionally z-standardizing based on negative controls
def calculateInteractions(fc_matrix, single_fc, singlesTable, fitFunction, zstandardize=True):
	emap1 = pd.DataFrame(np.zeros(fc_matrix.shape), index=fc_matrix.index, columns=fc_matrix.columns)
	variablePosition, fixedPosition = 'a','b'
	for i, sgRNA in enumerate(fc_matrix.index):
		xdata, ydata, bdata = getXYB(sgRNA, single_fc, fc_matrix, variablePosition, fixedPosition)
		fit = fitFunction(xdata, ydata, bdata)
		###reverse the deviation
		epistasis = -(ydata - fit(xdata))

		if zstandardize:
			emap1.loc[sgRNA,:] = epistasis / epistasis.loc[singlesTable['gene'] == 'negative_control'].std()
		else:
			emap1.loc[sgRNA,:] = epistasis 

	emap2 = pd.DataFrame(np.zeros(fc_matrix.shape), index=fc_matrix.index, columns=fc_matrix.columns)
	variablePosition, fixedPosition = 'b','a'
	for i, sgRNA in enumerate(fc_matrix.index):
		xdata, ydata, bdata = getXYB(sgRNA, single_fc, fc_matrix, variablePosition, fixedPosition)
		
		fit = fitFunction(xdata, ydata, bdata)
		###reverse the deviation
		epistasis = -(ydata - fit(xdata))

		if zstandardize:
			emap2.loc[sgRNA,:] = epistasis / epistasis.loc[singlesTable['gene'] == 'negative_control'].std()
		else:
			emap2.loc[sgRNA,:] = epistasis 

	emap12 = (emap1+emap2)/2
	
	emap_ave = (emap12 + emap12.T) / 2
	
	return emap1, emap2, emap_ave

### fit functions for calculating interactions and plotting
def linearFitForceIntercept(xdata, ydata, bdata):
	m1 = optimize.fmin(lambda m, x, y: ((m*x + bdata - y)**2).sum(), x0=0.1, args=(xdata, ydata), disp=0)[0]
	return lambda x1: m1*np.array(x1) + bdata

def plotSingleVsDouble_abba(sgRNA, fc_matrix, single_fc, emap, fitFunction=None, showXerr=True):
	fig, axes = plt.subplots(1,2,figsize=(5,6), gridspec_kw={'width_ratios':[6,1]})
	variablePosition, fixedPosition = 'a','b'
	xdata1, ydata1, bdata1, xerr1 = getXYB(sgRNA, single_fc, fc_matrix, variablePosition, fixedPosition, True)
	variablePosition, fixedPosition = 'b','a'
	xdata2, ydata2, bdata2, xerr2 = getXYB(sgRNA, single_fc, fc_matrix, variablePosition, fixedPosition, True)
	
	minVal = np.nanmin([np.nanmin(xdata1), np.nanmin(ydata1), np.nanmin(xdata2), np.nanmin(ydata2)])
	maxVal = np.nanmax([np.nanmax(xdata1), np.nanmax(ydata1), np.nanmax(xdata2), np.nanmax(ydata2)])
	
	#minVal *= 1.1
	#maxVal *= 1.5
	minVal = -3
	maxVal = 3
	axis = axes[0]
	
	if showXerr:
		axis.errorbar(xdata1, ydata1, xerr=xerr1, fmt='none', ecolor=almost_black, alpha=.1, lw=.5, capsize=2, zorder=1)
		
	colorRange = np.percentile(np.abs(np.reshape(emap.values, (len(emap)**2,))), 90)
	print(-1*colorRange)
	result = axis.scatter(xdata1, ydata1, c = emap.loc[:,sgRNA], s=8, alpha=.75, edgecolor = almost_black, cmap='RedGrey', vmin=-4*colorRange, vmax=colorRange)
	
	if fitFunction:
		fit = fitFunction(xdata1, ydata1, bdata1)
		fitrange = np.linspace(minVal,maxVal, 100)
		axis.plot(fitrange, fit(fitrange), color='#808000')

	axis.plot((minVal,maxVal), (minVal,maxVal), color='#BFBFBF',lw=.5)
	axis.plot((minVal,maxVal), (0,0), color='#BFBFBF',lw=.5)
	axis.plot((0,0), (minVal,maxVal), color='#BFBFBF',lw=.5)
	axis.plot((minVal,maxVal), (bdata1,bdata1), color=almost_black, lw=.5)
	axis.set_xlim((minVal, maxVal))
	axis.set_ylim((minVal, maxVal))
	
	axis.set_xlabel('sgRNA single phenotype', fontsize=8)
	axis.set_ylabel('sgRNA double phenotype with ' + sgRNA, fontsize=8)

	axis.xaxis.set_tick_params(labelsize=8)
	axis.yaxis.set_tick_params(labelsize=8)
	axis.yaxis.tick_left()
	axis.xaxis.tick_bottom()
	
	axis = axes[1]
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.spines['left'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.set_xticks([])
	axis.set_yticks([])
	
	cbar = fig.colorbar(result, ax=axis)
	cbar.ax.tick_params(labelsize=8)
	
	plt.tight_layout()
	
	return fig

def main():
	matrix_inf = sys.argv[1] ###each sgRNA pairs log2(fc) enrichment between d0 and d30
	out_name = sys.argv[2]  ###output the matrix of epistasis interaction
	fc_matrix, singlesTable, single_fc = generate_fc_matrix(matrix_inf)
	fc_matrix_abba, single_fc_abba = abbaAverageDepletionScore(fc_matrix, singlesTable)
	#emap1_l,emap2_l,emap_ave_l = calculateInteractions(fc_matrix, single_fc, singlesTable, linearFitForceIntercept)
	emap1_abba_l,emap2_abba_l,emap_ave_abba_l = calculateInteractions(fc_matrix_abba, single_fc_abba, singlesTable, linearFitForceIntercept)
	emap_ave_abba_l.to_csv(out_name + "_emap_ave_abba_l.txt", sep='\t')
	#fig1 = plotSingleVsDouble("Myc_Enhancer_3;CTGCCCAC", fc_matrix_abba, single_fc_abba, emap_ave_abba_l, fitFunction=linearFitForceIntercept,showXerr=False)
	#fig1.savefig('Myc_Enhancer_3;CTGCCCAC_l.pdf',format='pdf')
	

if __name__ == "__main__":
	main()
