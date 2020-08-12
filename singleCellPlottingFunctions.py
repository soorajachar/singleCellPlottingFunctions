#! /usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os,sys,math

from matplotlib import colors,ticker,cm
from scipy.signal import savgol_filter
from operator import itemgetter
from matplotlib.lines import Line2D
import matplotlib as mpl

#For more appropriate umap visualization; might need to install them if using google colab
import datashader as ds
import holoviews as hv
import colorcet as cc
import holoviews.operation.datashader as hd

import warnings
warnings.filterwarnings("ignore")

def returnTicks(xticksToUse):
    logxticks = [-1000,-100,-10,0,10,100,1000,10000,100000]
    logicleXTicks = [64, 212, 229, 231, 233, 251, 399, 684, 925]
    xtickValues = []
    xtickLabels = []
    
    for logxtick in xticksToUse:
        if(logxtick < 0):
            xtickLabels.append('$-10^'+str(int(np.log10(-1*logxtick)))+'$')
        elif(logxtick == 0):
            xtickLabels.append('0')
        else:
            xtickLabels.append('$10^'+str(int(np.log10(logxtick)))+'$')
    
    for tickval in xticksToUse:
        xtickValue = logicleXTicks[logxticks.index(tickval)]
        xtickValues.append(xtickValue)
    
    return xtickValues,xtickLabels

def returnLogYTicksAndLabels(yvals):
    miny = math.floor(yvals)
    maxy = math.ceil(yvals)
    allyticks = list(range(miny,maxy+1))
    allyticklabels = []
    for ytick in allyticks:
        allyticklabels.append('$10^{'+str(ytick)+'}$')
    minoryticks = np.log10(list(range(2,10)))
    allminoryticks = []
    for ytick in allyticks[:-1]:
        for minory in minoryticks:
            allminoryticks.append(ytick+minory)
    return allyticks,allyticklabels,allminoryticks

def is_number(s):
    try:
        complex(s) # for int, long, float and complex
    except ValueError:
        return False

    return True

def facetedSingleCellKDE(data=[],x='',hue='',hue_order='',size='',style='',row='',row_order='',col='',col_order='',col_wrap='',palette='',sharex=True,sharey=True,aspect=1,height=5,smooth_res=27,scaleToMode=False,logScale=False):
    """
    Wrapper for 1D smoothed histograms that follows same keyword conventions as seaborn figure level plots, with additional parameters for plotting single cell flow cytometery data
    
    
    Parameters:
    data (pd.DataFrame):tidy format pandas dataframe
    x (str):Name of dataframe level to plot on x axis
    hue (str):Name of dataframe level to color plot by
    hue_order (list):Fixes order of color assignment
    size (str):Name of dataframe level used to assign line thicknesses to different histograms 
    style (str):Name of dataframe level used to assign line plotting styles to different histograms 
    row (str):Name of dataframe level to assign to each row of subplots
    row_order (list):Fixes order of subplot row assignment
    col (str):Name of dataframe level to assign to each column of subplots
    col_order (list): Fixes order of subplot column assignment
    col_wrap (int):“Wrap” the column variable at this width, so that the column facets span multiple rows. Incompatible with a row facet.
    palette (str or list):Either the name of the list of colors to use, or an explicit list of colors to use in hex
    sharex (bool or str):True will share x axis across all subplots, row will share per row, col will share per column, False will not share at all
    sharey (bool or str):True will share y axis across all subplots, row will share per row, col will share per column, False will not share at all
    aspect (int):Width of figure; 1 is default
    height (int):Height of figure; 5 is default
    smooth_res (int):How much to smooth the histogram (to allow it to resemble a KDE). Lower values cause less smoothing
    scaleToMode (bool):Whether to scale the y axis of each subplot by the maximum count of each plot
    logScale (bool):Whether to scale the y axis of each subplot logarithmically 

    Returns:
    seaborn relplot facetgrid 
    """
    kwargIndices = []
    cols = list(data.columns)
    kwargDict = {}
    for kwargName,kwarg in zip(['hue','row','col','size','style'],[hue,row,col,size,style]):
        if kwarg != '':
            kwargIndices.append(cols.index(kwarg))
            kwargDict[kwargName] = kwarg
    secondaryKwargDict = {}
    for kwargName,kwarg in zip(['row_order','col_order','hue_order','col_wrap','palette','aspect','height'],[row_order,col_order,hue_order,col_wrap,palette,aspect,height]):
        if kwarg != '':
            secondaryKwargDict[kwargName] = kwarg 
        else:
            if kwargName == 'palette':
                if is_number(data[hue].iloc[0]):
                    secondaryKwargDict[kwargName] = sns.color_palette(sns.color_palette(),len(pd.unique(data[hue])))
    secondaryKwargDict['facet_kws'] = {'sharex':sharex,'sharey':sharey}

    kwargIndices = sorted(kwargIndices)
    uniqueKwargCombinations = [list(temp) for temp in set(tuple(temp) for temp in list(data.iloc[:,kwargIndices].values))]
    hist = [0]
    indexTuples = []
    for kwargCombo in uniqueKwargCombinations:
        selectDf = data.copy()
        for kwarg,kwargIndex in zip(kwargCombo,kwargIndices):
            selectDf = selectDf[selectDf[list(selectDf.columns)[kwargIndex]] == kwarg]
        subplotValueList = selectDf[x].values
        newvals = np.append(subplotValueList,[[0,1023]])
        temphist,_ = np.histogram(newvals, bins=256)
        temphist[0]-=1
        temphist[-1]-=1
        if scaleToMode:
            temphist = [temp/max(temphist) for temp in temphist]
        hist+=list(temphist)
        for i in range(len(list(temphist))):
            indexTuples.append(kwargCombo)

    hist = hist[1:]
    numUniquePlots = len(uniqueKwargCombinations)
    histBins = np.tile(list(range(1,1024,4)),numUniquePlots)
    if smooth_res % 2 == 0:
        smooth_res-=1
    smoothedHistBins = savgol_filter(hist, smooth_res, 2) 
    
    if not scaleToMode:
        if logScale:
            cutoff = 1
        else:
            cutoff = 0
    else:
        if logScale:
            cutoff = 0.001 
        else:
            cutoff = 0
    for i,val in enumerate(smoothedHistBins):
        if val < cutoff:
            smoothedHistBins[i] = cutoff
    
    trueNames = itemgetter(*kwargIndices)(list(data.columns))
    if type(trueNames) != list and type(trueNames) != tuple:
         trueNames = [trueNames]
    mi  = pd.MultiIndex.from_tuples(indexTuples,names=trueNames) 
    if scaleToMode:
        hist = [x*100 for x in hist]
        smoothedHistBins = [x*100 for x in smoothedHistBins]
        columns = [x,'% Max']
    else:
        columns = [x,'Count']

    #Construct correct array of data
    dataMatrix = np.matrix([histBins,smoothedHistBins]).T

    newdf = pd.DataFrame(dataMatrix,columns=columns,index=mi)
    data = newdf.reset_index()
    
    fg = sns.relplot(data=data,kind='line',x=x,y=columns[1],**kwargDict,**secondaryKwargDict)

    xtickValues,xtickLabels = returnTicks([-1000,1000,10000,100000])
    maxVal = max(data[x])
    minVal = min(data[x])
    if xtickValues[0] < minVal:
        minVal = xtickValues[0]
    if xtickValues[-1] > maxVal:
        maxVal = xtickValues[-1]
    for i,axis in enumerate(fg.axes.flat):
        axis.set_xticks(xtickValues)
        axis.set_xticklabels(xtickLabels)
        axis.set_xlim([minVal,maxVal])
    if logScale:
        fg.set(yscale='log')

    return fg

def createSingleDataShadedPlot(subsettedData,x='',y='',hue='',hue_order='',palette='',color_key='',label='',group='',spread_threshold=''):
    #No color plotting
    if hue == '':
        if spread_threshold == '':
            dataShadedPlot = hd.datashade(hv.Points(subsettedData, kdims=[x, y]))
        else:
            hd.dynspread(dataShadedPlot = hd.datashade(hv.Points(subsettedData, kdims=[x, y])), threshold=spread_threshold)
    else:
        #Categorical color plotting
        if type(subsettedData[hue].iloc[0]) == str:
            if spread_threshold == '':
                dataShadedPlot = hd.datashade(hv.Points(subsettedData, kdims=[x, y], vdims=[hue]), aggregator=ds.count_cat(hue),color_key=color_key).relabel(label=label,group=group)
            else:
                dataShadedPlot = hd.dynspread(hd.datashade(hv.Points(subsettedData, kdims=[x, y], vdims=[hue]), aggregator=ds.count_cat(hue),color_key=color_key).relabel(label=label,group=group), threshold=spread_threshold)
        #Continuous color plotting
        else:
            if spread_threshold == '':
                dataShadedPlot = hd.datashade(hv.Points(subsettedData, kdims=[x, y]), aggregator=ds.mean(hue),cmap=palette).relabel(label=label,group=group)
            else:
                dataShadedPlot = hd.dynspread(hd.datashade(hv.Points(subsettedData, kdims=[x, y]), aggregator=ds.mean(hue),cmap=palette).relabel(label=label,group=group), threshold=spread_threshold)

    return dataShadedPlot

def facetedSingleCellScatter(data=[],x='',y='',hue='',hue_order='',row='',row_order='',col='',col_order='',col_wrap='',palette='',fig_size=150,context='notebook',biExpXYScale = True, biExpHueScale = False,spread_threshold='',legend=True,dpi=120):
    """
    Wrapper for datashader plots that follows same keyword conventions as seaborn figure level plots, with additional parameters for plotting single cell flow cytometery data
    
    
    Parameters:
    data (pd.DataFrame):tidy format pandas dataframe
    x (str):Name of dataframe level to plot on x axis
    y (str):Name of dataframe level to plot on y axis
    hue (str):Name of dataframe level to color plot by
    hue_order (list):Fixes order of color assignment
    row (str):Name of dataframe level to assign to each row of subplots
    row_order (list):Fixes order of subplot row assignment
    col (str):Name of dataframe level to assign to each column of subplots
    col_order (list): Fixes order of subplot column assignment
    col_wrap (int):“Wrap” the column variable at this width, so that the column facets span multiple rows. Incompatible with a row facet.
    palette (str or list):Either the name of the list of colors to use, or an explicit list of colors to use in hex
    fig_size: (float):Scale the complete figure size by this percentage
    context: (str):Thickness of axis labels and elements. Order from thinnest to thickest is notebook,paper,talk,poster
    biExpXYScale (bool):Whether to use a special biexponential scale on the x/y axes for flow cytometry derived single cell data
    biExpHueScale (bool):Whether to use a special biexponential scale on the hue axis for flow cytometry derived single cell data
    spread_threshold (float):Value between 0 and 1 that determines how visible less dense regions are
    legend (bool):Whether to display legend or not

    Returns:
    matplotlib figure derived from datashaded bokeh figure
    """
    hv.extension("matplotlib")
    hv.output(backend="matplotlib")

    #TODO: Allow color to change even if hue is blank
    #Create palette/color_key
    #No color plotting
    if hue == '':
        color_key = ''
    else:
        #Categorical color plotting
        if type(data[hue].iloc[0]) == str:
            data[hue]=data[hue].astype("category")

            if hue_order == '':
                hue_order = list(pd.unique(data[hue]))
            if palette == '':
                palette = sns.color_palette(sns.color_palette(),len(hue_order)).as_hex()
            else:
                if type(palette) == str:
                    palette = sns.color_palette(palette,len(hue_order)).as_hex()
            color_key = {}
            for i,categoryValue in enumerate(hue_order):
                color_key[categoryValue] = palette[i]
        #Continuous color plotting
        else:
            if palette == '':
                palette = cc.fire
                mplPalette = cc.cm.fire
            else:
                if type(palette) == str:
                    palette = cm.get_cmap(palette, 256)
                mplPalette = palette
            color_key = ''

    shadeList = []
    #Create combinations of all row/column pairs
    if col == '':
        #If no row or column passed
        if row == '':
            dataShadedPlot = createSingleDataShadedPlot(data,x=x,y=y,hue=hue,hue_order=hue_order,palette=palette,color_key=color_key,label='',group='',spread_threshold=spread_threshold)
            shadeList = [dataShadedPlot]
        #If row passed
        else:
            rowValues = data[row]
            if row_order == '':
                row_order = list(pd.unique(rowValues))
            for rowVal in row_order:
                subsettedData = data[data[row] == rowVal]
                dataShadedPlot = createSingleDataShadedPlot(subsettedData,x=x,y=y,hue=hue,hue_order=hue_order,palette=palette,color_key=color_key,label=row+':',group=str(rowVal),spread_threshold=spread_threshold)
                shadeList.append(dataShadedPlot)
    else:
        #If column passed
        if row == '':
            colValues = data[col]
            if col_order == '':
                col_order = list(pd.unique(colValues))
            for colVal in col_order:
                subsettedData = data[data[col] == colVal]
                dataShadedPlot = createSingleDataShadedPlot(subsettedData,x=x,y=y,hue=hue,hue_order=hue_order,palette=palette,color_key=color_key,label=col+':',group=str(colVal),spread_threshold=spread_threshold)
                shadeList.append(dataShadedPlot)
        #If row and column passed
        else:
            rowValues = data[row]
            colValues = data[col]
            if row_order == '':
                row_order = list(pd.unique(rowValues))
            if col_order == '':
                col_order = list(pd.unique(colValues))
            uniquePairValues = list(pd.unique([str(rowValues[i])+','+str(colValues[i]) for i in range(data.shape[0])]))
            for rowVal in row_order:
                for colVal in col_order:
                    if str(rowVal)+','+str(colVal) in uniquePairValues:
                        subsettedData = data[(data[row] == rowVal) & (data[col] == colVal)]
                        dataShadedPlot = createSingleDataShadedPlot(subsettedData,x=x,y=y,hue=hue,hue_order=hue_order,palette=palette,color_key=color_key,label=row+': '+str(rowVal)+',',group=col+': '+str(colVal),spread_threshold=spread_threshold)
                        shadeList.append(dataShadedPlot)

    #Wrap columns appropriately
    if col != '':
        if col_wrap == '':
            subplotsPerRow = len(pd.unique(data[col]))
        else:
            subplotsPerRow = col_wrap
    else:
        subplotsPerRow = 1

    sns.set_context(context)
    hvLayout = hv.Layout(shadeList).cols(subplotsPerRow).opts(tight=True,sublabel_format="",fig_size=fig_size)
    fig = hv.render(hvLayout,dpi=dpi)

    #Make plots the same size (HACK; need to figure out how to make them the same size even with different scales)
    xmin = data[x].min()
    xmax = data[x].max()
    ymin = data[y].min()
    ymax = data[y].max()
    if biExpXYScale:
        ticks = [-1000,100,1000,10000,100000]
        xtickValues,xtickLabels = returnTicks(ticks)
        ytickValues,ytickLabels = returnTicks(ticks)
    for i,axis in enumerate(fig.axes):
        if biExpXYScale:
            axis.set_xticks(xtickValues)
            axis.set_xticklabels(xtickLabels)
            axis.set_yticks(ytickValues)
            axis.set_yticklabels(ytickLabels)
        axis.set_xlim([xmin,xmax])
        axis.set_ylim([ymin,ymax])
        #Also remove common x labels across columns and common y labels across rows
        if col_wrap != '':
            if subplotsPerRow*(math.ceil(len(col_order)/subplotsPerRow)-1) <= i < subplotsPerRow*math.ceil(len(col_order)/subplotsPerRow):
                pass
            else:
                axis.set_xlabel('')
        if row != '':
            if len(row_order) > 1:
                if subplotsPerRow*(len(row_order)-1) <= i < subplotsPerRow*len(row_order):
                    pass
                else:
                    axis.set_xlabel('')
        if col != '':
            if len(col_order) > 1:
                if i % subplotsPerRow == 0:
                    pass
                else:
                    axis.set_ylabel('')
    
    #Create legend
    if legend:
        if color_key == '':
            #No legend needed
            if palette == '':
                pass
            #Colorbar needed (continuous)
            else:
                norm = mpl.colors.Normalize(vmin=min(data[hue].values),vmax = max(data[hue].values))
                fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9)
                # add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
                barWidth = 0.02*(2/subplotsPerRow)
                barHeight = 0.8*(1/math.ceil(len(shadeList)/subplotsPerRow))
                cbar_ax = fig.add_axes([0.95, 0.5-(0.1+barHeight/2)+0.1, barWidth, barHeight])
                if biExpHueScale:
                    hueticks = [-1000,100,1000,10000,100000]
                    hueTickValues,hueTickLabels = returnTicks(hueticks)
                    cbl = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=mplPalette),orientation='vertical',ticks=hueTickValues,cax=cbar_ax)
                    cbl.ax.set_yticklabels(hueTickLabels)
                else:
                    cbl = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=mplPalette),orientation='vertical',cax=cbar_ax)
                cbl.set_label(hue)
        #Regular legend needed (categorical)
        else:
            legendElemList = [Line2D([0], [0], marker='o', color='w', label=hue,markerfacecolor='w', markersize=15)]
            for categoryValue in color_key:
                legend_element = Line2D([0], [0], marker='o', color='w', label=categoryValue,markerfacecolor=color_key[categoryValue], markersize=10)
                legendElemList.append(legend_element)
            fig.legend(handles=legendElemList,frameon=False, loc='center left', bbox_to_anchor=(0.95, 0.5), bbox_transform=fig.transFigure)

    return fig
