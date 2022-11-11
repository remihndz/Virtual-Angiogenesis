# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys, os
from math import acos
from scipy.linalg import solve
import vtk


from utils import *

sns.set_theme(context='talk', style='darkgrid')
plt.style.use('fivethirtyeight')

# %%
dataDir    = '../Data/'
resultsDir = './'

cols = ['Flow', 'Pressure', 'FlowVelocity', 'PressureDrop',
        'Boundary', 'Viscosity', 'WSR', 'WSS', 'Radius', 'Length', 
        'Bifurcation', 'Stage', 
        'Diameter', 'Resistance',
        'Simulation']
TotalDF = pd.DataFrame([], columns=cols)


for k, CCOFile in enumerate(sys.argv[1:]):

    imagesDir  = resultsDir + CCOFile.split('/')[-1][:-4]+'/'

    try:
        # Create target Directory
        os.mkdir(imagesDir)
        print("Directory " , imagesDir,  " Created ") 
    except FileExistsError:
        print("Directory " , imagesDir,  " already exists")

    # %%
    C, Radii, Length, D, df = ReadCCO(CCOFile, project=True)

    # %%
    df = SolveFlow(C, Radii, Length, D, df, viscosityModel='Pries', Qperf=15)
    df = df.assign(Simulation = k)
    TotalDF = pd.concat([df[cols], TotalDF])

    # %%
    '''
    Units should be:
        - Length (and volumes) in microns
        - Time in seconds
        - Pressure in mmHg
        - Wall shear rate in s^-1
        - Wall shear stress in dyn.cm^2 (=0.1Pa)
    '''


    scale=False
    if scale:
        df.Flow = df.Flow/df.Flow.max()
        
        
    fig, axes = plt.subplots(3,2, figsize=(16,12),sharex='none', sharey='none')
    #fig.tight_layout()
    axes = axes.ravel()
    for ax in axes:
        ax.set(yscale='linear')
    sns.scatterplot(data=df, x='Radius', y='Flow', hue='Boundary', ax=axes[0])
    sns.scatterplot(data=df, x='Radius', y='Pressure', hue='Boundary', ax=axes[2])
    sns.scatterplot(data=df, x='Bifurcation', y='Diameter', hue='Boundary', ax=axes[1])
    #sns.scatterplot(data=df, x='Radius', y='Pressure', hue='Boundary', ax=axes[3])
    axes[3].set(yscale='linear')
    sns.histplot(data=df, y='PressureDrop', x='Bifurcation', discrete=(True, False), cbar=True, ax=axes[3])

    sns.scatterplot(data=df, x='Radius', y='WSR', hue='Bifurcation', ax=axes[4])
    sns.scatterplot(data=df, x='Radius', y='WSS', hue='Bifurcation', ax=axes[5])

    fig.tight_layout()
    fig.savefig(imagesDir + 'Haemodynamics.jpg', dpi=300)

    # %% [markdown]
    # # Comparing 1<sup>st</sup> and 2<sup>nd</sup> order arterioles
    # Data from https://iovs.arvojournals.org/article.aspx?articleid=2182855

    # %%
    ## !!! The data has been scaled to match the synthetic vasculature's values
    dataVisc = pd.read_csv(dataDir + 'Viscosity-vs-ShearRate.csv', skipinitialspace=True)
    dataVisc.loc[dataVisc.Study == 'Nagaoka (2006)', 'Viscosity'] = dataVisc.loc[dataVisc.Study=='Nagaoka (2006)', 'Viscosity']/dataVisc['Viscosity'].max()*df['Viscosity'].max()
    dataVisc['Study'] = 'Nagaoka (2006)'

    dataVisc = pd.concat((dataVisc,
                        pd.DataFrame([[wsr, visc, 'Synthethic Vasculature'] 
                                    for wsr,visc in zip(df.WSR.values, df.Viscosity.values)], 
                                    columns=dataVisc.columns)))

    fig, ax = plt.subplots(figsize=(10,6))
    sns.scatterplot(data=dataVisc, x='WSR', y='Viscosity', hue='Study', ax=ax)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

    fig.tight_layout()
    fig.savefig(imagesDir + 'Viscosity-vs-WSR.jpg', dpi=300)

    # %%
    dataHaemodynamics = pd.read_csv(dataDir + 'HaemodynamicsFirstAndSecondOrder.csv', skipinitialspace=True)
    dataToPlot = ['Diameter', 'FlowVelocity', 'Flow', 'Resistance', 'WSR', 'Viscosity', 'WSS']

    # for measurementData, measurementSynth in zip(dataHaemodynamics.columns.values[2:], dataToPlot):
    #     print(measurementData)
    #     print('\tOrder\tData\tSynthetic vessels')
    #     for bif in [1,2]:
    #         data = dataHaemodynamics.loc[(dataHaemodynamics.Bifurcation==bif-1)
    #                                     & (dataHaemodynamics['Vessel type']=='Artery')][measurementData].values[0]
    #         synth = df.loc[df.Bifurcation==bif][measurementSynth].mean()
    #         nVessels = df.loc[df.Bifurcation==bif][measurementSynth].shape[0]
    #         print(f'\t{bif+1}\t{data}\t{synth} +- {df.loc[df.Bifurcation==bif][measurementSynth].std()} (mean of {nVessels} vessels)')

    # %% [markdown]
    # # Blood flow 
    # Data from:
    # 
    #     - https://opg.optica.org/boe/fulltext.cfm?uri=boe-5-2-630&id=277856 (FD Dopler OCT)
    #     - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4465795/ (fluorescent RBC - rats)
    #     - https://iovs.arvojournals.org/article.aspx?articleid=2159754 (bidirectional Dopler velocimetry + fundus photograph)
    #     - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9456494/ (Dopler OCT)
    #     - https://diabetesjournals.org/diabetes/article-abstract/doi/10.2337/db22-0219/147603/Retinal-Oxygen-Metabolism-in-Patients-with-Type-II?redirectedFrom=fulltext (Dopler OCT)
    #     - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8525846/ (Doplet OCT)

    # %%
    totalRBF = pd.read_csv(dataDir+'TotalRetinalBloodFlow.csv', skipinitialspace=True)
    # totalRBF.at[totalRBF.shape[0]] = [df.Flow.sum(), 'Synthetic Vasculature (sum of all flows)', np.nan]
    totalRBF.loc[totalRBF.shape[0]] = [df.Flow.max(), 'Synthetic Vasculature', np.nan]

    fig, ax = plt.subplots(figsize=(8,4))
    for i, study in enumerate(totalRBF.Study.value_counts().index.to_list()):
        mean = totalRBF.loc[totalRBF.Study==study].mean(numeric_only=True).values
        std  = totalRBF.loc[totalRBF.Study==study].std(numeric_only=True).values
        ax.bar(i, mean[0], yerr=std[0], label=study)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.set_xticks([])
    ax.set_ylabel('Total RBF (microliter/min)')

    fig.tight_layout()
    fig.savefig(imagesDir + 'TotalRetinalBF.jpg', dpi=300)

    # %%
    dataFlow = pd.read_csv(dataDir+'Flow-vs-Diameter.csv', skipinitialspace=True)
    dataFlow = pd.concat((dataFlow.loc[dataFlow.Study!='Kornfield (2015)'],
                        pd.DataFrame([[d,f,'Synthethic Vasculature']
                                for d,f in zip(df.Diameter.values, df.Flow.values)], columns=dataFlow.columns)),
                        ignore_index=True)
    CFD = pd.read_csv(dataDir + 'CFD_Liu_2008.csv', skipinitialspace=True)

    fig, axes = plt.subplots(1,2, figsize=(16,6))
    axes = axes.ravel()

    sns.scatterplot(data=dataFlow, x=dataFlow.columns[0], y=dataFlow.columns[1], hue='Study', ax=axes[0])
    #axes[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

    sns.scatterplot(data=CFD, x='Diameter (microns)', y='Ratio of inlet flow', ax=axes[1])
    sns.scatterplot(x=df['Diameter'].to_numpy(), y=df['Flow']/df['Flow'].max(), ax=axes[1])
    axes[1].legend(['Liu (2008) CFD simulations', 'This work'])

    fig.tight_layout()
    fig.savefig(imagesDir + 'Flow_vs_Diameter.jpg', dpi=300)

    # %% [markdown]
    # # Branching analysis
    # Data from:
    # 
    #     - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2215753/

    # %%
    # Compute the branching angles and caliber ratios

    dldp, dsdl, dsdp = [], [], []
    theta_s, theta_l = [], []

    for idx, v in df.loc[df.Boundary!=-1].iterrows():

        branchesIdx = v.BranchesId
        branchesRadius = df.loc[branchesIdx, 'Radius'].to_numpy()         
        dldp.append(branchesRadius.max()/v.Radius)
        dsdp.append(branchesRadius.min()/v.Radius)
        dsdl.append(branchesRadius.min()/branchesRadius.max())    

        largerDescendentIdx  = branchesIdx[np.argmax(branchesRadius)]
        smallerDescendentIdx = branchesIdx[np.argmin(branchesRadius)]

        parentSegment  = (df.at[idx, 'xDist'] - df.at[idx, 'xProx'])/df.at[idx, 'Length']
        largerSegment  = ((df.at[largerDescendentIdx, 'xDist'] - df.at[largerDescendentIdx, 'xProx'])
                        /df.at[largerDescendentIdx, 'Length'])
        smallerSegment = ((df.at[smallerDescendentIdx, 'xDist'] - df.at[smallerDescendentIdx, 'xProx'])
                        /df.at[smallerDescendentIdx, 'Length'])

        # Check we have normalized properly
        if np.linalg.norm(parentSegment) != 1.0:
            parentSegment /= np.linalg.norm(parentSegment)
        if np.linalg.norm(largerSegment) != 1.0:
            largerSegment /= np.linalg.norm(largerSegment)
        if np.linalg.norm(smallerSegment) != 1.0:
            smallerSegment /= np.linalg.norm(smallerSegment)
        
        # Sometimes, a 'branch' is just the continuation of a vessel and has the same direction (theta < 1e-3 rad)
        dotProdLarger  = parentSegment.dot(smallerSegment)
        dotProdSmaller = parentSegment.dot(largerSegment)
        
        np.seterr(all='raise')
        try:
            theta = np.arccos(dotProdLarger)
        except FloatingPointError:
            theta = np.nan
        if theta < 1e-3:
            theta = np.nan
        theta_l.append(theta)

        try:
            theta = np.arccos(dotProdSmaller)
        except FloatingPointError:
            theta = np.nan
        if theta < 1e-3:
            theta = np.nan
        theta_s.append(theta) 

    # df['dsdl'] = dsdl
    # df['dldp'] = dldp
    # df['ds_dp'] = dsdp
    # df['theta larger branch'] = theta_l
    # df['theta smaller branch'] = theta_s
 
    branching = pd.read_csv(dataDir+'Branching.csv', skipinitialspace=True)
    branching = pd.concat((branching, 
                        pd.DataFrame([[ds_dl, dl_dp, ds_dl, ds_dp, ds_dl, 
                                        thetal, ds_dl, thetas, 'Synthetic Vasculature']
                                        for ds_dl, dl_dp,
                                        ds_dp, thetal, thetas in zip(dsdl, dldp, dsdp, theta_l, theta_s)], 
                                    columns=branching.columns.values)))

    # Sort the data into bins for plotting
    branching['bins_thetal'] = pd.cut(branching['ds/dl_theta_l'], 50).astype(str)
    branching['bins_thetas'] = pd.cut(branching['ds/dl_theta_s'], 50).astype(str)

    # %%
    fig, axes = plt.subplots(2,2, figsize=(16,12), sharex='none', sharey='none')
    axes = axes.ravel()

    #sns.scatterplot(data=df, x='Radius', y='Length', ax=ax[0])
    alpha = 0.3
    # sns.scatterplot(data=branching, x='ds/dl_theta_l', y='theta_l', hue='Study',
    #              ax=axes[0], legend=False)
    # sns.scatterplot(data=branching, x='ds/dl_theta_s', y='theta_s', hue='Study',
    #              ax=axes[1], legend=False)

    axes[0].scatter(branching.loc[branching.Study!='Synthetic Vasculature', 'ds/dl_theta_l'].to_numpy(),
                    branching.loc[branching.Study!='Synthetic Vasculature', 'theta_l'].to_numpy(), alpha=1)
    axes[0].scatter(branching.loc[branching.Study=='Synthetic Vasculature', 'ds/dl_theta_l'].to_numpy(),
                    branching.loc[branching.Study=='Synthetic Vasculature', 'theta_l'].to_numpy(), alpha=alpha)

    axes[1].scatter(branching.loc[branching.Study!='Synthetic Vasculature', 'ds/dl_theta_s'].to_numpy(),
                    branching.loc[branching.Study!='Synthetic Vasculature', 'theta_s'].to_numpy(), alpha=1)
    axes[1].scatter(branching.loc[branching.Study=='Synthetic Vasculature', 'ds/dl_theta_s'].to_numpy(),
                    branching.loc[branching.Study=='Synthetic Vasculature', 'theta_s'].to_numpy(), alpha=alpha)


    ### !!!! Shifts the values of the data to the left axis
    # sns.pointplot(data=branching, x='bins_thetal', y='theta_l', hue='Study', errorbar='sd', 
    #               markers='.', ax=axes[0], legend=False, linestyles='none')
    # sns.pointplot(data=branching, x='bins_thetas', y='theta_s', hue='Study', errorbar='sd', 
    #               markers='.', ax=axes[1], legend=False, linestyles='none')
    # axes[0].get_legend().remove()
    # axes[1].get_legend().remove()
    # axes[0].set_xticklabels('')
    # axes[1].set_xticklabels('')


    sns.scatterplot(data=branching, x='ds/dl_dl/dp', y='dl/dp', hue='Study', ax=axes[2], legend=True)
    axes[2].legend(bbox_to_anchor=(.02, 0), loc='lower left', borderaxespad=0)
    sns.scatterplot(data=branching, x='ds/dl_ds/dp', y='ds/dp', hue='Study', ax=axes[3], legend=False)

    fig.text(0.5, 0.0, 'Branches assymetry ratio', ha='center', va='center', fontsize=20)
    axes[0].set_ylabel(r'$\theta_{larger\,branch}$', fontsize=25)
    axes[1].set_ylabel(r'$\theta_{smaller\,branch}$', fontsize=25)
    axes[2].set_ylabel(r'$\frac{r_{larger\,branch}}{r_{parent}}$', fontsize=25)
    axes[3].set_ylabel(r'$\frac{r_{smaller\,branch}}{r_{parent}}$', fontsize=25)
    for ax in axes:
        ax.set_xlabel('')

    fig.tight_layout()
    fig.savefig(imagesDir + 'BranchingBehaviour.jpg', dpi=300)

    # %%
    fig = plt.figure(figsize=(12,12))
    for idx, v in df.iterrows():
        # Print the connections
        xIn,yIn, _   = v.xProx
        xOut,yOut, _ = v.xDist
        if v.Boundary==1:
            plt.plot([xIn,xOut], [yIn, yOut], c='blue', label='Inlet')
        elif v.Boundary==-1:
            plt.plot([xIn,xOut], [yIn, yOut], c='red', label='Outlet')
        else:
            plt.plot([xIn,xOut], [yIn, yOut], c='black')
            
    plt.legend(['Inlet', 'Outlet'])
    fig.tight_layout()
    fig.savefig(imagesDir + 'Visu.jpg')

    # %%
    df.to_csv(imagesDir + CCOFile.split('/')[-1][:-4] + '.csv')
    # %%

    vtkTree = vtk.vtkPolyData()
    points  = vtk.vtkPoints()
    lines   = vtk.vtkCellArray()
    nodeDataRadius = vtk.vtkDoubleArray()
    cellDataRadius = vtk.vtkDoubleArray()
    cellDataFlow   = vtk.vtkDoubleArray()
    cellDataStage  = vtk.vtkDoubleArray()

    nodeDataRadius.SetName('Radius')
    cellDataRadius.SetName('Radius')
    cellDataFlow.SetName('Flow')
    cellDataStage.SetName('Stage')

    for idx, v in df.iterrows():
        idProx = points.InsertNextPoint(v.xProx)
        nodeDataRadius.InsertNextValue(v.Radius)
        idDist = points.InsertNextPoint(v.xDist)
        nodeDataRadius.InsertNextValue(v.Radius)

        vtkSegment = vtk.vtkLine()
        vtkSegment.GetPointIds().SetId(0, idProx)
        vtkSegment.GetPointIds().SetId(1, idDist)

        lines.InsertNextCell(vtkSegment)
        cellDataFlow.InsertNextValue(v.Flow)
        cellDataRadius.InsertNextValue(v.Radius)
        cellDataStage.InsertNextValue(v.Stage)

    vtkTree.SetPoints(points)
    vtkTree.SetLines(lines)

    vtkTree.GetPointData().AddArray(nodeDataRadius)

    vtkTree.GetCellData().AddArray(cellDataFlow)
    vtkTree.GetCellData().AddArray(cellDataRadius)
    vtkTree.GetCellData().AddArray(cellDataStage)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(imagesDir + CCOFile.split('/')[-1][:-4]+'.vtp')
    writer.SetInputData(vtkTree)
    writer.SetDataModeToBinary()
    writer.Write()

    print("Output writen in " + writer.GetFileName() + '\n')


TotalDF.to_csv(resultsDir + f'{k+1}_simulations.csv')
print('Merged trees saved in ', resultsDir + f'{k+1}_simulations.csv')