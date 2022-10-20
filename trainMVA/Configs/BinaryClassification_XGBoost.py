# In this file you can specify the training configuration

#####################################################################
######Do not touch this
import numpy as np
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dropout
from tensorflow.keras.callbacks import EarlyStopping
#####################################################################
####Start here
#####################################################################
OutputDirName = 'BinaryClassification_XGBoost_barrel' #All plots, models, config file will be stored here
branches=['scE',
          'scEt',
          'scEta',
          'scPhi',
          'scEtaWidth',
          'scMinDrWithGsfElectornSC_',
          'scFoundGsfMatch_',
          'scR9',
          'scSigmaIetaIeta',
          'scSigmaIetaIphi',
          'scEMaxRatio',
          'scSwissCross',
          'scE2x5_MaxRatio',
          'scE2ndRatio' ,       
          'scPFPhoIso1',
          'scPFPhoIso2',
          'scPFPhoIso3',
          'scPFPhoIso4',
          'scPFPhoIso5',
          'scPFChIso1',
          'scPFChIso2',
          'scPFChIso3',
          'scPFChIso4',
          'scPFChIso5'          
         ]
baseInputPath="/eos/user/r/rchudasa/SWAN_projects/photonID-pp-Run2/makeInputTrees/workarea/"
Debug=True # If True, only a small subset of events/objects are used for either Signal or background #Useful for quick debugging
processes = [


]
#Branches to read #Should be in the root files #Only the read branches can be later used for any purpose

SaveDataFrameCSV,loadfromsaved=True,False #If loadfromsaved=True, dataframe stored in OutputDirName will be read

Classes,ClassColors = ['IsolatedSignal','NonIsolated'],['#377eb8', '#ff7f00']
#Remeber: For binary classification, first class of the Classes argument should be signal, otherwise, plots might not make sense.

processes = [
    {'Class':'Signal','fileName':'mc_signal_bsMMG_barrel.root',
    'treeName':'genMatchedBMMGSCTree', 'selection':'(scEt>4) & (scEt<15) & abs(scEta)<1.4442', 'process':'bsMMG', 'category':0},
    {'Class':'Backg','fileName':'mc_qcd20To30EmEnriched_barrel.root',
    'treeName':'mergedPi0_SCTree', 'selection':'(scEt>4) & (scEt<15) & abs(scEta)<1.4442', 'process':'QCD','category':1}, 
    {'Class':'Backg','fileName':'mc_bkg_flat_pi0_barrel.root',
    'treeName':'mergedPi0_SCTree', 'selection':'(scEt>4) & (scEt<15) & abs(scEta)<1.4442', 'process':'Flat_pi0','category':1}, 
    {'Class':'Backg','fileName':'data_2018D_barrel.root',
    'treeName':'dataAllSCTree', 'selection':'(scEt>4) & (scEt<15) & abs(scEta)<1.4442', 'process':'Data','category':1}     
]


#MVAs to use as a list of dictionaries
MVAs = [
    #can add as many as you like: For MVAtypes XGB and DNN are keywords, so names can be XGB_new, DNN_old etc. 
    #But keep XGB and DNN in the names (That is how the framework identifies which algo to run
    
    {"MVAtype":"XGB_1", #Keyword to identify MVA method.
     "Color":"green", #Plot color for MVA
     "Label":"XGB try1", # label can be anything (this is how you will identify them on plot legends)
     "features":["ele_fbrem", "ele_deltaetain", "ele_deltaphiin", "ele_oldsigmaietaieta", 
                 "ele_oldhe", "ele_ep", "ele_olde15", "ele_eelepout",
                 "ele_kfchi2", "ele_kfhits", "ele_expected_inner_hits","ele_dr03TkSumPt",
                 "ele_dr03EcalRecHitSumEt","ele_dr03HcalTowerSumEt","ele_gsfchi2","scl_eta","ele_pt",
                 'ele_nbrem','ele_deltaetaseed','ele_hadronicOverEm','ele_olde25max','ele_olde55'],
     "feature_bins":[100 for i in range(22)], #same length as features
     #Binning used only for plotting features (should be in the same order as features), does not affect training
     'Scaler':"MinMaxScaler", #Scaling for features before passing to the model training
     'UseGPU':False, #If you have a GPU card, you can turn on this option (CUDA 10.0, Compute Capability 3.5 required)
     "XGBGridSearch":{'min_child_weight': [5]} #All standard XGB parameters supported
    },
]

ptbins = [4,6,8,10,13,15]
etabins = [-1.4442,-1.1,-0.8,-0.5,0,0.5,0.8,1.1,1.4442]
