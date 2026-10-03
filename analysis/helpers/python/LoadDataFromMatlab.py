#%%
import math
import matplotlib.pyplot as plt  # for plotting
import random
import numpy as np
from scipy.io import loadmat
import importlib
from functools import partial
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.collections import LineCollection

#  load data
Class_and_function_directory=r"/path/to/python_helpers"

import linear_transforms
import utility_functions as utf
PCA = linear_transforms.PCA

# import sys
# sys.path.append(Class_and_function_directory)
# import linear_transforms
# importlib.reload(linear_transforms)
# PCA=linear_transforms.PCA

class load_cell_SP_posture:
    def __init__(self,
                 dir_BrainRegions=r'/path/to/data/BrainRegions.mat',
                 dir_AllMiceDataFile=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat",
                 dir_ListOfDirectoriesForAllMice=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_for_python.mat"):
        self.dir_BrainRegions = dir_BrainRegions
        self.dir_AllMiceDataFile=dir_AllMiceDataFile
        self.dir_ListOfDirectoriesForAllMice=dir_ListOfDirectoriesForAllMice
        self.already_loaded = False
        self.PlaneNames=np.array([
        "Animal_01 left",  "Animal_01 right",
        "Animal_02 left",  "Animal_02 right",
        "Animal_03 left",  "Animal_03 right",
        "Animal_04 left",  "Animal_04 right",
        "Animal_05 left",  "Animal_05 right",
        "Animal_06 left",  "Animal_06 right",
        "Animal_07 left",  "Animal_07 right",
        "Animal_08 left",  "Animal_08 right",
        "Animal_09 left",  "Animal_09 right",
        "Animal_10 left",  "Animal_10 right"
        ])
    def load_directory_list(self):
        # Path to your modified MAT file
        mat_file = self.dir_ListOfDirectoriesForAllMice
        # Load the MAT file
        matdata = loadmat(mat_file, squeeze_me=True, simplify_cells=True)
        # Assuming your MAT file contains variables 'newallmicedir' and 'allmiceNames'
        self.allmicedir = matdata['char_allmicedir']
        self.allmiceNames = matdata['char_allmiceNames']
        self.already_loaded = True
############################
    def load_data(self,plane,rot_rev_exists=False):
        self.load_directory_list()
        self.load_cellac(self.allmicedir[plane]['cellac'])
        self.load_posture(self.allmicedir[plane]['DLCdata_3DpointsOnly'])
        self.load_posture_info(self.allmicedir[plane]['DLCdata_info'])
        self.load_org(self.allmicedir[plane]['SyncAndSP'],rot_rev_exists=rot_rev_exists)
        self.load_skewness(self.allmicedir[plane]['skewness'])
        self.plane=plane
        if rot_rev_exists:
            self.PlaneNames = np.array([
                "Animal_05 left",
                "Animal_05 right",
                "Animal_06 left",
                "Animal_06 right",
                "Animal_07 left",
                "Animal_07 right",
                "Animal_08 left",
                "Animal_08 right",
            ])
    def load_cellac(self,directory):
        matdata=loadmat(directory)
        cellac=matdata['cellac']    # cellac is a 2D array of shape (n_timepoints, n_cells)
        self.cellac=self.zscore_normalize_along_1stDim(cellac)  # normalize cellac data
    def load_skewness(self,directory):
        matdata=loadmat(directory)
        self.skewness=matdata['skew']
    def load_org(self,directory,rot_rev_exists=False):
        matdata=loadmat(directory)
        org=matdata['org'][0,0]
        self.org=org
        self.f2p_t=org['f2p_t']-1       # 1 indexing to 0 indexing
        self.t_f2p=org['t_f2p']-1       # 1 indexing to 0 indexing
        self.fdlc_t=org['fdlc_t']-1     # 1 indexing to 0 indexing
        self.t_fdlc=org['t_fdlc']-1     # 1 indexing to 0 indexing
        self.t_strial=org['t_strial'].squeeze()-1 # 1 indexing to 0 indexing
        self.t_etrial=org['t_etrial'].squeeze()-1 # 1 indexing to 0 indexing
        self.sp_t=org['sp_t']
        self.cossp_t=org['cossp_t']
        self.sinsp_t=org['sinsp_t']
        self.motionLog=org['motionLog'][0]-1    # 1 indexing to 0 indexing
        self.SPparams_mtypeind=org['SPparams_mtypeind'] 
        self.sp_param_ind=np.array(["theta","phi","tau","x","y","z"])
        if rot_rev_exists:
            self.mtypes= np.array([
            "rot",    "rot_rev",   "tilt",    "yaw",    "zvert",
            "p8_dir1","p8_dir2","p8_dir3",
            "p8_dir4","p8_dir5","p8_dir6",
            "p8_dir7","p8_dir8"
            ])
        else:
            self.mtypes= np.array([
                "rot",    "tilt",    "yaw",    "zvert",
                "p8_dir1","p8_dir2","p8_dir3",
                "p8_dir4","p8_dir5","p8_dir6",
                "p8_dir7","p8_dir8"
            ])
        self.Preprocess_sp_t()
    def Preprocess_sp_t(self):
        new_sp_t=[]
        for i,param in enumerate(self.sp_param_ind):
            if param in ["theta","phi","tau"]:
                temp=np.arctan2(self.sinsp_t[:,i],self.cossp_t[:,i])
                temp=np.mod(temp,2*math.pi)
            else:
                temp=self.sp_t[:,i]
            new_sp_t.append(temp)
        self.mod_sp_t=np.stack(new_sp_t,axis=1)
        # adjust theta and phi
        theta=self.sp_t[:,self.ind_sp_param('theta')]
        phi=self.sp_t[:,self.ind_sp_param('phi')]
        neglog=(phi<0)
        theta[neglog]=-theta[neglog]
        phi[neglog]=-phi[neglog]
        theta=np.mod(theta,2*np.pi)
        self.sp_t[:,self.ind_sp_param('theta')]=theta
        self.sp_t[:,self.ind_sp_param('phi')]=phi
        R,p=self.make_trans_mat_from_sp_t(self.mod_sp_t)
        self.R_t=R
        self.p_t=p
        self.euler_angles_t=self.rotation_matrices_to_euler_angles(R, seq='zyx', degrees=False)
        for trial_ind in range(len(self.motionLog)):
            if np.isnan(self.motionLog[trial_ind]):
                continue
            s,e=self.t_strial[trial_ind].round().astype(int),self.t_etrial[trial_ind].round().astype(int)
            motion_num=self.motionLog[trial_ind].round().astype(int)
            if self.mtypes[motion_num]=='rot':
                theta=self.sp_t[s:e,0]
                theta=utf.force_monotonic_increase(theta)
                self.sp_t[s:e,0]=theta
            if self.mtypes[motion_num]=='rot_rev':
                theta=self.sp_t[s:e,0]
                theta=-utf.force_monotonic_increase(-theta)
                self.sp_t[s:e,0]=theta
            if self.mtypes[motion_num]=='p8':
                theta=self.sp_t[s:e,0]
                mtheta=np.median(theta,axis=0,keepdims=True)
                theta=np.broadcast_to(mtheta,theta.shape)
                self.sp_t[s:e,0]=theta
    
    def make_trans_mat_from_sp_t(self,sp_t):
        theta=sp_t[:,self.ind_sp_param("theta")]
        phi=sp_t[:,self.ind_sp_param("phi")]
        tau=sp_t[:,self.ind_sp_param("tau")]
        x=sp_t[:,self.ind_sp_param("x")]
        y=sp_t[:,self.ind_sp_param("y")]
        z=sp_t[:,self.ind_sp_param("z")]
        zvec=[]
        zvec.append(np.cos(theta)*np.sin(phi))
        zvec.append(np.sin(theta)*np.sin(phi))
        zvec.append(np.cos(phi))
        zvec=np.stack(zvec,axis=1)                  # (n_timepoints, 3)
        pre_x_vec=[]
        pre_x_vec.append(np.cos(phi))
        pre_x_vec.append(np.zeros(theta.shape))
        pre_x_vec.append(-np.cos(theta)*np.sin(phi))
        pre_x_vec=np.stack(pre_x_vec,axis=1)
        x_vec=pre_x_vec/np.linalg.norm(pre_x_vec,axis=1)[:,np.newaxis]
        y_vec=np.cross(zvec,x_vec,axis=1)
        tau=tau[:,np.newaxis]
        xvec=np.cos(tau)*x_vec+np.sin(tau)*y_vec    # (n_timepoints, 3)
        yvec=-np.sin(tau)*x_vec+np.cos(tau)*y_vec   # (n_timepoints, 3)
        R=[xvec,yvec,zvec]
        R=np.stack(R,axis=2)                        # (n_timepoints, 3, 3
        R=np.transpose(R,(1,2,0))   # (3,3,n_timepoints)
        p=[x,y,z]
        p=np.stack(p,axis=1)  # (n_timepoints, 3)
        return R,p
    def rotation_matrices_to_euler_angles(self,R_matrices, seq='zyx', degrees=False):
        """
        Converts a series of 3x3 rotation matrices to Euler angles in a vectorized manner.
        
        Parameters:
            R_matrices (np.ndarray): A numpy array of shape (3, 3, n) where n is the number
                                    of timepoints (rotation matrices).
            seq (str): Euler angle sequence to use (default is 'zyx' corresponding to 
                    yaw-pitch-roll order). Change as needed.
            degrees (bool): If True, returns the angles in degrees. Default is False (radians).
        
        Returns:
            np.ndarray: An array of shape (n, 3) with the Euler angles for each rotation matrix.
        """
        from scipy.spatial.transform import Rotation as R
        # Rearrange the array to shape (n, 3, 3) for vectorized processing.
        R_matrices_vec = np.transpose(R_matrices, (2, 0, 1))
        # Convert all rotation matrices to Euler angles at once.
        euler_angles = R.from_matrix(R_matrices_vec).as_euler(seq, degrees=degrees)
        return euler_angles
    def ind_sp_param(self,sp_param):
        is_scalar = isinstance(sp_param, str)
        sp_param_array = np.atleast_1d(sp_param)
        indices = []
        for param in sp_param_array:
            idx = np.where(self.sp_param_ind == param)[0]
            if idx.size == 0:
                raise ValueError(f"'{param}' not found")
            indices.append(int(idx[0]))
        # Return scalar if input was a scalar, otherwise a NumPy array
        return indices[0] if is_scalar else np.array(indices)
    def load_posture(self,directory):
        matdata=loadmat(directory)
        self.posture=matdata['Q_bp']  # posture is a 2D array of shape (n_timepoints, n_points*3)
        self.catposture=np.concatenate(self.posture[0],axis=1)  # concatenate the (1,14) array to get a (n_timepoints, 42) array
    def load_posture_info(self,directory):
        matdata=loadmat(directory)['rec3d'][0,0]
        self.triError=np.concatenate(matdata['offset'][0,:],axis=1)  # concatenate the (1,14) array to get a (n_timepoints, 14) array
        self.bp_ind= np.array([
            "RightHand", "RightElbow", "RightToe",
            "RightHeel", "RightKnee", "BackU",
            "BackM", "BackL", "Tail", "LeftHand",
            "LeftElbow", "LeftToe", "LeftHeel",
            "LeftKnee"
        ])
        self.edgeList= np.array([
            ["RightHand", "RightElbow"],
            ["RightElbow", "BackU"],
            ["LeftHand",   "LeftElbow"],
            ["LeftElbow",  "BackU"],
            ["RightToe",   "RightHeel"],
            ["RightHeel",  "RightKnee"],
            ["RightKnee",  "BackM"],
            ["LeftToe",    "LeftHeel"],
            ["LeftHeel",   "LeftKnee"],
            ["LeftKnee",   "BackM"],
            ["BackU",      "BackM"],
            ["BackM",      "BackL"],
            ["BackL",      "Tail"]
        ])
        self.jointList=np.array([
            ["RightElbow", "RightHand",  "BackU"],
            ["RightHeel",  "RightToe",   "RightKnee"],
            ["RightKnee",  "RightHeel",  "BackM"],
            ["BackU",      "RightElbow", "LeftElbow"],
            ["BackU",      "RightElbow", "BackM"],
            ["BackU",      "LeftElbow",  "BackM"],
            ["BackM",      "RightKnee",  "LeftKnee"],
            ["BackM",      "RightKnee",  "BackU"],
            ["BackM",      "RightKnee",  "BackL"],
            ["BackM",      "LeftKnee",   "BackU"],
            ["BackM",      "LeftKnee",   "BackL"],
            ["BackM",      "BackU",      "BackL"],
            ["BackL",      "BackM",      "Tail"],
            ["LeftElbow",  "LeftHand",   "BackU"],
            ["LeftHeel",   "LeftToe",    "LeftKnee"],
            ["LeftKnee",   "LeftHeel",   "BackM"]
        ])
    def ind_bp(self,targets):
        is_scalar = isinstance(targets, str)
        targets_array = np.atleast_1d(targets)
        indices = []
        for target in targets_array:
            idx = np.where(self.bp_ind == target)[0]
            if idx.size == 0:
                raise ValueError(f"'{target}' not found in bp_ind.")
            indices.append(int(idx[0]))
        # Return scalar if input was a scalar, otherwise a NumPy array
        return indices[0] if is_scalar else np.array(indices)
    def zscore_normalize_along_1stDim(self,arr):
        mean = np.mean(arr, axis=0, keepdims=True)
        std = np.std(arr, axis=0, keepdims=True)
        std[std == 0] = 1
        return (arr - mean) / std
    def zscore_normalize_along_1stDim_meanstd(self,arr):
        mean = np.mean(arr, axis=0, keepdims=True)
        std = np.std(arr, axis=0, keepdims=True)
        std[std == 0] = 1
        return (arr - mean) / std, mean, std
############################
    def preprocessData(self,cutoff_freq=10,SamplingRate_Hz=10,order=4,triError_threshold_in_mm=4):
        filtered_posture=self.filter_posture(self.catposture,triError_threshold_in_mm=triError_threshold_in_mm)
        self.smoothed_posture=self.smooth_posture(filtered_posture,cutoff_freq,order)
        self.posture_vel=np.diff(self.smoothed_posture,n=1,axis=0)
        self.posture_vel=np.concatenate((self.posture_vel,self.posture_vel[-1:,:]),axis=0)
        self.posture_acc=np.diff(self.smoothed_posture,n=2,axis=0)
        self.posture_acc=np.concatenate((self.posture_acc,self.posture_acc[-1:,:]),axis=0)
        self.joint_angles=self.make_joint_angle(self.smoothed_posture)
        self.joint_vel=np.diff(self.joint_angles,n=1,axis=0)
        self.joint_vel=np.concatenate((self.joint_vel,self.joint_vel[-1:,:]),axis=0)
        self.joint_acc=np.diff(self.joint_angles,n=2,axis=0)
        self.joint_acc=np.concatenate((self.joint_acc,self.joint_acc[-1:,:]),axis=0)
       
        t_trial=self.make_t_downSampled_trial(SamplingRate_Hz)
        cellac_trial=self.sort_cellac_to_trials_AndDownsample(self.cellac,self.f2p_t,t_trial)
        sp_trial=self.sort_sp_to_trials_AndDownsample(self.sp_t,t_trial)
        posture_trial=self.sort_posture_to_trials_AndDownsample(self.smoothed_posture,self.fdlc_t,t_trial)
        posture_vel_trial=self.sort_posture_to_trials_AndDownsample(self.posture_vel,self.fdlc_t,t_trial)
        posture_acc_trial=self.sort_posture_to_trials_AndDownsample(self.posture_acc,self.fdlc_t,t_trial)
        
        cellac_mtypes=self.classify_to_motionTypes(cellac_trial)
        sp_mtypes=self.classify_to_motionTypes(sp_trial)
        posture_mtypes=self.classify_to_motionTypes(posture_trial)
        posture_vel_mtypes=self.classify_to_motionTypes(posture_vel_trial)
        posture_acc_mtypes=self.classify_to_motionTypes(posture_acc_trial)
        
        self.tckd_mtypes=self.combine_p8(cellac_mtypes)
        self.tpkd_mtypes=self.combine_p8(sp_mtypes)
        self.posture_mtypes=self.combine_p8(posture_mtypes)
        self.posture_vel_mtypes=self.combine_p8(posture_vel_mtypes)
        self.posture_acc_mtypes=self.combine_p8(posture_acc_mtypes)
        
        self.ktc_mtypes=self.transpose_allFields(self.tckd_mtypes)
        self.ktp_mtypes=self.transpose_allFields(self.tpkd_mtypes)
        self.kt_posture_mtypes=self.transpose_allFields(self.posture_mtypes)
        self.kt_posture_vel_mtypes=self.transpose_allFields(self.posture_vel_mtypes)
        self.kt_posture_acc_mtypes=self.transpose_allFields(self.posture_acc_mtypes)
    def make_t_downSampled_trial(self,SamplingRate_Hz):
        t_trial=[]
        t_trial=[np.arange(self.t_strial[trial],self.t_etrial[trial],1000/SamplingRate_Hz).round().astype(int) for trial in range(self.t_strial.shape[0])]
        return t_trial
    def sort_cellac_to_trials_AndDownsample(self,cellac,f2p_t,t_trial):
        t_trial_cat=np.concat(t_trial,axis=0)
        cellac_interp=self.interpolate_first_dim(cellac, np.arange(0,cellac.shape[0],1), f2p_t[t_trial_cat])
        t_trial_length=np.array([currtrial.shape[0] for currtrial in t_trial])
        trial_edges=np.cumsum(t_trial_length)[0:-1]
        cellac_trial=np.split(cellac_interp,trial_edges,axis=0)
        return cellac_trial
    def sort_sp_to_trials_AndDownsample(self,sp_t,t_trial):
        t_trial_cat=np.concat(t_trial,axis=0)
        sp_trial_cat=sp_t[t_trial_cat,:]
        t_trial_length=np.array([currtrial.shape[0] for currtrial in t_trial])
        trial_edges=np.cumsum(t_trial_length)[0:-1]
        sp_trial=np.split(sp_trial_cat,trial_edges,axis=0)
        return sp_trial
    def sort_posture_to_trials_AndDownsample(self,posture,fdlc_t,t_trial):
        t_trial_cat=np.concat(t_trial,axis=0)
        posture_trial_cat=posture[fdlc_t[t_trial_cat].round().squeeze().astype(int),:]
        t_trial_length=np.array([currtrial.shape[0] for currtrial in t_trial])
        trial_edges=np.cumsum(t_trial_length)[0:-1]
        posture_trial=np.split(posture_trial_cat,trial_edges,axis=0)
        return posture_trial
    def classify_to_motionTypes(self,data_trial):
        data_mtypes=dict()
        for trial in range(self.motionLog.shape[0]):
            currmotion=self.motionLog[trial]
            if np.isnan(currmotion):
                continue
            currmotion=currmotion.round().astype(int)
            currfield=self.mtypes[currmotion]
            if currfield not in data_mtypes:
                data_mtypes[currfield]=[]
            data_mtypes[currfield].append(data_trial[trial])
        for motion in data_mtypes:
            data_mtypes[motion]=np.stack(data_mtypes[motion],axis=2)
        return data_mtypes
    def combine_p8(self,data_mtypes):
        p8_keys=np.array([f"p8_dir{i}" for i in range(1,9)])
        minsize=np.inf
        for motion in p8_keys:
            currsize=data_mtypes[motion].shape[2]   # size of the third dimension (trial)
            minsize=np.min([currsize,minsize])
        minsize=minsize.round().astype(int)
        p8_data=[]
        for motion in p8_keys:
            p8_data.append(data_mtypes[motion][:,:,:minsize])
        p8_data=np.stack(p8_data,axis=3)    # stack direction to 4th dimension
        
        # remove individual p8_dir and replace them with a single p8 field
        for motion in p8_keys:
            data_mtypes.pop(motion,None)
        data_mtypes['p8']=p8_data
        return data_mtypes
    def transpose_allFields(self,data):
        for key, value in data.items():
            if key != 'p8':
                data[key]=data[key].transpose([2,0,1])
        data['p8']=data['p8'].transpose([3,2,0,1])
        return data
    
    def smooth_posture(self,posture,cutoff_freq,order=4):
        fps=1000/np.mean(np.diff(self.t_fdlc, axis=0))
        smoothed_posture = self.lowpass_filter(posture, cutoff_freq, fps, order=order)
        return smoothed_posture
    def filter_posture(self,posture,triError_threshold_in_mm=4):
        new_posture=[]
        for i in range(self.triError.shape[1]):
            limbs=self.index_posture(posture, i)
            trierror=self.triError[:,i]
            valid_trierror_ind=np.where(trierror<=triError_threshold_in_mm)[0]
            valid_limbs=limbs[valid_trierror_ind,:]
            query_points=np.arange(0,limbs.shape[0],1)
            new_posture.append(self.interpolate_first_dim (valid_limbs, valid_trierror_ind,query_points))
        new_posture=np.concatenate(new_posture,axis=1)
        return new_posture
    def lowpass_filter(self,data, cutoff, fs, order=4):
        from scipy.signal import butter, filtfilt
        """
        Applies a Butterworth low-pass filter along the first dimension of the data.
        Parameters:
            data   : np.ndarray
                    Input data with time as the first dimension.
            cutoff : float
                    The cutoff frequency of the filter in Hz.
            fs     : float
                    The sampling rate in Hz.
            order  : int, optional
                    The order of the filter (default is 4).        
        Returns:
            filtered_data : np.ndarray
                            The filtered data.
        """
        nyquist = 0.5 * fs
        normal_cutoff = cutoff / nyquist
        # Create a Butterworth low-pass filter
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        # Apply zero-phase filtering
        filtered_data = filtfilt(b, a, data, axis=0)
        return filtered_data
    def index_posture(self,array, index):
        if np.mod(array.shape[1], 3) != 0:
            raise ValueError("Array length must be a multiple of 3.")
        newarray=array[:,index*3:index*3+3]
        return newarray
    def make_joint_angle(self,posture):
        """
        Computes joint angles from the catposture data.
        The angles are computed based on the jointList defined in the class.
        """
        joint_angles = []
        for joint in self.jointList:
            joint1 = self.ind_bp(joint[0])
            joint2 = self.ind_bp(joint[1])
            joint3 = self.ind_bp(joint[2])
            vec1 = posture[:, joint2 * 3:joint2 * 3 + 3] - posture[:, joint1 * 3:joint1 * 3 + 3]
            vec2 = posture[:, joint3 * 3:joint3 * 3 + 3] - posture[:, joint2 * 3:joint2 * 3 + 3]
            angle = np.arccos(np.clip(np.sum(vec1 * vec2, axis=1) / (np.linalg.norm(vec1, axis=1) * np.linalg.norm(vec2, axis=1)), -1.0, 1.0))
            joint_angles.append(angle)
        joint_angles = np.stack(joint_angles, axis=1)  # shape (n_timepoints, n_joints)
        return joint_angles
################################# load continuous data
    def preprocessData_continuous(self,cutoff_freq=10,SamplingRate_Hz=10,order=4,triError_threshold_in_mm=4,UseCachedData:bool=True):
        if not UseCachedData:
            filtered_posture=self.filter_posture(self.catposture,triError_threshold_in_mm=triError_threshold_in_mm)
            self.smoothed_posture=self.smooth_posture(filtered_posture,cutoff_freq,order)
            self.posture_vel=np.diff(self.smoothed_posture,n=1,axis=0)
            self.posture_vel=np.concatenate((self.posture_vel,self.posture_vel[-1:,:]),axis=0)
            self.posture_acc=np.diff(self.smoothed_posture,n=2,axis=0)
            self.posture_acc=np.concatenate((self.posture_acc,self.posture_acc[-1:,:]),axis=0)
            self.joint_angles=self.make_joint_angle(self.smoothed_posture)
            self.joint_vel=np.diff(self.joint_angles,n=1,axis=0)
            self.joint_vel=np.concatenate((self.joint_vel,self.joint_vel[-1:,:]),axis=0)
            self.joint_acc=np.diff(self.joint_angles,n=2,axis=0)
            self.joint_acc=np.concatenate((self.joint_acc,self.joint_acc[-1:,:]),axis=0)
        timestamps=self.make_downsampled_t(SamplingRate_Hz, self.f2p_t.shape[0])
        self.cellac_cont=self.downsample_cellac(self.cellac, self.f2p_t, timestamps)
        self.sp_cont=self.downsample_sp(self.sp_t, timestamps)
        self.posture_cont=self.downsample_posture(self.smoothed_posture, self.fdlc_t, timestamps)
        self.posture_vel_cont=self.downsample_posture(self.posture_vel, self.fdlc_t, timestamps)
        self.posture_acc_cont=self.downsample_posture(self.posture_acc, self.fdlc_t, timestamps)
        self.joint_cont=self.downsample_posture(self.joint_angles, self.fdlc_t, timestamps)
        self.joint_vel_cont=self.downsample_posture(self.joint_vel, self.fdlc_t, timestamps)
        self.joint_acc_cont=self.downsample_posture(self.joint_acc, self.fdlc_t, timestamps)
    def make_downsampled_t(self,SamplingRate_Hz,totalTimePoints):
        timestamps=np.arange(0,totalTimePoints,1000/SamplingRate_Hz).round().astype(int)
        return timestamps
    def downsample_cellac(self,cellac,f2p_t,timestamps):
        cellac_interp=self.interpolate_first_dim(cellac, np.arange(0,cellac.shape[0],1), f2p_t[timestamps])
        return cellac_interp
    def downsample_sp(self,sp_t,timestamps):
        sp_t_downsampled=sp_t[timestamps,:]
        return sp_t_downsampled
    def downsample_posture(self,posture,fdlc_t,timestamps):
        posture_downsampled=posture[fdlc_t[timestamps].round().squeeze().astype(int),:]
        return posture_downsampled
################################# set current data
#########################################################
    def sync_and_interpolate_all_data(self,sampling_rate=10):
            timepoints=np.arange(0,len(self.fdlc_t),1000/sampling_rate)
            self.sampled_cellac=self.interpolate_first_dim(self.cellac, self.t_f2p, timepoints)
            self.sampled_catposture=self.interpolate_first_dim(self.catposture, self.t_fdlc, timepoints)
            self.sampled_smoothed_catposture=self.interpolate_first_dim(self.smoothed_catposture, self.t_fdlc, timepoints)
            self.sampled_smoothed_joint_angles=self.interpolate_first_dim(self.joint_angles, self.t_fdlc, timepoints)
            self.sampled_sp_t=self.mod_sp_t[np.round(timepoints).astype(int),:]
            self.sampled_R_t=self.R_t[:,:,np.round(timepoints).astype(int)]
            self.sampled_p_t=self.p_t[np.round(timepoints).astype(int),:]
            self.sampled_euler_angles_t=self.euler_angles_t[np.round(timepoints).astype(int),:]
            ######## the posture data to use for most of the analysis is the smoothed sampled posture below
            self.sampled_posture=self.sampled_smoothed_catposture
    def interpolate_first_dim(self,data, sample_points, query_points, left=None, right=None):
        """
        Interpolates along the first dimension of 'data'.
        
        Parameters:
            data (ndarray): Input array of shape (n, ...).
            sample_points (1D array): Known sample points corresponding to the first axis.
            query_points (1D array): Points where interpolation is computed.
            left (optional): Value for query_points below the minimum sample point.
            right (optional): Value for query_points above the maximum sample point.
            
        Returns:
            ndarray: An array with the first dimension replaced by len(query_points).
        """
        sample_points=np.squeeze(sample_points)
        query_points=np.squeeze(query_points)
        if data.ndim < 2:
            raise ValueError("Data must have at least 2 dimensions.")
        if len(sample_points) != data.shape[0]:
            raise ValueError("Length of sample_points must match the first dimension of data.")
        if len(query_points) == 0:
            raise ValueError("query_points cannot be empty.")
        flat_data=data.reshape(data.shape[0], -1)
        interp_data=[]
        for i in range(flat_data.shape[1]):
            interp_data.append(np.interp(query_points, sample_points, flat_data[:,i], left=left, right=right))
        interp_data=np.stack(interp_data,axis=1)
        interp_data=interp_data.reshape(len(query_points),*data.shape[1:])
        if np.isnan(interp_data).any():
            print("There are NaN values in the result array.")
        # new_shape = list(data.shape)
        # new_shape[0] = len(query_points)
        # result = np.empty(new_shape, dtype=data.dtype)

        # for idx in np.ndindex(data.shape[1:]):
        #     result[(slice(None),) + idx] = np.interp(query_points, sample_points, data[(slice(None),) + idx], left=left, right=right)
        return interp_data
####################
    def plot3d_dktc_p8(self,dktc,dktp,colormap:str='hsv',ax=None,linewidth=2,alpha=0.5,trials_to_combine=1):
        from mpl_toolkits.mplot3d import Axes3D                   # registers 3D projection
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        from scipy.stats import mode
        
        dktc=utf.jumping_mean(dktc, dim=1, num_to_combine=trials_to_combine)

        preTKDc=utf.multiBlockCat(dktc[:,:,:-1,:], dim2squish=[1,0], dim2concat=2).squeeze()
        postTKDc=utf.multiBlockCat(dktc[:,:,1:,:], dim2squish=[1,0], dim2concat=2).squeeze()
        segments=np.stack([preTKDc,postTKDc],axis=1)

        # turn sp parameter to color index
        theta=dktp[:,:,:,0]
        theta=utf.jumping_mean(theta, dim=1, num_to_combine=trials_to_combine)
        m,_=mode(theta,axis=2,keepdims=True)
        theta=np.broadcast_to(m,theta.shape)       # replace theta values with its mode
        theta=utf.multiBlockCat(theta[:,:,:-1], dim2squish=[1,0], dim2concat=2).squeeze()
        cmap=self.make_circular_colormap(theta,colormap_name=colormap)

        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(projection='3d')

        # overlay the loop edges
        lc3d = Line3DCollection(segments,
                                linewidths=linewidth,
                                alpha=alpha,
                                colors=cmap)
        ax.add_collection3d(lc3d,autolim=False)
        
        # hook the depth‐sorter, binding in our data
        callback = partial(
            self.depth_sort,
            ax=ax,
            segments=segments,
            colors=cmap,
            collection=lc3d)
        ax.figure.canvas.mpl_connect('draw_event', callback)

        # autoscale so everything fits
        ax.auto_scale_xyz(dktc[:,:,:,0],dktc[:,:,:,1],dktc[:,:,:,2],had_data=False)
        self.equalize_data_scale(ax)
        # ax.set_title(f"{CycleIND+1}th Longest lived loop. {csp.currPlaneName}, p8\nPCA to {pcadims2use} dims, cPCA to {cPCAdims} dims using p8, weight={weight}")
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        # plt.show()
    def plot2d_dktc_p8(self, dktc, dktp, colormap='hsv', ax=None, linewidth=2, alpha=0.5, trials_to_combine=1):
        from matplotlib.collections import LineCollection
        from scipy.stats import mode
        dktc = utf.jumping_mean(dktc, dim=1, num_to_combine=trials_to_combine)
        preTKDc  = utf.multiBlockCat(dktc[:,:,:-1,:], dim2squish=[1,0], dim2concat=2).squeeze()
        postTKDc = utf.multiBlockCat(dktc[:,:,1: ,:], dim2squish=[1,0], dim2concat=2).squeeze()
        segments = np.stack([preTKDc, postTKDc], axis=1)
        theta = utf.jumping_mean(dktp[:,:,:,0], dim=1, num_to_combine=trials_to_combine)
        m, _ = mode(theta, axis=2, keepdims=True)
        theta = np.broadcast_to(m, theta.shape)
        theta = utf.multiBlockCat(theta[:,:,:-1], dim2squish=[1,0], dim2concat=2).squeeze()
        cmap = self.make_circular_colormap(theta, colormap_name=colormap)
        if ax is None: fig, ax = plt.subplots()
        lc2d = LineCollection(segments, linewidths=linewidth, alpha=alpha, colors=cmap)
        ax.add_collection(lc2d)
        ax.autoscale_view()
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
        # plt.show()
    def scatter3d_dktc_p8(self,dktc,dktp,colormap:str='hsv',dotsize=2,alpha=0.5,ax=None,trials_to_combine=1):
        from mpl_toolkits.mplot3d import Axes3D                   # registers 3D projection
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        from scipy.stats import mode

        dktc=utf.jumping_mean(dktc, dim=1, num_to_combine=trials_to_combine)
        TKDc=utf.multiBlockCat(dktc, dim2squish=[1,0], dim2concat=2).squeeze()
        
        # turn sp parameter to color index
        theta=dktp[:,:,:,0]
        theta=utf.jumping_mean(theta, dim=1, num_to_combine=trials_to_combine)
        m,_=mode(theta,axis=2,keepdims=True)
        theta=np.broadcast_to(m,theta.shape)       # replace theta values with its mode
        theta=utf.multiBlockCat(theta, dim2squish=[1,0], dim2concat=2).squeeze()
        cmap=self.make_circular_colormap(theta,colormap_name=colormap)

        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(projection='3d')
        ax.scatter(TKDc[:,0],TKDc[:,1],TKDc[:,2],c=cmap,s=dotsize,alpha=alpha)

        # autoscale so everything fits
        ax.auto_scale_xyz(dktc[:,:,:,0],dktc[:,:,:,1],dktc[:,:,:,2])
        self.equalize_data_scale(ax)
        # ax.set_title(f"{CycleIND+1}th Longest lived loop. {csp.currPlaneName}, p8\nPCA to {pcadims2use} dims, cPCA to {cPCAdims} dims using p8, weight={weight}")
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        # plt.show()
    def scatter3d_ktc_rot(self,ktc,ktp,colormap:str='hsv',dotsize=2,alpha=0.5,ax=None,trials_to_combine=1):
        from mpl_toolkits.mplot3d import Axes3D                   # registers 3D projection
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        
        ktc=utf.jumping_mean(ktc, dim=0, num_to_combine=trials_to_combine)
        TKc=utf.multiBlockCat(ktc, dim2squish=0, dim2concat=1).squeeze()
        kt_theta=ktp[:,:,0]
        kt_color=self.make_circular_colormap(kt_theta,colormap)
        kt_color=utf.jumping_mean(kt_color,dim=0,num_to_combine=trials_to_combine)
        TK_color=utf.multiBlockCat(kt_color,dim2squish=0,dim2concat=1).squeeze()

        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(projection='3d')
        ax.scatter(TKc[:,0],TKc[:,1],TKc[:,2],c=TK_color,s=dotsize,alpha=alpha)
        
        # autoscale so everything fits
        ax.auto_scale_xyz(ktc[:,:,0],ktc[:,:,1],ktc[:,:,2])
        self.equalize_data_scale(ax)
        # ax.set_title(f"{CycleIND+1}th Longest lived loop. {csp.currPlaneName}, p8\nPCA to {pcadims2use} dims, cPCA to {cPCAdims} dims using p8, weight={weight}")
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        # plt.show()
    def plot3d_ktc_rot(self,ktc,ktp,colormap:str='hsv',linewidth=2,alpha=0.5,ax=None,trials_to_combine=1):
        from mpl_toolkits.mplot3d import Axes3D                   # registers 3D projection
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        
        ktc=utf.jumping_mean(ktc, dim=0, num_to_combine=trials_to_combine)        
        preTKc=utf.multiBlockCat(ktc[:,:-1,:], dim2squish=0, dim2concat=1).squeeze()
        postTKc=utf.multiBlockCat(ktc[:,1:,:], dim2squish=0, dim2concat=1).squeeze()
        segments=np.stack([preTKc,postTKc],axis=1)
                
        kt_theta=ktp[:,1:,0]
        kt_color=self.make_circular_colormap(kt_theta,colormap)
        kt_color=utf.jumping_mean(kt_color,dim=0,num_to_combine=trials_to_combine)
        TK_color=utf.multiBlockCat(kt_color,dim2squish=0,dim2concat=1).squeeze()

        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(projection='3d')

        # overlay the loop edges
        lc3d = Line3DCollection(segments[...,:3],
                                 linewidths=linewidth,
                                 alpha=alpha,
                                 colors=TK_color)
        ax.add_collection3d(lc3d,autolim=False)

        # autoscale so everything fits
        ax.auto_scale_xyz(ktc[:,:,0],ktc[:,:,1],ktc[:,:,2],had_data=False)
        self.equalize_data_scale(ax)
        # hook the depth‐sorter, binding in our data
        callback = partial(
            self.depth_sort,
            ax=ax,
            segments=segments,
            colors=TK_color,
            collection=lc3d)
        ax.figure.canvas.mpl_connect('draw_event', callback)

        # ax.set_title(f"{CycleIND+1}th Longest lived loop. {csp.currPlaneName}, p8\nPCA to {pcadims2use} dims, cPCA to {cPCAdims} dims using p8, weight={weight}")
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        # plt.show()
    def plot2d_ktc_rot(self, ktc, ktp, colormap='hsv', linewidth=2, alpha=0.5, ax=None, trials_to_combine=1):
        from scipy.stats import mode
        from matplotlib.collections import LineCollection
        ktc = utf.jumping_mean(ktc, dim=0, num_to_combine=trials_to_combine)
        preTKc  = utf.multiBlockCat(ktc[:,:-1,:], dim2squish=0, dim2concat=1).squeeze()
        postTKc = utf.multiBlockCat(ktc[:,1: ,:], dim2squish=0, dim2concat=1).squeeze()
        segments = np.stack([preTKc, postTKc], axis=1)
        kt_theta = ktp[:,1:,0]
        kt_color = self.make_circular_colormap(kt_theta, colormap_name=colormap)
        kt_color = utf.jumping_mean(kt_color, dim=0, num_to_combine=trials_to_combine)
        TK_color = utf.multiBlockCat(kt_color, dim2squish=0, dim2concat=1).squeeze()
        if ax is None: fig, ax = plt.subplots()
        lc2d = LineCollection(segments, linewidths=linewidth, alpha=alpha, colors=TK_color)
        ax.add_collection(lc2d); 
        ax.autoscale_view()
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('PC1'); ax.set_ylabel('PC2'); plt.show()
    def scatter3d_ktc(self,ktc,kt_param,param_min_max=None,colormap:str='hsv',dotsize=2,alpha=0.5,ax=None,trials_to_combine=1,code_direction:bool=False):
        from mpl_toolkits.mplot3d import Axes3D                   # registers 3D projection
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        
        ktc=utf.jumping_mean(ktc, dim=0, num_to_combine=trials_to_combine)
        TKc=utf.multiBlockCat(ktc, dim2squish=0, dim2concat=1).squeeze()
        
        kt_param=utf.jumping_mean(kt_param, dim=0, num_to_combine=trials_to_combine)
        if code_direction:
            kt_color=self.make_colormap(kt_param,colormap_name=colormap,min_max=param_min_max,dim_of_direction_to_code=1)
        else:
            kt_color=self.make_colormap(kt_param,colormap_name=colormap,min_max=param_min_max)
        TK_color=utf.multiBlockCat(kt_color, dim2squish=0, dim2concat=1).squeeze()

        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(projection='3d')
        ax.scatter(TKc[:,0],TKc[:,1],TKc[:,2],c=TK_color,s=dotsize,alpha=alpha)
        
        # autoscale so everything fits
        ax.auto_scale_xyz(ktc[:,:,0],ktc[:,:,1],ktc[:,:,2])
        self.equalize_data_scale(ax)
        # ax.set_title(f"{CycleIND+1}th Longest lived loop. {csp.currPlaneName}, p8\nPCA to {pcadims2use} dims, cPCA to {cPCAdims} dims using p8, weight={weight}")
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        # plt.show()
    def plot3d_ktc(self,ktc,kt_param,param_min_max=None,colormap:str='hsv',linewidth=2,alpha=0.5,ax=None,trials_to_combine=1,code_direction:bool=False):
        from mpl_toolkits.mplot3d import Axes3D                   # registers 3D projection
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        
        ktc=utf.jumping_mean(ktc, dim=0, num_to_combine=trials_to_combine)        
        preTKc=utf.multiBlockCat(ktc[:,:-1,:], dim2squish=0, dim2concat=1).squeeze()
        postTKc=utf.multiBlockCat(ktc[:,1:,:], dim2squish=0, dim2concat=1).squeeze()
        segments=np.stack([preTKc,postTKc],axis=1)

        kt_param=utf.jumping_mean(kt_param[:,1:], dim=0, num_to_combine=trials_to_combine)
        if code_direction:
            kt_color=self.make_colormap(kt_param,colormap_name=colormap,min_max=param_min_max,dim_of_direction_to_code=1)
        else:
            kt_color=self.make_colormap(kt_param,colormap_name=colormap,min_max=param_min_max)
        TK_color=utf.multiBlockCat(kt_color, dim2squish=0, dim2concat=1).squeeze()

        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(projection='3d')
        lc3d = Line3DCollection(segments[...,:3],
                                 linewidths=linewidth,
                                 alpha=alpha,
                                 colors=TK_color)
        ax.add_collection3d(lc3d,autolim=False)
        
        # autoscale so everything fits
        ax.auto_scale_xyz(ktc[:,:,0],ktc[:,:,1],ktc[:,:,2],had_data=False)
        self.equalize_data_scale(ax)
        # hook the depth‐sorter, binding in our data
        callback = partial(
            self.depth_sort,
            ax=ax,
            segments=segments,
            colors=TK_color,
            collection=lc3d)
        ax.figure.canvas.mpl_connect('draw_event', callback)
        
        # ax.set_title(f"{CycleIND+1}th Longest lived loop. {csp.currPlaneName}, p8\nPCA to {pcadims2use} dims, cPCA to {cPCAdims} dims using p8, weight={weight}")
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        # plt.show()
    def plot2d_ktc(self, ktc, kt_param, param_min_max=None, colormap='hsv',
               linewidth=2, alpha=0.5, ax=None, trials_to_combine=1, code_direction=False):
        ktc = utf.jumping_mean(ktc, dim=0, num_to_combine=trials_to_combine)
        preTKc = utf.multiBlockCat(ktc[:,:-1,:], dim2squish=0, dim2concat=1).squeeze(); \
        postTKc = utf.multiBlockCat(ktc[:,1: ,:], dim2squish=0, dim2concat=1).squeeze()
        segments = np.stack([preTKc, postTKc], axis=1)
        ktp = utf.jumping_mean(kt_param[:,1:], dim=0, num_to_combine=trials_to_combine)
        color = (self.make_colormap(ktp, colormap_name=colormap, min_max=param_min_max, dim_of_direction_to_code=1)
                 if code_direction else
                 self.make_colormap(ktp, colormap_name=colormap, min_max=param_min_max))
        TK_color = utf.multiBlockCat(color, dim2squish=0, dim2concat=1).squeeze()
        if ax is None: fig, ax = plt.subplots()
        lc2d = LineCollection(segments, linewidths=linewidth, alpha=alpha, colors=TK_color)
        ax.add_collection(lc2d)
        ax.autoscale_view()
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('PC1'); ax.set_ylabel('PC2'); plt.show()
    def make_circular_colormap(self,theta,colormap_name:str='hsv'):
        # returns an array with the same dimension as theta but the last dimension correponds to color (size 3)
        int_theta=np.round(255*np.mod(theta,2*np.pi)/(2*np.pi)).astype(int)
        cmap=utf.getColorMap(colormap_name, size=256)
        theta_color=cmap[int_theta]
        return theta_color
    def make_colormap(self,param,colormap_name:str='hsv',min_max=None,dim_of_direction_to_code=None):
        param=np.asarray(param)
        if min_max is None:
            min_max=[param.min(),param.max()]
        min_max=np.asarray(min_max)
        param=np.clip(param,min_max[0],min_max[1])
        param=(param-min_max[0])/(min_max[1]-min_max[0])
        min_max=np.array([0,1])
        if dim_of_direction_to_code is not None:
            dparam=np.diff(param,n=1,axis=dim_of_direction_to_code)
            dparam=np.concatenate([dparam,np.take(dparam,[-1],axis=dim_of_direction_to_code)],axis=dim_of_direction_to_code)
            neg_log=(dparam<0)
            param[neg_log]=-param[neg_log]
            min_max[0]=-min_max[1]
            param=(param-min_max[0])/(min_max[1]-min_max[0])
        int_param=np.round(255*param).astype(int)
        cmap=utf.getColorMap(colormap_name,size=256)
        param_color=cmap[int_param]
        return param_color
    def depth_sort(self, event, ax, segments, colors, collection):
        """Draw‐event callback to depth‐sort a 3D Line3DCollection."""
        if event.canvas.figure is not ax.figure:
            return
    
        from mpl_toolkits.mplot3d import proj3d
        import numpy as np
    
        # compute each segment’s midpoint
        mids = segments.mean(axis=1)   # (n_segs, 3)
        # get the projection matrix
        proj = ax.get_proj()
        # project to display coords and grab Z
        _, _, zs = proj3d.proj_transform(
            mids[:, 0], mids[:, 1], mids[:, 2], proj
        )
        # sort so that farthest → first, nearest → last
        order = np.argsort(zs)[::-1]  
    
        collection.set_segments(segments[order])
        collection.set_color   (colors[order])
    def equalize_data_scale(self,ax):
        """
        After auto-scaling, expand all 3D axis limits so that
        X, Y and Z each span the same numeric range.
        """
        # grab current limits
        x0, x1 = ax.get_xlim3d()
        y0, y1 = ax.get_ylim3d()
        z0, z1 = ax.get_zlim3d()
    
        # compute centers and half-spans
        xc, yc, zc = np.mean([x0, x1]), np.mean([y0, y1]), np.mean([z0, z1])
        hx, hy, hz = (x1 - x0)/2, (y1 - y0)/2, (z1 - z0)/2
    
        # pick the largest half-span
        H = max(hx, hy, hz)
    
        # set each axis to center ± H
        ax.set_xlim(xc - H, xc + H)
        ax.set_ylim(yc - H, yc + H)
        ax.set_zlim(zc - H, zc + H)
####################
    def zscore_posture(self):
        self.zscore_catposture,self.zscore_catposture_mean,self.zscore_catposture_std=self.zscore_normalize_along_1stDim_meanstd(self.catposture)
        self.zscore_smoothed_catposture,self.zscore_smoothed_catposture_mean,self.zscore_smoothed_catposture_std=self.zscore_normalize_along_1stDim_meanstd(self.smoothed_catposture)
        self.zscore_sampled_catposture,self.zscore_sampled_catposture_mean,self.zscore_sampled_catposture_std=self.zscore_normalize_along_1stDim_meanstd(self.sampled_catposture)
        self.zscore_sampled_smoothed_catposture,self.zscore_sampled_smoothed_catposture_mean,self.zscore_sampled_smoothed_catposture_std=self.zscore_normalize_along_1stDim_meanstd(self.sampled_smoothed_catposture)
        self.sampled_posture_mean=self.zscore_sampled_smoothed_catposture_mean
        self.sampled_posture_std=self.zscore_sampled_smoothed_catposture_std
    def make_eigen_posture(self):
        self.zscore_posture()
        pca=PCA(data=self.sampled_posture)
        self.pca=pca
        self.eigen_posture=pca.scores_
        self.eigen_posture_varRatio=pca.explained_variance_ratio_
        return self.eigen_posture
    def make_pos_from_eigen_posture(self,eigen_posture):
        self.posture_eigen_reconstruct=self.pca.inverse_transform(eigen_posture)
        return self.posture_eigen_reconstruct
###################### downsample and classify all data to trials
    def downsample_and_classify_all_data(self, SamplingRate_Hz=10):
        self.SamplingRate_Hz=SamplingRate_Hz
        downsampled_t=self.make_downsampled_t_4eachTrials(SamplingRate_Hz)
        self.cellac_trial=self.sort_cellac_to_trials_and_downsample(self.cellac, downsampled_t)
    def sort_cellac_to_trials_and_downsample(self,cellac, downsampled_t):
        vec_t = np.concatenate(downsampled_t)
        # assume org.f2p_t is an array-like mapping; if it's callable, use org.f2p_t(vec_t) instead
        f2p_ind = self.f2p_t[vec_t]
        # Interpolate each column of cellac onto f2p_ind (MATLAB’s 1-based indexing)
        xp = np.arange(1, cellac.shape[0] + 1)
        cols = cellac.shape[1]
        interpcellac=[np.interp(f2p_ind, xp, cellac[:, j]) for j in range(cols)]
        interpcellac = np.stack(interpcellac, axis=1)
        trial_length=[len(arr) for arr in downsampled_t]
        split_indices = np.cumsum(trial_length)[:-1]
        cellac_trials = np.split(interpcellac, split_indices, axis=0)
        return cellac_trials
    def classify_trials(self,data_tp,downsampled_t):
        # data_tp: 2D array of shape (n_timepoints, n_features)
        data_tp_trials = []
        
    def make_downsampled_t_4eachTrials(self, sampling_rate_hz): # 0-based indexing!
        downsampled_t = []
        step = 1000.0 / sampling_rate_hz
        for start, end in zip(self.t_strial, self.t_etrial):
            t = np.arange(start, end + 1, dtype=float)
            N = t.size
            xq = np.arange(0.0, N + 0.0, step)  # query points
            xp = np.arange(0, N , dtype=float)
            td = np.interp(xq, xp, t)
            downsampled_t.append(np.round(td).astype(int))
        return downsampled_t
    

######################
#%%
if __name__ == "__main__":
    csp=load_cell_SP_posture()
    csp.load_directory_list()
    csp.load_data(0)
    csp.filter_posture(cutoff_freq=10,triError_threshold_in_mm=4)
    #%%
    csp.downsample_and_classify_all_data(SamplingRate_Hz=10)
    # csp.sync_and_interpolate_all_data(sampling_rate=10)
    #%%
    csp.zscore_posture()