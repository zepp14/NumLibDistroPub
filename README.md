# SDMI Function Library Documention

Simple overview of use/purpose.

## Description

An in-depth paragraph about your project and overview of use.

## Installation 
Install python dependancies:
```
python -m pip install numpy scipy sklearn matplotlib
```
Clone and Unzip Repo:
```
git clone https://github.com/zepp14/NumLibDistroPub.git
cd NumLibDistroPub
unzip SDMI_Base_Folder_V2.zip
cd SDMI_Base_Folder
```
Build and Install:

```
python setup.py build
sudo python setup.py install
```

Test:
```
import sdmiFuncLib as sfl
sfl.generalParaEq(0.2, [1, 3], [0,0], 4, 0, 0)
``` 

Expected Output: (1.909796863652831, 3.1844249545367442)
## Examples

## API Reference

### Dominant Region Parametric Generation Functions
> #### generalParaEq(theta, Invader, Defender, alpha, r, T)
>> ***Description:***  This function generates a single X,Y coordinate point on the boundry of the dominant region based on the parametric angular coordinate ___theta___.
>> 
>> ***Parameters:***
>> - theta: scalar, parametric angular coordinate of point on dominant region boundry
>> - Invader: [X,Y], position of invader 
>> - Defender: [X,Y], position of defender
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - r: scalar, capture range of defender 
>> - T: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> 
>> ***Output:***
>> - [X,Y]: coordinate point on the boundry of the dominant region based on the parametric angular coordinate


> #### paraEqPath(Invader, Defender, alpha, r, T, [bounds], [num]):
>> ***Description:***  This function generates an array of (#)num points on dominant region boundry, evenly distributed between bounds placed on the parametric angular coordinate.
>> 
>> ***Parameters:***
>> - Invader: [X,Y], position of invader 
>> - Defender: [X,Y], position of defender
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - r: scalar, capture range of defender 
>> - T: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - bounds [optional]: [LowerBound, Upperbound ], default = [-pi, pi]
>> - num [optional]: number of points generated, default = 314
>> 
>> ***Output:***
>> - [[X, Y],.. [Xn, Yn]]: Array Size [num,2]

### Convex Target Region Function (Assumed Circular)
> #### f(Pos, Invader, Defender, alpha, r, T, G_circ)
>> ***Description:***  Convex Target Region function f(p) >= 0 is out of region and f(p) < 0 inside region.
>> 
>> ***Parameters:***
>> - Pos: [X,Y], Any 2-D point on game field
>> - Invader: [X,Y], position of invader 
>> - Defender: [X,Y], position of defender
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - r: scalar, capture range of defender 
>> - T: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> 
>> ***Output:***
>> - g: Value of convex target region evaluated at point Pos

### Simulated Trial Plotting Tools
> #### createInvDict(I_array)
>> ***Description:***  Takes in 2D array of Invader location, [[X, Y],.. [Xn, Yn]], creates python dictionary object with invader data.
>> 
>> ***Parameters:***
>> - I_array: [[X, Y],.. [Xn, Yn]], positions of invader
>> 
>> ***Output:***
>> - I_dict: Dictionary Object with invader location
>> - I_list: List of Invader names

> #### plotTrial(Defender, I_dict, Capt_Order, G_circ, alpha, r, axisSz=[-10,10,-5,10], plotFlag=1, T = 0, effMethod = 0, plt=plt)
>> ***Description:*** This function plays out an SDMI game from initial condition with a predefined capture order defined by the order of Capt_Order
>> It can plot out the trial, and returns some key stats about the trial
>> 
>> ***Parameters:***
>> - Defender: [X,Y], position of defender
>> - I_dict: Dictionary Object with invader location
>> - Capt_Order: ex. ['I1', 'I3', 'I2'] Ordered List of Invaders (Capture Order)
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - r: scalar, capture range of defender 
>> - axisSz [optional]: size of gamefield plotted
>> - plotFlag [optional]: set to 1(default) to plot to 0 to suppress
>> - T: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - effMethod: scalar, effeciency calculation method, 2 = Distance Traveled, 1 = Distance from Target region, 0 = hybrid effeciency method [Han Fu, 2020]
>> 
>> ***Output:***
>> - TotEff: Total efficency of Capture order plan
>> - P_star: 2D array of optimal capture points
>> - T: Time duration of game 

> #### loadTrialData(filename, trialNo)
>> ***Description:*** Load trial from dataset (.npz) file and retrieve key info on trial results 
>> 
>> ***Parameters:***
>> - filename: str, name of npz file
>> - trialNo: scalar, Number ID of trial in question
>> 
>> ***Output:***
>> - BestOrder: Total efficency of Capture order plan
>> - BestScore: 2D array of optimal capture points
>> - Iarr: Time duration of game
>> - NumCount: Number of SDSI sub-problems solved
>> - meanVect: Mean of the aggregated performance of dataset
>> - StdData: Standard deviation of aggregated perfromance of the dataset 

 ### Trial Scoring Tools
> #### flightPathScore(P_star, Defender)
>> ***Description:***  Computes a 'fliability' score of the flight path based on how sharp turns in the flight path are.
>> 
>> ***Parameters:***
>> - P_star: 2D array of optimal capture points
>> - Defender: [X,Y], position of defender
>> 
>> ***Output:***
>> - Score: score of capture order plan

> #### numInTarget(P_star, G_circ, Capt_Order, I_dict)
>> ***Description:***  Calculates how many invaders have made it past the defences
>> 
>> ***Parameters:***
>> - P_star: 2D array of optimal capture points
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - Capt_Order: ex. ['I1', 'I3', 'I2'] Ordered List of Invaders (Capture Order)
>> - I_dict: Dictionary Object with invader location
>> 
>> ***Output:***
>> - Score: Number of invaders that have reached the target region

### Capture Order Decision Tools

> #### computeBestOrderEnum(I_array, r, alpha, Defender, G_circ, dist=0, eMethod = 0 ):
>> ***Description:***  Computes the best capture order through brute force enumeration
>> 
>> ***Parameters:***
>> - I_array: [[X, Y],.. [Xn, Yn]], 2D array of Invader Locations
>> - r: scalar, capture range of defender 
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - Defender: [X,Y], position of defender
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - dist: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - eMethod: scalar, effeciency calculation method, 2 = Distance Traveled, 1 = Distance from Target region, 0 = hybrid effeciency method [Han Fu, 2020]
>> 
>> ***Output:***
>> - Best Order: Str Array, Best Capture Order
>> - P_star: 2D array of optimal capture points
>> - meanVect: Mean of the aggregated performance of dataset
>> - StdData: Standard deviation of aggregated perfromance of the dataset 
>> - T: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
 
> #### plotEnumerations(I_array, r, alpha, Defender, G_circ, dist=0, eMethod = 0, plt = plt)
>> ***Description:***  This function computes stats on the total set of feasible solution, used for research purposes
>> 
>> ***Parameters:***
>> - I_array: [[X, Y],.. [Xn, Yn]], 2D array of Invader Locations
>> - r: scalar, capture range of defender 
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - Defender: [X,Y], position of defender
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - dist: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - eMethod: scalar, effeciency calculation method, 2 = Distance Traveled, 1 = Distance from Target region, 0 = hybrid effeciency method [Han Fu, 2020]
>> - plt: plot figure ID
>> 
>> ***Output:***
>> - Best Order: Str Array, Best Capture Order
>> - P_star: 2D array of optimal capture points
>> - meanVect: Mean of the aggregated performance of dataset
>> - StdData: Standard deviation of aggregated perfromance of the dataset 

> #### chooseCluster(r, alpha, Defender, G_circ, C1, C2, dist=0, effMethod=0)
>> ***Description:***  Facilitates The virtual invader approximation and picks the best cluster to visit first
>> 
>> ***Parameters:***
>> - r: scalar, capture range of defender 
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - Defender: [X,Y], position of defender
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - C1: [X,Y], Centroid of invader cluster 1
>> - C2: [X,Y], Centroid of invader cluster 2
>> - dist: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - effMethod: scalar, effeciency calculation method, 2 = Distance Traveled, 1 = Distance from Target region, 0 = hybrid effeciency method [Han Fu, 2020]
>> 
>> ***Output:***
>> - BestIdx: Best Option between C1 and C2


> #### orderByCluster(I_dict, I_list, r, alph, Defender, G_circ, dist=0, eMethod=0)
>> ***Description:***  Facilitates KMeans Clustering of the Invaders and outputs Invader Clusters
>> 
>> ***Parameters:***
>> - I_dict: Dictionary Object with invader location
>> - I_list: List of Invader names
>> - r: scalar, capture range of defender 
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - Defender: [X,Y], position of defender
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - dist: scalar, Time offset before defender takes optimal pursuit, used for SDMI cases
>> - eMethod: scalar, effeciency calculation method, 2 = Distance Traveled, 1 = Distance from Target region, 0 = hybrid effeciency method [Han Fu, 2020]
>> 
>> ***Output:***
>> - invInClust0: Str Array, Cluster of Invaders1
>> - invInClust1: Str Array Cluster of Invaders2
>> - ClusterLabel: Labels of Clusters

> #### adaptiveClusteringOrder (I_array, Defender, G_circ, alpha, r, typeEff=0, InitMode=True, reInit = False) 
>> ***Description:***  This function implements the adaptive clustering approach to find a good capture.
>> 
>> ***Parameters:***
>> - I_array: [[X, Y],.. [Xn, Yn]], 2D array of Invader Locations
>> - Defender: [X,Y], position of defender
>> - G_circ: [X_c, Y_c, R], Circular Function Parameters, [X_c, Y_c] = Center of target region, R = Radius of dominant region
>> - alpha: scalar, ratio of speeds, V_defender / V_invader
>> - r: scalar, capture range of defender 
>> - typeEff: scalar, effeciency calculation method, 2 = Distance Traveled, 1 = Distance from Target region, 0 = hybrid effeciency method [Han Fu, 2020]
>> - InitMode [optional]: Initialization Flag (Keep True for standard use)
>> - reInit [optional]: Allow Reinitialization Flag (Default=False: If true, it can be less effecient but the results are more reliable)


>> 
>> ***Output:***
>> - predictedOrder: Str Array, Approximated Optimal Capture Order


