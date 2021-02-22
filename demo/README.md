# HaselREC 
HAzard-based SELection of RECords

# Demo

Some demos are currently available to run the different execution modes.

## Demo 1 (execution modes 1 + 4)

* To test the execution mode 1 - record selection only

Run:

```
cd ..
./haselrec.py demo/job_selection_1.ini --run-selection
```

The output files are stored in *demo/Output_1* 

Four subfolders are created:
- PGA-site_1-poe-1   
- PGA-site_1-poe-2  
- SA(0.2)-site_1-poe-1  
- SA(0.2)-site_1-poe-2

The folder name has the structure: *IM-site_{#num_site}-poe-{#num_poe}*
*IM* is the required intensity measure
*#num_site* is the site number
*#num_poe* is the probability of exceedance number

Each folder cointains 5 files (3 plots and 2 txt files). 
The 2 txt files have the extension:
- *summary_selection.txt* - Contains the summary of record selection
- *_CS.txt* - Contains the CS

* To test the execution mode 4 - identification of missing NGA-West2 records

Run:
```
cd ..
./haselrec.py demo/job_selection_1.ini --check-NGArec
```

The output files are stored in *demo/Output_1* 

A file called *missing_NGArec.txt* is produced
It contains the ID of all missing NGArecords

## Demo 2 (execution modes 1 + 2)

* To test the execution mode 1 - record selection only

Run:
```
cd ..
./haselrec.py demo/job_selection_2.ini --run-selection
```
The output files are stored in *demo/Output_2*

One subfolder is created:
- PGA-site_1-poe-1

The folder name has the structure: *IM-site_{#num_site}-poe-{#num_poe}*
*IM* is the required intensity measure
*#num_site* is the site number
*#num_poe* is the probability of exceedance number

The folder cointains 5 files (3 plots and 2 txt files). 
The 2 txt files have the extension:
- *summary_selection.txt* - Contains the summary of record selection
- *_CS.txt* - Contains the CS

* To test the execution mode 2 - record scaling only
An internet connection is required, along with the *token.txt* file 
to download ESM recordings. 

Run:

```
cd ..
./haselrec.py demo/job_selection_2.ini --run-scaling
```
10 files (*nGMx2*) called *GMR_time_scaled_acc_#GMnum_#comp.txt* are produced
*GMnum* is a sequential number ranging from 1 to *nGM*
*#comp* can be *1* or *2* to indicate the component of motion 
They contains the selected scaled accelerograms

## Demo 3

Tests the execution mode 3 - record selection and scaling.
An internet connection is required, along with the *token.txt* file 
to download ESM recordings.

Run:

```
cd ..
./haselrec.py demo/job_selection_3.ini --run-complete
```

The output files are stored in *demo/Output_3*.
They are the same as *demo/Output_2*



