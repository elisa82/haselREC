# HaselREC - Demos

Some demos are currently available to run the different execution modes.

## Demo 1 (execution modes 1 + 4)

* To test the execution mode 1 - record selection only

Run:

```
python -m haselrec demo/job_selection_1.ini --run-selection
```

The output files are stored in `demo/Output_1` 

* To test the execution mode 4 - identification of missing NGA-West2 records

Run:
```
python -m haselrec demo/job_selection_1.ini --check-NGArec
```

The output files are stored in `demo/Output_1` 

## Demo 2 (execution modes 1 + 2)

* To test the execution mode 1 - record selection only

Run:
```
python -m haselrec demo/job_selection_2.ini --run-selection
```
The output files are stored in `demo/Output_2`

* To test the execution mode 2 - record scaling only
An internet connection is required, along with the `token.txt` file 
to download ESM recordings. 

Run:
```
python -m haselrec demo/job_selection_2.ini --run-scaling
```
The output files are stored in `demo/Output_2`

## Demo 3

Tests the execution mode 3 - record selection and scaling.
An internet connection is required, along with the `token.txt` file 
to download ESM recordings.

Run:
```
python -m haselrec demo/job_selection_3.ini --run-complete
```
The output files are stored in `demo/Output_3`.
