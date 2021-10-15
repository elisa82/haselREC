#run all demos
echo 'demo 1'
python -m haselrec demo/job_selection_1.ini --run-selection
echo 'demo 2'
python -m haselrec demo/job_selection_1.ini --check-NGArec
echo 'demo 3'
python -m haselrec demo/job_selection_2.ini --run-selection
echo 'demo 4'
python -m haselrec demo/job_selection_2.ini --run-scaling
echo 'demo 5'
python -m haselrec demo/job_selection_3.ini --run-complete
