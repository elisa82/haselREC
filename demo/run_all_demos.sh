#run all demos
echo 'demo 1'
python3 -m haselrec demo/job_selection_1.ini --run-selection
echo 'demo 2'
python3 -m haselrec demo/job_selection_1.ini --check-NGArec
echo 'demo 3'
python3 -m haselrec demo/job_selection_2.ini --run-selection
echo 'demo 4'
python3 -m haselrec demo/job_selection_2.ini --run-scaling
echo 'demo 5'
python3 -m haselrec demo/job_selection_3.ini --run-complete
echo 'demo 6'
python3 -m haselrec demo/job_selection_4.ini --run-selection
echo 'demo 7'
python3 -m haselrec demo/job_selection_5.ini --run-selection
echo 'demo 8'
python3 -m haselrec demo/job_selection_6.ini --run-selection
