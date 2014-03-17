# get the ncbi id_names from the summary file
cat summary.txt  | awk '{print $1}' > u

# find all the available in summary.txt
for s in `cat u`; do find . -name *${s%.[0-9]}* >> files; done &
