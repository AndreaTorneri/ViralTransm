import csv

rows = csv.DictReader(open("antiviral.csv"))

qsub = "qsub -A ${budget} -batch antiviral.pbs "

for row in rows:
    vars = []
    for key in row:
        vars.append(key + "=" + row[key])
    print(qsub + "-v \"" + ",".join(vars) + "\"")
