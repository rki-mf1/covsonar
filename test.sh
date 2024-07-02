#!/bin/bash

# set bash mode to strict
set -euo pipefail
IFS=$'\n\t'

CWD=$(pwd)
DIR=$(cd "$(dirname "$0")"; pwd)

#conda env
conda env create --name covsonar_ci --file $DIR/sonar.env.yml
source $(conda info --root)/etc/profile.d/conda.sh
set +u
conda activate covsonar_ci
set -u

#doctest
cd $DIR/lib
./sonardb.py
cd "$CWD"

#end-to-end test
db=$(mktemp /tmp/sonar.XXXXXX)
exec 3>"$db"
rm "$db"
: ...
echo foo >&3
out=$(mktemp /tmp/sonar.XXXXXX)
exec 3>"$out"
rm "$out"
: ...
echo foo >&3
$DIR/sonar.py add --db "$db" -f "$DIR/test/test.fasta"
$DIR/sonar.py update --db "$db" --pangolin "$DIR/test/test_pangolin.csv"
$DIR/sonar.py update --db "$db" --tsv "$DIR/test/test.tsv.gz" --fields accession=accessions zip=regions date=dates gisaid=gisaid ena=ena lab=lab source=source collection=collection technology=technology platform=platform chemistry=chemistry material=material ct=ct software=software version=software_version
$DIR/sonar.py optimize --db "$db"

$DIR/sonar.py match --db "$db" | (sed -u 1q; sort) > "$out"
if diff -q $DIR/test/expected.csv "$out"; then
    echo 'match 1 as expected'
else
    echo 'match 1 not as expected'
		diff $DIR/test/expected.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --date 2021-01-01:2021-01-31 -i C241T --db "$db" > "$out"
if diff -q $DIR/test/expected2.csv "$out"; then
    echo 'match 2 as expected'
else
    echo 'match 2 not as expected'
		diff $DIR/test/expected2.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match -i G1820A --db "$db" > "$out"
if diff -q $DIR/test/expected2.csv "$out"; then
    echo 'match 3 as expected'
else
    echo 'match 3 not as expected'
		diff $DIR/test/expected2.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match -e C1348T -i G1820A --db "$db" > "$out"
if diff -q $DIR/test/expected2.csv "$out"; then
    echo 'match 4 as expected'
else
    echo 'match 4 not as expected'
		diff $DIR/test/expected2.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --acc test2 --db "$db" > "$out"
if diff -q $DIR/test/expected2.csv "$out"; then
    echo 'match 5 as expected'
else
    echo 'match 5 not as expected'
		diff $DIR/test/expected2.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --zip 0177 --db "$db" > "$out"
if diff -q $DIR/test/expected2.csv "$out"; then
    echo 'match 6 as expected'
else
    echo 'match 6 not as expected'
		diff $DIR/test/expected2.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --acc '^NC_045512.2' --lineage '^B.1.1.297' --db "$db" > "$out"
if diff -q $DIR/test/expected2.csv "$out"; then
    echo 'match 7 as expected'
else
    echo 'match 7 not as expected'
		diff $DIR/test/expected2.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

n=$($DIR/sonar.py match --count --lab l3 --source sentinel --db "$db")
if [ "$n" -eq 0 ]; then
    echo 'match 8 as expected'
else
    echo "match 8 not as expected (expected 0 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

echo "updating database ..."
$DIR/sonar.py add --db "$db" -f "$DIR/test/test2.fasta.xz" --quiet
$DIR/sonar.py update --db "$db" --tsv "$DIR/test/test.tsv.gz" --fields accession=accessions zip=regions date=dates gisaid=gisaid ena=ena lab=lab source=source collection=collection technology=technology platform=platform chemistry=chemistry material=material ct=ct software=software version=software_version

n=$($DIR/sonar.py match --count --lab l3 --source sentinel --db "$db")
if [ "$n" -eq 1 ]; then
    echo 'match 9 as expected'
else
    echo "match 9 not as expected (expected 1 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

n=$($DIR/sonar.py match --count --collection random --technology Illumina --db "$db")
if [ "$n" -eq 2 ]; then
    echo 'match 10 as expected'
else
    echo "match 10 not as expected (expected 2 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

n=$($DIR/sonar.py match --count --zip 033 --lab l3 --source sentinel --collection RANDOM --technology illumina --platform nextseq --chemistry flex cleanplex --material swap --min_ct 30 --max_ct 34 --software covpipe --version 3.0.5 --db "$db")
if [ "$n" -eq 1 ]; then
    echo 'match 11 as expected'
else
    echo "match 11 not as expected (expected 1 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

n=$($DIR/sonar.py match --count --min_ct 10 --max_ct 30 --db "$db")
if [ "$n" -eq 2 ]; then
    echo 'match 12 as expected'
else
    echo "match 12 not as expected (expected 2 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

n=$($DIR/sonar.py match --count --no_frameshifts --db "$db")
if [ "$n" -eq 3 ]; then
    echo 'match 13 as expected'
else
    echo "match 13 not as expected (expected 3 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

$DIR/sonar.py match --only_frameshifts --db "$db" --ambig > "$out"
if diff -q $DIR/test/expected3.csv "$out"; then
    echo 'match 14 as expected'
else
    echo 'match 14 not as expected'
		diff $DIR/test/expected3.csv "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

$DIR/sonar.py restore --db "$db" --acc test3 > "$out"
if diff -q $DIR/test/test2.fasta "$out"; then
    echo 'match 15 as expected'
else
    echo 'match 15 not as expected'
		diff $DIR/test/test2.fasta "$out"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

$DIR/sonar.py add --db "$db" -f "$DIR/test/test2.fasta.xz" --source TEST --noprogress

n=$($DIR/sonar.py match --count --source TEST --db "$db")
if [ "$n" -eq 1 ]; then
    echo 'match 16 as expected'
else
    echo "match 16 not as expected (expected 1 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

$DIR/sonar.py remove --db "$db" --acc "test3"

n=$($DIR/sonar.py match --count --source TEST --db "$db")
if [ "$n" -eq 0 ]; then
    echo 'match 17 as expected'
else
    echo "match 17 not as expected (expected 1 got $n)"
		conda deactivate
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

#cleaning
conda deactivate
conda env remove -y --name covsonar_ci
