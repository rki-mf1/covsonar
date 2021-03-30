#!/bin/bash

# set bash mode to strict
set -euo pipefail
IFS=$'\n\t'

CWD=$(pwd)
DIR=$(cd "$(dirname "$0")"; pwd)

#conda env
conda env create --name covsonar_ci --force --file $DIR/sonar.env.yml
source $(conda info --root)/etc/profile.d/conda.sh
set +u
conda activate covsonar_ci
set -u

#doctest
cd $DIR/lib
#./sonardb.py
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
$DIR/sonar.py update --db "$db" --tsv "$DIR/test/test.tsv" --fields accession=accessions zip=regions date=dates gisaid=gisaid ena=ena
$DIR/sonar.py optimize --db "$db"
$DIR/sonar.py match --db "$db" | sort > "$out"
if cmp -s $DIR/test/expected.csv "$out"; then
    echo 'match 1 as expected'
else
    echo 'match 1 not as expected'
		diff $DIR/test/expected.csv "$out"
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --date 2021-01-01:2021-01-31 -i C241T --db "$db" > "$out"
if cmp -s $DIR/test/expected2.csv "$out"; then
    echo 'match 2 as expected'
else
    echo 'match 2 not as expected'
		diff $DIR/test/expected2.csv "$out"
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match -i G1820A --db "$db" > "$out"
if cmp -s $DIR/test/expected2.csv "$out"; then
    echo 'match 3 as expected'
else
    echo 'match 3 not as expected'
		diff $DIR/test/expected2.csv "$out"
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match -e C1348T -i G1820A --db "$db" > "$out"
if cmp -s $DIR/test/expected2.csv "$out"; then
    echo 'match 3 as expected'
else
    echo 'match 3 not as expected'
		diff $DIR/test/expected2.csv "$out"
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --acc test2 --db "$db" > "$out"
if cmp -s $DIR/test/expected2.csv "$out"; then
    echo 'match 4 as expected'
else
    echo 'match 4 not as expected'
		diff $DIR/test/expected2.csv "$out"
	  conda env remove -y --name covsonar_ci
	  exit 1
fi
$DIR/sonar.py match --zip 177 --db "$db" > "$out"
if cmp -s $DIR/test/expected2.csv "$out"; then
    echo 'match 5 as expected'
else
    echo 'match 5 not as expected'
		diff $DIR/test/expected2.csv "$out"
	  conda env remove -y --name covsonar_ci
	  exit 1
fi

#cleaning
conda deactivate
conda env remove -y --name covsonar_ci
