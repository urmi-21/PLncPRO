#!/bin/sh

testBuild() {
	python build.py -p sample_data/train/mrna.fa -n sample_data/train/lncrna.fa -o tests/out/build_out -m sample_model -d /home/usingh/work/urmi/hoap/pldb/mydb -t 5 -v -r
}

testPred() {
    assertEquals 1 2
}

. shunit2-2.1.6/src/shunit2