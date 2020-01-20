#!/bin/sh

testBuild() {
	python build.py -p sample_data/train/mrna.fa -n sample_data/train/lncrna.fa -o tests/out/build_out -m sample_model -d tests/pldb/mydb -t 5 -v -r
}

testPred() {
    python prediction.py -i sample_data/test/test.fa -o tests/out/output_dir -p output_file_name -t 6 -d tests/pldb/mydb -m tests/out/build_out/sample_model -v -r
}

. shunit2-2.1.6/src/shunit2

