#!/bin/bash
load $HOME/bin/bats-support/load.bash
load $HOME/bin/bats-assert/load.bash
load $HOME/bin/bats-file/load.bash

# create some variables for tests to use later
setup_file() {

    # a timestamp directory to contain test results
    datetime=$(date --iso-8601=seconds)
    test_dir=$PWD
    test_dir_mod=$PWD/tests/FastQC
    test_time_dir=$test_dir_mod/$datetime
    mkdir $test_time_dir
    test_input=$test_dir_mod/Inputs/
    expected_outdir=$test_dir_mod/Outputs
    
    # export these variables for use in other tests
    export test_time_dir
    export test_dir
    export test_dir_mod
    export expected_outdir
    export test_input
}


setup() {
    mkdir $test_time_dir/${BATS_TEST_NUMBER}_$BATS_TEST_DESCRIPTION
    mkdir $test_time_dir/${BATS_TEST_NUMBER}_$BATS_TEST_DESCRIPTION/results
    outdir=$test_time_dir/${BATS_TEST_NUMBER}_$BATS_TEST_DESCRIPTION/results
    export outdir
}

# Run nextflow version
@test "RunNextflowVersion" {
    run $HOME/nextflow -version 
    [ "$status" -eq 0 ]
    echo "$output" >$outdir/RunNextflowVersion.txt
    assert_file_not_empty $outdir/RunNextflowVersion.txt
}

@test "TestFastQCModule" {
    run $HOME/nextflow run $test_dir/modules/FastQC/fastqc.nf -c $test_dir/configs/FastQC/fastqc.config --qcInputs=$test_input --verOutputs=$outdir --qcOutputs=$outdir
    [ "$status" -eq 0 ]
    assert_exists $outdir/
    for dir in $outdir/* 
    do
        for file in $dir
        do
            assert_file_exists "$file"
            assert_file_not_empty "$file"
            fullname=$(wc -c $file | cut -f2 --delimiter=' ')
            name=$(basename $fullname)
            bytesize=$(wc -c $file | cut -f1 --delimiter=' ' | cut -f10 --delimiter='/')
            assert_file_size_equals ${expected_outdir}/${name} $bytesize
        done
    done
}