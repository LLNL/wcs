#!/bin/bash

###############################################################################
#                            SSA integration tests
###############################################################################

# At the first round, this scrpipt runs ssa simulations and generates reference
# results. Subsequent runs compare the new results to the reference to see if
# they differ. If there is difference, that means the current version of the
# code differs from the previous version with which the reference results were
# generated.

# Set WCS_INSTALL_DIR to the path where the 'bin/ssa' executable can be found.
#
# These tests requires "WCS_WITH_EXPRKT=ON"
#
# The outline of this script is as follows
# - Define uility functions
# - Setup the path to reference results
#   - Some tests rely on runing simulations and comparing the results against
#     the past results of the same simulations.
# - Setup Global test parameters
#   - Defines two probems: one is simple and the other is complex
#   - Multiple random seeds to use
#   - SSA methods
# - Define SSA tests for correctness, sampling, and performance.
#   - Test correctnes by comparing a new tracing result with the reference one
#     and see if they differ. If so, there is likely an error.
#   - Test if sampling results are consistent with the tracing results.
#   - Test performance by measuring the wall clock time multiple times and
#     compare the average against the past performance.


###############################################################################
#                          Define utility functions
###############################################################################

# Check if the name of a variable that stores a path is defined. If so,
# check if the path stored in the variable exists.
# The argument of this function is the name of the variable that is
# supposed to store an accessible directory.

function check_dir () {
    if [ $# -ne 1 ] ; then
        echo "function ${0}() requires one arguments:" 1>&2
        echo "  1. variable_name_for_dir: name of the variable of a dir" 1>&2
        exit 1
    fi

    if [ -z ${!1+x} ] ; then
        echo -e "${1} is not defined!" 1>&2
        exit 1
    elif [ -z "${!1}" ] ; then
        echo -e "${1} is current working directory!" 1>&2
        eval ${1}='.'
    elif [ ! -d "${!1}" ] ; then
        echo -e "${1} '${!1}' does not exist!" 1>&2
        exit 1
    fi
}

# Convert a relative path to an absolute path.
# The target path mut be accessible for this function to work.

function absolute_path {
    if [ $# -ne 1 ] ; then
        echo "function ${0}() requires one arguments:" 1>&2
        echo "  1. path" 1>&2
        exit 1
    fi

    local ap_rel_path="$1"
    local ap_bname="$(basename "${ap_rel_path}")"

    if [ "${ap_rel_path}" == "." ]; then
        echo "$(pwd)"
    elif [ "${ap_rel_path}" == ".." ]; then
        echo "$(dirname "$(pwd)")"
    else
        if [ -z ${ap_rel_path} ] || [ ! -d ${ap_rel_path} ] ; then
            echo "${ap_rel_path} does not exist!" 1>&2
            exit 1
        else
            local ap_cwd="$(pushd "${ap_rel_path}" > /dev/null; pwd)"
            case ${ap_bname} in
                ..) echo "$(absolute_path "${ap_cwd}")";;
                 .) echo "${ap_cwd}";;
                 *) echo "${ap_cwd}";;
            esac
        fi
    fi
}

# Compare the simulation results of test runs against the reference results
function compare_output () {
    if [ $# -lt 2 ] || [ $# -gt 3 ] ; then
        echo "function ${0}() requires 2 or 3 arguments:" 1>&2
        echo " 1. my_dir: path to the new results to check" 1>&2
        echo " 2. ref_dir: path to the reference results to compare against" 1>&2
        echo " 3. [prob_type]: Optional argument that indentifies which type of" 1>&2
        echo "                 problem to use. E.g., prob_simple or prob_complex" 1>&2
        echo "                 When omitted, both are used." 1>&2
        exit 1
    fi

    local co_my_dir=$1
    local co_ref_dir=`absolute_path $2`
    local co_more_files=""
    local co_diff=0

    if [ $# -eq 3 ] ; then
        if [ -z ${!3+x} ] ; then
            echo "${3} is not defined!" 1>&2
            exit 1
        elif [ -z "${!3}" ] ; then
            echo "${3} is empty!" 1>&2
            exit 1
        fi
        local co_probs_to_compare="${!3}"
    else
        local co_probs_to_compare="${problems}"
    fi

    check_dir co_my_dir
    check_dir co_ref_dir

    for p in ${co_probs_to_compare}
    do
        local co_prob=`basename ${p} | sed 's/\.graphml//'`
        local co_prob_list="${co_prob_list}${co_prob} "
    done

    pushd ${co_my_dir} > /dev/null

    local co_files_to_check=$(echo ${co_prob_list} \
          | awk '{for (i=1; i<=NF;i++) {printf("%s*.txt ", $i);}}')
    local co_fl=`ls -1 ${co_files_to_check} 2> /dev/null`

    for co_f in ${co_fl}
    do
        if [ ! -f ${co_f} ] ; then
            echo "${co_f} does not exists!" 1>&2
            exit 1
        elif [ ! -f ${co_ref_dir}/${co_f} ] ; then
            echo "${co_ref_dir}/${co_f} does not exists!" 1>&2
            co_more_files="${co_more_files} ${co_f}"
        else
            local co_cnt=`diff ${co_f} ${co_ref_dir}/${co_f} | wc -l`
            if [ ${co_cnt} -gt 0 ] ; then
                echo "Inconsistent result in ${co_f}" 1>&2
                co_diff=1
            fi
        fi
    done

    # Detect the files that do not exist in the set of referece results
    if [ "${co_more_files}" != "" ] ; then
        co_more_files="Unchecked extra files: ${co_more_files}."
        echo "${co_more_files}" 1>&2
    fi

    if [ ${co_diff} -eq 0 ] ; then
        echo "OK"
    fi
    popd > /dev/null
}

function verify_samples () {
    if [ $# -ne 2 ] ; then
        echo "function ${0}() requires two arguments:" 1>&2
        echo " 1. sample_file: simulation sampling result file" 1>&2
        echo " 2. trace_file: tracing result file of the same simulation" 1>&2
    fi

    local vs_sfile=${1}
    local vs_tfile=${2}

    if [ ! -f ${vs_sfile} ] ; then
        echo "Sampling file '${vs_sfile}' does not exist!" 1>&2
        exit 1
    fi

    if [ ! -f ${vs_tfile} ] ; then
        echo "Tracing file '${vs_tfile}' does not exist!" 1>&2
        exit 1
    fi

    local vs_ok=$(awk \
       'BEGIN { ok = 1; }
        {
            if (FNR == 1) {
                fcnt ++
                if (($1 == "num_species") && ($4 == "num_reactions")) {
                    num_species = $3;
                    num_reactions = $6;
                }
            } else if (FNR > 2) {
                if (fcnt == 1) {
                    tid[$1] = FNR - 2;
                    for (i = 1; i <= num_species; i++) {
                        species[$1" "i] = $(i+1);
                    }
                    r = 0;
                    for (i = num_species + 2; i <= NF; i++) {
                        r = r + 1;
                        reactions[$1" "r] = $(i);
                    }
                    num_samples = num_samples + 1;
                } else {
                    if ($1 in tid) {
                        t = $1;
                        for (i = 1; i <= num_species; i++) {
                            if (species[t" "i] != $(i+1)) {
                                err_sample = tid[t];
                                ok = 0;
                                exit;
                            }
                        }
                        r = 0;
                        for (i = num_species + 3; i <= NF; i++) {
                            r = r + 1;
                            if (reactions[t" "r] != $(i)) {
                                err_sample = tid[t];
                                ok = 0;
                                exit;
                            }
                        }
                    }
                }
            }
        }

        END {
            if (ok) {
                printf("OK\n");
            } else {
                printf("NOT OK\n");
                printf("%d-th sample is incorrect\n", err_sample) > "/dev/stderr";
            }
        }' "${vs_sfile}" "${vs_tfile}")
    echo ${vs_ok}
}

# Compare the wall clock timing results recorded in a log file against those
# in a reference log file, and see if the average timing difference is within
# a given threshold.

function compare_performance_log ()
{
    if [ $# -ne 6 ] ; then
        echo "function ${0}() requires two arguments:" 1>&2
        echo " 1. log1: the current performance log file" 1>&2
        echo " 2. log2: the the reference log file" 1>&2
        echo " 3. threshold: permitted margin for average wallclock increase" 1>&2
        echo " 4. number of measurements: number of measurements in each log" 1>&2
        echo "    It is used to compute the average difference." 1>&2
        echo " 5. prob: the problem used to run tests. This is for annotating" 1>&2
        echo "          comparison." 1>&2
        echo " 6. method: the SSA method used to run tests. This is for" 1>&2
        echo "            annotating comparison." 1>&2
        exit 1
    fi

    local pl_logf1=${1}
    local pl_logf2=${2}
    local pl_threshold=${3}
    local pl_n_measure=${4}
    local pl_prob="${5}"
    local pl_method="${6}"

    if [ ! -f ${pl_logf1} ] ; then
        echo "${pl_logf1} does not exist!" 1>&2
        exit 1
    fi

    if [ ! -f ${pl_logf2} ] ; then
        echo "${pl_logf2} does not exist!" 1>&2
        exit 1
    fi

        diff ${pl_logf1} ${pl_logf2} \
    | awk -v threshold=${pl_threshold} -v nm=${pl_n_measure} \
          -v prob="${pl_prob}" -v method="${pl_method}" \
    'BEGIN { ok = 1; }
    {
        if ($1 == "---") {
            next;
        } else if (($2 != "Wall") || ($3 != "clock")) {
            for (i=1; i <= NF; i++) {
                if (index($i, prob) != 0) {
                    $1 = "";
                    printf("Problem may differ: %s\n", $0) > "/dev/stderr";
                    next;
                }
            }
            if ($1 == "<") {
                $1 = "";
                printf("Incorrect simulation result:\n%s\n", $0) > "/dev/stderr";
                ok = 0;
                exit
            }
        } else {
            if (($1 == "<") && ($(NF) == "(sec)")) {
                n++;
                t_new[n] = $(NF-1);
            } else if (($1 == ">") && ($(NF) == "(sec)")) {
                t_ref[n] = $(NF-1);
            }
        }
    }
    END {
        if (!ok) exit;
        for (i=1; i <= n; i++) {
            diff[i] = double (t_new[i] - t_ref[i])/t_ref[i];
            #printf("%u-th new time (%f) differs the reference (%f) by %.2f%\n",
            #       i, t_new[i], t_ref[i], 100*diff[i]) > "/dev/stderr";
            #if (diff[i] > threshold) {
            #    printf("exceeding the threshold %.2f\n", threshold) > "/dev/stderr";
            #}
            sum += diff[i];
        }
        if (sum >  threshold*nm) {
            printf("Average difference (%.2f) exceeds the threshold %.2f for method %u\n",
                   sum/nm, threshold, method) > "/dev/stderr";
            ok = 0;
        } else {
            printf("Average difference (%.2f) is within the threshold %.2f for method %u\n",
                   sum/nm, threshold, method) > "/dev/stderr";
        }
        printf("%s\n", ok? "OK" : "NOT OK");
    }'
}

###############################################################################
#                       Setup common test environment
###############################################################################

if [ -z "${WCS_INSTALL_DIR}" ] ; then
    WCS_INSTALL_DIR="$(absolute_path '../../install')"
fi
if [ -z "${WCS_SRC_DIR}" ] ; then
    WCS_SRC_DIR="$(absolute_path '../..')"
fi

WCS_TEST_DIR=${WCS_SRC_DIR}/tests
WCS_REF_OUTPUT_DIR=${WCS_TEST_DIR}/integration/expected_output

echo "========================================================================"
echo "WCS_INSTALL_DIR=${WCS_INSTALL_DIR}"
echo "WCS_SRC_DIR=${WCS_SRC_DIR}"
echo "WCS_TEST_DIR=${WCS_TEST_DIR}"
echo "WCS_REF_OUTPUT_DIR=${WCS_REF_OUTPUT_DIR}"
echo "========================================================================"


###############################################################################
#                       Setup common test parameters
###############################################################################

seeds="47 147 1147"
methods="0 1 2"
frag_sz=0

prob_simple="${WCS_TEST_DIR}/problem/Gillespie/eq29-exprtk.graphml"
prob_complex="${WCS_TEST_DIR}/problem/Birtwhistle/birtwistle.rate.graphml"
problems="${prob_simple} ${prob_complex}"


###############################################################################
#                       Perform sanity Checks
###############################################################################

function sanity_checks () {
    check_dir WCS_SRC_DIR
    check_dir WCS_INSTALL_DIR
    check_dir WCS_TEST_DIR

    if [ "${GENERATE_REF}" != "1" ] ; then
        check_dir WCS_REF_OUTPUT_DIR
    fi

    for p in ${problems} ; do
        if [ ! -f $p ] ; then
            echo -e "Problem input file '$p' does not exist!"
            exit 1
        fi
    done

   local  compile_check=$(grep '#define WCS_HAS_EXPRTK' ${WCS_INSTALL_DIR}/include/wcs_config.hpp)
    if [ "${compile_check}" == "" ] ; then
        echo "These tests requires the 'WCS_WITH_EXPRTK=ON' cmake option set because"
        echo "the reference results are produced with the binary compiled with it."
        exit 1
    fi
}

###############################################################################
#                       Define SSA correctnes test
###############################################################################

# Generate tracing results using both simple and complex problems for a certain
# number of iterations with different seeds and different SSA methods.
# Compare the results against the reference ones to see if they are identical.

function ssa_correctness () {
    local sc_tname=${FUNCNAME[0]}
    local sc_exec="${WCS_INSTALL_DIR}/bin/ssa"
    local sc_n_step=20  # number of simulation steps

    mkdir -p ${sc_tname}

    echo "================ ${sc_tname} ..." 1>&2

    for sc_p in ${problems}
    do
        local sc_prob=`basename ${sc_p} | sed 's/\.graphml//'`

        for sc_m in ${methods}
        do
            for sc_s in ${seeds}
            do
                local sc_of=${sc_tname}/${sc_prob}.m${sc_m}.s${sc_s}.i${sc_n_step}.txt
                local sc_arg="-d -i ${sc_n_step} -m ${sc_m} -s ${sc_s} -o ${sc_of} -f ${frag_sz}"
                echo "ssa ${sc_arg} ${sc_prob}.graphml"
                eval "${sc_exec} ${sc_arg} ${sc_p}"
                echo ""
            done
        done
    done > ${sc_tname}/log.txt

    if [ "${GENERATE_REF}" == "1" ] ; then
        echo "OK"
    else
        echo `compare_output ${sc_tname} ${WCS_REF_OUTPUT_DIR}/${sc_tname}`
    fi
}


###############################################################################
#                       Define SSA sampling test
###############################################################################

# Check if the sampling produces correct results.
# Here we only check if the new sampling results are identical to the
# reference sampling results beucase the consistency of the reference
# sampling results is pre-verified against tracing by the other function
# 'check_sample_file()' in advance.
# This test only uses the simple problem.

function ssa_sampling () {
    local ss_tname=${FUNCNAME[0]}
    local ss_exec="${WCS_INSTALL_DIR}/bin/ssa"
    mkdir -p ${ss_tname}
    local ss_n_step=20  # number of simulation steps
    local ss_ismpl='i3' # cmd line option for sampling by iteration
    local ss_tsmpl='t0.0002' # cmd line option for sampling by time

    echo "================ ${ss_tname} ..." 1>&2

    for ss_p in ${prob_simple}
    do
        local ss_prob=`basename ${ss_p} | sed 's/\.graphml//'`

        for ss_m in ${methods}
        do
            for ss_s in `echo ${seeds} | cut -d ' ' -f 1`
            do
                local ss_of  # outfile
                local ss_arg
                ss_of=${ss_tname}/${ss_prob}.m${ss_m}.s${ss_s}.i${ss_n_step}.r${ss_ismpl}.txt
                ss_arg="-i ${ss_n_step} -m ${ss_m} -s ${ss_s} -o ${ss_of} -f ${frag_sz} -r ${ss_ismpl}"
                echo "ssa ${ss_arg} ${ss_prob}.graphml"
                eval "${ss_exec} ${ss_arg} ${ss_p}" 2>&1
                echo ""
                ss_of=${ss_tname}/${ss_prob}.m${ss_m}.s${ss_s}.i${ss_n_step}.r${ss_tsmpl}.txt
                ss_arg="-i ${ss_n_step} -m ${ss_m} -s ${ss_s} -o ${ss_of} -f ${frag_sz} -r ${ss_tsmpl}"
                echo "ssa ${ss_arg} ${ss_prob}.graphml"
                eval "${ss_exec} ${ss_arg} ${ss_p}" 2>&1
                echo ""
            done
        done
    done > ${ss_tname}/log.txt

    if [ "${GENERATE_REF}" == "1" ] ; then
        echo "OK"
    else
        echo `compare_output ${ss_tname} ${WCS_REF_OUTPUT_DIR}/${ss_tname} prob_simple`
    fi
}

# Check if a file that constains a sequence of samples of the simulation state
# is consistent with the tracing results. This is done by comparing every sample,
# for example, taken at time t, has the same state as that recorded at time t by
# tracing for the same simulation.

function check_sample_file () {
    if [ $# -ne 4 ] && [ $# -ne 7 ] ; then
        echo "function ${0}() requires 5 or 7 arguments:" 1>&2
        echo " 1. SSA method: e.g., 0, 1, or 2." 1>&2
        echo " 2. seed: the random number seed used." 1>&2
        echo " 3. n_iter: the number of simulation steps." 1>&2
        echo " 4. interval: the sampling interval used. e.g., i3 or t0.0002." 1>&2
        echo " 5. trace_dir: the directory that contains tracing results" 1>&2
        echo "               such as 'ssa_correctness'." 1>&2
        echo " 6. sample_dir: the directory that contains sampling results" 1>&2
        echo "                such as 'ssa_sampling'." 1>&2
        echo " 7. problem: the input model" 1>&2
        exit 1
    fi

    local cs_m=${1}     # SSA method used
    local cs_s=${2}     # seed used
    local cs_n=${3}     # number of simulation steps ran
    local cs_i=${4}     # sampling interval option used
    local cs_tdir="ssa_correctness" # tracing result dir
    local cs_sdir="ssa_sampling"   # sampling result dir
    local cs_p=$(basename "${prob_simple}" | sed 's/\.graphml//')

    if [ $# -eq 7 ] ; then
        cs_tdir=${5}
        cs_sdir=${6}
        cs_p=$(basename "${7}" | sed 's/\.graphml//')
    fi

    local cs_tfile=${cs_tdir}/${cs_p}.m${cs_m}.s${cs_s}.i${cs_n}.txt
    local cs_sfile=${cs_sdir}/${cs_p}.m${cs_m}.s${cs_s}.i${cs_n}.r${cs_i}.txt

    echo "$(verify_samples "${cs_sfile}" "${cs_tfile}")"
}

# Check if the sampling produces correct results.
# This is similar to 'ssa_sampling()' above. However, this does not require
# reference sampling results. Instead, it makes direct comparison to the
# tracing results, which can be slower as the trace files are larger in
# general. So, this function can be used once at the first time. Then,
# for subsequent verifications, 'ssa_sampling()' can be used instead.

function ssa_sampling_0 () {
    local s0_tname=${FUNCNAME[0]}
    local s0_tdir="ssa_correctness" # tracing result dir
    local s0_sdir="ssa_sampling"    # sampling result dir

    local s0_n_step=20  # number of simulation steps
    local s0_ismpl='i3' # cmd line option for sampling by iteration
    local s0_tsmpl='t0.0002' # cmd line option for sampling by time
    local s0_ok="OK"

    if [ "${GENERATE_REF}" == "1" ] ; then
        exit
    fi
    echo "================ ${s0_tname} ..." 1>&2

    for s0_p in ${prob_simple}
    do
        for s0_m in ${methods}
        do
            for s0_s in `echo ${seeds} | cut -d ' ' -f 1`
            do
                s0_ok=$(check_sample_file ${s0_m} ${s0_s} ${s0_n_step} \
                          ${s0_ismpl} ${s0_tdir} ${s0_sdir} ${s0_p})
                if [ "${s0_ok}" != "OK" ] ; then
                    break 3
                fi
                s0_ok=$(check_sample_file ${s0_m} ${s0_s} ${s0_n_step} \
                          ${s0_tsmpl} ${s0_tdir} ${s0_sdir} ${s0_p})
                if [ "${s0_ok}" != "OK" ] ; then
                    break 3
                fi
            done
        done
    done
    echo "${s0_ok}"
}


###############################################################################
#                       Define SSA performance test
###############################################################################

# Measure the average wall clock times to run a simulation, and see if it is
# comparable to the past performance. i.e., Check if the increase in time is
# within a permitted margin in terms of the fraction of time.
# This test only uses the complex problem. The first argument is the threshold
# and the second argument is the number of measurements for each seed and each
# SSA method. Finally, all the measurements are averaged per method for
# comparison against the reference performance.

function ssa_performance () {
    if [ $# -ne 1 ] && [ $# -ne 2 ] ; then
        echo "function ${0}() requires one or two arguments:" 1>&2
        echo " 1. threshold: permitted margin for average wallclock increase" 1>&2
        echo " 2. number of measurements: repeat the measurement as many as" 1>&2
        echo "    this per combination of seed and SSA method." 1>&2
        echo "    Finally, compute the average timing for each SSA method" 1>&2
        echo "    and compare it against the reference timing to see if it" 1>&2
        echo "    is within the threshold." 1>&2
        exit 1
    fi

    local sp_tname=${FUNCNAME[0]}
    local sp_exec="${WCS_INSTALL_DIR}/bin/ssa"
    mkdir -p ${sp_tname}

    local sp_n_step=5000    # number of simulation steps
    local sp_threshold=${1} # permitted upper margin of performance variation
    local sp_n_repeat=1     # number of runs of the same simulation.
                            #i.e., with the same method and the same seed

    if [ $# -eq 2 ] ; then
        sp_n_repeat=${2}
    fi

    echo "================ ${sp_tname} ..." 1>&2

    for sp_m in ${methods}
    do
        local sp_logf=${sp_tname}/log.m${sp_m}.i${sp_n_step}.txt

        local sp_p=${prob_complex}
        local sp_prob=`basename ${sp_p} | sed 's/\.graphml//'`
        local sp_n_measure=$(echo ${seeds} | awk -v nr=${sp_n_repeat} '{ print NF*nr; }')

        for sp_c in `seq 1 ${sp_n_repeat}`
        do
            for sp_s in ${seeds}
            do
                local sp_arg="-i ${sp_n_step} -m ${sp_m} -s ${sp_s}"
                echo "ssa ${sp_arg} ${sp_prob}.graphml"
                eval "${sp_exec} ${sp_arg} ${sp_p}"
                echo ""
                local sp_of=${sp_tname}/${sp_prob}.m${sp_m}.s${sp_s}.i${sp_n_step}.c${sp_c}.tr.txt
                local sp_arg_tr="${sp_arg} -f ${frag_sz} -d -o ${sp_of}"
                echo "ssa ${sp_arg_tr} ${sp_prob}.graphml"
                eval "${sp_exec} ${sp_arg_tr} ${sp_p}"
                echo ""
            done
        done > ${sp_logf}

        if [ "${GENERATE_REF}" == "1" ] ; then
            continue
        fi
        if [ ! -f ${WCS_REF_OUTPUT_DIR}/${sp_logf} ] ; then
            echo "${WCS_REF_OUTPUT_DIR}/${sp_logf} does not exist!" 1>&2
            continue
        fi

        compare_performance_log ${sp_logf} ${WCS_REF_OUTPUT_DIR}/${sp_logf} \
                  ${sp_threshold} ${sp_n_measure} "${prob_complex}" ${sp_m}
    done
}


###############################################################################
#                             Run SSA tests
###############################################################################

if [ ! -d expected_output ] ; then
    if [ -f expected_output.tgz ] ; then
        tar zxvf expected_output.tgz
    else
        GENERATE_REF=1
    fi
fi

sanity_checks

ssa_correctness
ssa_sampling
#ssa_sampling_0
ssa_performance 0.05 4

if [ "${GENERATE_REF}" == "1" ] ; then
    mkdir -p expected_output
    mv ssa_correctness expected_output
    mv ssa_sampling expected_output
    mv ssa_performance expected_output
    tar zcvf expected_output.tgz expected_output
fi
