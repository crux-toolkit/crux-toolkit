#!/bin/bash
#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=32.0G
#$ -l h_rt=128:0:0
#$ -q noble-long.q

hostname
date
echo PID=$$

if [[ -e /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
  module load modules modules-init modules-gs modules-noble
fi

# These three options make it harder to for some part of a script to
# fail without you recognizing it. nounset means that references to an
# unset variable result in an error. This means that you can no longer
# do stuff like this if VAR is potentially unset, because "$VAR" returns
# an error rather than "":
#
# if [ "$VAR" ]; then
#
# fi
#
# To explicitly indicate that you are OK with the variable potentially
# being empty you can instead use ${VAR:-}.

# To avoid errexit for a single command, use "|| true", e.g.,
#    diff foo foobar || true

set -o nounset
set -o pipefail
set -o errexit 
set -o xtrace

# Location of the data
ms2_file=../performance-tests/051708-worm-ASMS-10.ms2
fasta_file=../performance-tests/worm+contaminants.fa

#CRUX=../../src/crux
CRUX=../../build.release/src/crux

# Avoid NFS overhead by using a scratch directory.
scratch_dir=/scratch
if [[ ! -e $scratch_dir ]]; then
    scratch_dir=.
else
    cp $ms2_file $scratch_dir
    ms2_file=`basename $ms2_file`
    ms2_file="$scratch_dir/$ms2_file"
    cp $fasta_file $scratch_dir
    fasta_file=`basename $fasta_file`
    fasta_file="$scratch_dir/$fasta_file"
fi

# Build the index.
index=$scratch_dir/my_index
if [[ ! -e $index ]]; then
    $CRUX tide-index --decoy-format none \
	  --output-dir $index \
	  $fasta_file $index
fi

# Convert the MS2 to spectrumrecords
spectrum_records=$scratch_dir/my_spectra
if [[ ! -e $spectrum_records ]]; then
    $CRUX tide-search \
	  --output-dir tmp \
	  --store-spectra $spectrum_records \
	  $ms2_file $index
    rm -r tmp
fi

html=results.html
echo "<html><body><pre>" > $html
for precursor in 3 10; do

    comet_params="--decoy_search 0 --num_output_lines 1"
    tide_params="--top-match 1"

    if [[ $precursor == 10 ]]; then
	comet_precursor="--peptide_mass_units 2 --peptide_mass_tolerance 10"
	tide_precursor="--precursor-window 10 --precursor-window-type ppm"
    else
	comet_precursor=""
	tide_precursor=""
    fi
    
    for fragment in 1 02; do
	for threads in 1 4; do
	    comet_threads="--num_threads $threads"
	    tide_threads="--num-threads $threads"
	    
	    for engine in tide-xcorr tide-p tide-combined comet; do
		root=$engine.pre=$precursor.frag$fragment.threads$threads

                if [[ $fragment == 02 ]]; then
                    comet_fragment="--fragment_bin_tol 0.02"
                    if [[ $engine == "tide-p" ]]; then
                        continue # Can't do high-res p-values.
                    elif [[ $engine == "tide-combined" ]]; then
	                tide_fragment="--mz-bin-width 1.0005079 --fragment-tolerance 0.02"
                    else
	                tide_fragment="--mz-bin-width 0.02"
                    fi
	        else
	            comet_fragment="--fragment_bin_tol 1.0005079"
	            tide_fragment="--mz-bin-width 1.0005079 --fragment-tolerance 1.0005079"
	        fi
	
		# Select among the four different search engines.
		log_file=$scratch_dir/$root/tide-search.log.txt
		if [[ $engine == "tide-xcorr" ]]; then
		    search_command="tide-search"
		elif [[ $engine == "tide-p" ]]; then
		    search_command="tide-search --exact-p-value T"
		elif [[ $engine == "tide-combined" ]]; then
		    search_command="tide-search --exact-p-value T --score-function both"
		elif [[ $engine == "comet" ]]; then
		    log_file=$scratch_dir/$root/comet.log.txt
		    search_command="comet"
		fi

		# Run the actual search.
		if [[ ! -e $log_file ]]; then
		    if [[ $engine == "comet" ]] ; then
  			$CRUX $search_command $comet_params \
			      $comet_precursor \
			      $comet_fragment \
			      $comet_threads \
			      --output-dir $scratch_dir/$root --overwrite T \
			      $ms2_file $fasta_file
		    else
  			$CRUX $search_command $tide_params \
			      $tide_precursor \
			      $tide_fragment \
			      $tide_threads \
			      --output-dir $scratch_dir/$root --overwrite T \
			      $spectrum_records $index
		    fi
		fi
		echo -n "$root " >> $html
		awk -F ":" '$2 == " Elapsed time" {print $3}' $log_file >> $html
	    done
	done
    done
done
echo "</pre></body></html>" >> $html		
