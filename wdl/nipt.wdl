## Copyright NGeneBio, 2018
## This WDL 
## Requrements/expectations :
## Outputs :
## Cromwell version support
## LICENSING :

# WORKFLOW DEFINITION
workflow NIPTWf {
   File sample_bam
   String bedtools
   File? hg19_chrom_size
   File? hg19_bins
   File? hg19_gc
   File? coverage_bins
   String? window_size_override
   Boolean? createbin
   Boolean? createcov
   Boolean createbin_or_default = select_first([createbin, false])
   Boolean createcov_or_default = select_first([createcov, false])
   String default_window_size = select_first([window_size_override, "100000"])
   File nipt_python
  
   # build bin size
   if (createbin_or_default) {
      call CreateBin {
         input:
            hg19_chrom_size = hg19_chrom_size,
            window_size = default_window_size,
            bedtools = bedtools
      }      
   }

   if (createcov_or_default){
      call BuildCoverage {
         input:
            sample_bam = sample_bam,
            hg19_bins = CreateBin.out_bins,
            #hg19_bins = hg19_bins,
            bedtools = bedtools
      }      
   }

   call RunNipt {
      input:
          #bam_coverage_bins = coverage_bins,
          bam_coverage_bins = BuildCoverage.coverage_bins,
          hg19_gc = hg19_gc,
          nipt_python = nipt_python
   }
   
   output {
      File out_bins = "hg19_bins"
      File out_coverage_bins = "cov.bins"
      File nipt_png = RunNipt.nipt_png
      File nipt_png2 = RunNipt.nipt_png2
   } 
}

#TASK DEFINITIONS
task RunNipt {
   File bam_coverage_bins
   File nipt_python
   File hg19_gc
   command {
      source activate env
      cd /cromwell_root
      cp ${nipt_python} ${bam_coverage_bins} ${hg19_gc} /cromwell_root
      python nipt2.py
   }
   output {
      File nipt_png = "nipt.png"
      File nipt_png2 = "nipt_recal.png"
   }
   runtime {
      docker: "hongiiv/hongiiv-nipt:1"
      cpu: "1"
      memory: "5GB"
      preemptible: 3
   }
}

task BuildCoverage {
   File sample_bam
   File hg19_bins
   String bedtools
   command {
      ${bedtools} coverage -abam ${sample_bam} -b ${hg19_bins} > cov.bins
   }
   output {
      File coverage_bins = "cov.bins"   
   }
   runtime{
      docker: "hongiiv/hongiiv-nipt:1"
      cpu: "1"
      memory: "5GB"
      disks: "local-disk " + 100 + " HDD"
      preemptible: 3
   }
}
task CreateBin {
   # Command parameters
   File hg19_chrom_size
   String window_size
   String bedtools

   # Runtime paramters
      
   command {
      ${bedtools} makewindows -g ${hg19_chrom_size} -w ${window_size} > hg19_bins 
   }
   output {
      File out_bins = "hg19_bins"
   }
   runtime {
      docker: "hongiiv/hongiiv-nipt:1"
      cpu: "1"
      memory: "1GB"
      preemptible: 3
   }
}
